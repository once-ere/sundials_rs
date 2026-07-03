/* -----------------------------------------------------------------
 * Translated from src/sundials/sundials_hashmap.c and
 * src/sundials/sundials_hashmap_impl.h (SUNDIALS 7.7.0).
 *
 * A simple hashmap implementation for string keys, using linear
 * probing (without wrap-around, exactly like the C code) to resolve
 * collisions. The map owns its values and frees them on destruction.
 *
 * Design notes for the Rust translation:
 *  - The map is generic over the value type: `SUNHashMap<V>` (C used
 *    `void*` values plus a `destroyKeyValue` callback). A generic
 *    parameter was chosen over `Box<dyn Any>` because both intended
 *    clients store one homogeneous value type per map — the logger
 *    stores output streams (`Box<dyn Write>`) and the profiler stores
 *    timer records — and a generic map gives them typed access with
 *    no downcasting. The C `destroyKeyValue` callback argument is
 *    dropped: `Drop` on `V` performs the equivalent cleanup.
 *  - The C backing store `SUNStlVector_SUNHashMapKeyValue` (a vector
 *    of nullable pointers kept with size == capacity) becomes
 *    `Vec<Option<SUNHashMapKeyValue<V>>>`. `SUNHashMap_Capacity`
 *    therefore returns `buckets.len()` (the C vector is always
 *    pushed full of NULL entries so its size equals its capacity).
 *  - `SUNHashMap_Iterate`'s `void* ctx` argument is subsumed by a
 *    closure (`FnMut`) that captures its context.
 *  - Out-pointer arguments (`SUNHashMap*`, `void**`) become
 *    `&mut Option<...>` out-parameters so the exact C names and
 *    return-flag conventions are preserved.
 *  - `SUNHashMap_Values` drops the C `value_size` byte-count argument
 *    (it only fed `malloc`) and `SUNHashMap_Sort`'s qsort comparator
 *    becomes a `fn(..) -> Ordering`. Rust's `sort_by` is stable while
 *    qsort's order for equal elements is unspecified, so a stable
 *    order is a valid refinement.
 *  - NULL-argument checks on pointers that are references in Rust
 *    (`map == NULL`, `key == NULL`, ...) drop out; the corresponding
 *    error returns are unreachable here.
 *  - In `SUNHashMap_Remove` the C code dereferences the probed bucket
 *    unconditionally (undefined behavior if the key is absent); the
 *    Rust version returns SUNHASHMAP_KEYNOTFOUND instead.
 * -----------------------------------------------------------------*/

use std::cmp::Ordering;
use std::io::Write;

use crate::sundials_errors::*;

pub const SUNHASHMAP_ERROR: i64 = -99;
pub const SUNHASHMAP_KEYNOTFOUND: i64 = -1;
pub const SUNHASHMAP_DUPLICATE: i64 = -2;

const HASH_PRIME: u64 = 14695981039346656037;
const HASH_OFFSET_BASIS: u64 = 1099511628211;

/* stl/sunstl_vector.h: growth factor applied when a resize is needed */
const SUNSTLVECTOR_GROWTH_FACTOR: f64 = 1.5;

/// struct SUNHashMapKeyValue_ (sundials_hashmap_impl.h)
pub struct SUNHashMapKeyValue<V> {
    pub key: String,
    pub value: V,
}

/// struct SUNHashMap_ (sundials_hashmap_impl.h). The `destroyKeyValue`
/// callback member is dropped (`Drop` handles cleanup).
pub struct SUNHashMap<V> {
    pub buckets: Vec<Option<SUNHashMapKeyValue<V>>>,
}

/*
  This is a 64-bit implementation of the 'a' modification of the
  Fowler-Noll-Vo hash (i.e., FNV1-a), matching the C source exactly
  (including its naming of the two constants). The C loop reads
  signed `char`s, so bytes >= 0x80 are sign-extended before the xor.
*/
fn fnv1a_hash(s: &str) -> u64 {
    let mut hash: u64 = HASH_OFFSET_BASIS;
    for b in s.bytes() {
        hash = (hash ^ ((b as i8) as i64 as u64)).wrapping_mul(HASH_PRIME);
    }
    hash
}

/* We want the index to be in [0, SUNHashMap_Capacity(map)) */
fn sunHashMapIdxFromKey<V>(map: &SUNHashMap<V>, key: &str) -> i64 {
    let end = SUNHashMap_Capacity(map) - 1;
    if end == 0 {
        end
    } else {
        (fnv1a_hash(key) % (end as u64)) as i64
    }
}

/*
  This function creates a new SUNHashMap object allocated to hold
  up to 'capacity' entries (C also took a `destroyKeyValue` callback,
  replaced by `Drop`).

  **Returns:**
    * A SUNErrCode indicating success or a failure
*/
pub fn SUNHashMap_New<V>(capacity: i64, map: &mut Option<SUNHashMap<V>>) -> SUNErrCode {
    if capacity <= 0 {
        return SUN_ERR_ARG_OUTOFRANGE;
    }

    *map = None;

    /* Initialize all buckets to NULL */
    let mut buckets: Vec<Option<SUNHashMapKeyValue<V>>> = Vec::with_capacity(capacity as usize);
    for _ in 0..capacity {
        buckets.push(None);
    }

    *map = Some(SUNHashMap { buckets });

    SUN_SUCCESS
}

/*
  This function returns the capacity of the hashmap.
*/
pub fn SUNHashMap_Capacity<V>(map: &SUNHashMap<V>) -> i64 {
    map.buckets.len() as i64
}

/*
  This function frees the SUNHashMap object.
*/
pub fn SUNHashMap_Destroy<V>(map: &mut Option<SUNHashMap<V>>) -> SUNErrCode {
    *map = None;
    SUN_SUCCESS
}

/*
  This function iterates the map over the range [start, N]. N is either the
  index at which ``yieldfn`` indicates the iteration should stop, or the max
  entries in the map.

  ``yieldfn`` should return SUNHASHMAP_ERROR to continue the iteration, or
  [0, SUNHASHMAP_KEYNOTFOUND] to stop; the C `void* ctx` argument is
  subsumed by the closure's captures.

  **Returns:**
    * ``SUNHASHMAP_ERROR`` -- an error occurred
    * ``capacity`` -- iterated the whole map
    * ``>=0`` -- the index at which the iteration stopped
*/
pub fn SUNHashMap_Iterate<V>(
    map: &SUNHashMap<V>,
    start: i64,
    mut yieldfn: impl FnMut(i64, &Option<SUNHashMapKeyValue<V>>) -> i64,
) -> i64 {
    let mut i = start;
    while i < map.buckets.len() as i64 {
        let retval = yieldfn(i, &map.buckets[i as usize]);
        if retval == SUNHASHMAP_ERROR {
            /* keep looking */
        } else {
            return retval; /* yieldfn indicates the loop should break */
        }
        i += 1;
    }

    SUNHashMap_Capacity(map)
}

fn sunHashMapLinearProbeInsert<V>(idx: i64, kv: &Option<SUNHashMapKeyValue<V>>) -> i64 {
    /* find the next open spot */
    if kv.is_none() {
        return idx; /* open spot found at idx */
    }
    SUNHASHMAP_ERROR /* keep looking */
}

fn sunHashMapResize<V>(map: &mut SUNHashMap<V>) -> SUNErrCode {
    let old_capacity = SUNHashMap_Capacity(map);
    let new_capacity = if old_capacity == 0 {
        2
    } else {
        ((old_capacity as f64) * SUNSTLVECTOR_GROWTH_FACTOR).ceil() as i64
    };

    /* Set all buckets of the new store to NULL */
    let mut new_buckets: Vec<Option<SUNHashMapKeyValue<V>>> =
        Vec::with_capacity(new_capacity as usize);
    for _ in 0..new_capacity {
        new_buckets.push(None);
    }
    let mut old_buckets = std::mem::replace(&mut map.buckets, new_buckets);

    /* Rehash and reinsert (from the highest old index down, like the C code) */
    for i in (0..old_capacity).rev() {
        if let Some(kvp) = old_buckets[i as usize].take() {
            /* the C code ignores the insert status here as well */
            let _ = SUNHashMap_Insert(map, &kvp.key, kvp.value);
        }
    }

    SUN_SUCCESS
}

/*
  This function creates a key-value pair and attempts to insert it into the map.
  Will use linear probing if there is a collision (the map keeps a copy of the
  key and takes ownership of the value).

  **Returns:**
    * ``0`` -- success
    * ``SUNHASHMAP_ERROR`` -- an error occurred
    * ``SUNHASHMAP_DUPLICATE`` -- duplicate key
*/
pub fn SUNHashMap_Insert<V>(map: &mut SUNHashMap<V>, key: &str, value: V) -> i64 {
    let mut idx = sunHashMapIdxFromKey(map, key);

    /* Check if the bucket is already filled (i.e., we might have had a collision) */
    if let Some(kvp) = &map.buckets[idx as usize] {
        /* Determine if key is actually a duplicate (not allowed) */
        if key == kvp.key {
            return SUNHASHMAP_DUPLICATE;
        }

        /* OK, it was a real collision, so find the next open spot */
        let retval = SUNHashMap_Iterate(map, idx + 1, sunHashMapLinearProbeInsert);
        if retval == SUNHASHMAP_ERROR {
            /* an error occurred */
            return retval;
        } else if retval == SUNHashMap_Capacity(map) {
            /* the map is out of empty buckets, so we grow it */
            let err = sunHashMapResize(map);
            if err != SUN_SUCCESS {
                return err as i64;
            }
            return SUNHashMap_Insert(map, key, value);
        }

        idx = retval;
    }

    /* Create and insert the key-value pair */
    map.buckets[idx as usize] = Some(SUNHashMapKeyValue {
        key: key.to_string(),
        value,
    });

    /* SUNStlVector_SUNHashMapKeyValue_Set returns SUN_SUCCESS (0) here */
    0
}

fn sunHashMapLinearProbeGet<V>(idx: i64, kv: &Option<SUNHashMapKeyValue<V>>, key: &str) -> i64 {
    /* find the matching entry */
    match kv {
        None => SUNHASHMAP_ERROR, /* keep looking since this bucket is empty */
        Some(kvp) => {
            if kvp.key == key {
                idx /* found it at idx */
            } else {
                SUNHASHMAP_ERROR /* keep looking */
            }
        }
    }
}

/*
  Shared probing sequence of SUNHashMap_GetValue / SUNHashMap_Remove (the two
  C functions contain this identical code). Returns the bucket index holding
  ``key``, SUNHASHMAP_ERROR, or an out-of-range index (== capacity) when the
  key was not found anywhere.
*/
fn sunHashMapProbeIdx<V>(map: &SUNHashMap<V>, key: &str) -> i64 {
    let mut idx = sunHashMapIdxFromKey(map, key);

    let kvp = &map.buckets[idx as usize];

    /* Check for a collision (an empty bucket means there was a collision at one
    point, but the colliding key has since been removed) */
    let collision = match kvp {
        Some(kv) => kv.key != key,
        None => true,
    };

    /* Resolve a collision via linear probing */
    if collision {
        let retval = SUNHashMap_Iterate(map, idx + 1, |i, kv| {
            sunHashMapLinearProbeGet(i, kv, key)
        });
        if retval == SUNHASHMAP_ERROR {
            /* the key was either not found anywhere or an error occurred */
            return retval;
        } else {
            idx = retval;
        }
    }

    idx
}

/*
  This function gets the value for the given key, returning a shared
  reference to it through the ``value`` out-parameter (the C code returned
  the aliased `void*`). On failure ``value`` is left untouched, like C.

  **Returns:**
    * ``0`` -- success
    * ``SUNHASHMAP_ERROR`` -- an error occurred
    * ``SUNHASHMAP_KEYNOTFOUND`` -- key not found
*/
pub fn SUNHashMap_GetValue<'a, V>(
    map: &'a SUNHashMap<V>,
    key: &str,
    value: &mut Option<&'a V>,
) -> i64 {
    let idx = sunHashMapProbeIdx(map, key);
    if idx == SUNHASHMAP_ERROR {
        return idx;
    }

    /* Return a reference to the value only */
    match map.buckets.get(idx as usize).and_then(|b| b.as_ref()) {
        Some(kvp) => {
            *value = Some(&kvp.value);
            0
        }
        None => SUNHASHMAP_KEYNOTFOUND,
    }
}

impl<V> SUNHashMap<V> {
    /// Rust-side mutable counterpart of `SUNHashMap_GetValue`. The C clients
    /// (e.g. the logger) simply alias the returned `void*` and mutate the
    /// pointee; safe Rust needs an explicit `&mut` accessor for that.
    pub fn get_mut(&mut self, key: &str) -> Option<&mut V> {
        let idx = sunHashMapProbeIdx(self, key);
        if idx < 0 {
            return None;
        }
        self.buckets
            .get_mut(idx as usize)
            .and_then(|b| b.as_mut())
            .map(|kvp| &mut kvp.value)
    }
}

/*
  This function removes the key-value pair, returning ownership of the value
  through the ``value`` out-parameter.

  **Returns:**
    * ``0`` -- success
    * ``SUNHASHMAP_ERROR`` -- an error occurred
    * ``SUNHASHMAP_KEYNOTFOUND`` -- key not found
*/
pub fn SUNHashMap_Remove<V>(map: &mut SUNHashMap<V>, key: &str, value: &mut Option<V>) -> i64 {
    let idx = sunHashMapProbeIdx(map, key);
    if idx == SUNHASHMAP_ERROR {
        return idx;
    }

    /* Clear the bucket, handing the value back to the caller (the C code
    dereferences the bucket unconditionally here; a missing key is reported
    as SUNHASHMAP_KEYNOTFOUND instead of invoking undefined behavior). */
    match map.buckets.get_mut(idx as usize) {
        Some(bucket) => match bucket.take() {
            Some(kvp) => {
                *value = Some(kvp.value);
                /* SUNStlVector_SUNHashMapKeyValue_Set returns SUN_SUCCESS (0) */
                0
            }
            None => SUNHASHMAP_KEYNOTFOUND,
        },
        None => SUNHASHMAP_KEYNOTFOUND,
    }
}

/*
  This function fills ``sorted`` with the map's buckets (including the empty
  ones, exactly like the C code) ordered by ``compar``, leaving the map
  unchanged. The C qsort comparator (returning an int) becomes a comparator
  returning `Ordering`.

  **Returns:**
    * A SUNErrCode indicating success or a failure
*/
pub fn SUNHashMap_Sort<'a, V>(
    map: &'a SUNHashMap<V>,
    sorted: &mut Vec<Option<&'a SUNHashMapKeyValue<V>>>,
    compar: fn(&Option<&SUNHashMapKeyValue<V>>, &Option<&SUNHashMapKeyValue<V>>) -> Ordering,
) -> SUNErrCode {
    sorted.clear();

    /* Copy the buckets into a new array */
    for i in 0..SUNHashMap_Capacity(map) {
        sorted.push(map.buckets[i as usize].as_ref());
    }

    sorted.sort_by(|a, b| compar(a, b));

    SUN_SUCCESS
}

/*
  This function fills ``values`` with references to just the values of the
  map (the C `value_size` byte-count argument only fed malloc and is dropped).

  **Returns:**
    * A SUNErrCode indicating success or a failure
*/
pub fn SUNHashMap_Values<'a, V>(map: &'a SUNHashMap<V>, values: &mut Vec<&'a V>) -> SUNErrCode {
    values.clear();

    /* Copy the values into a new array */
    for i in 0..SUNHashMap_Capacity(map) {
        if let Some(kvp) = &map.buckets[i as usize] {
            values.push(&kvp.value);
        }
    }

    SUN_SUCCESS
}

/// SUNHashMap_PrintKeys (C `FILE*` becomes `&mut dyn Write`; fprintf status
/// is ignored exactly like the C code).
pub fn SUNHashMap_PrintKeys<V>(map: &SUNHashMap<V>, file: &mut dyn Write) -> SUNErrCode {
    /* Print keys into a new array */
    let _ = write!(file, "[");
    for i in 0..SUNHashMap_Capacity(map) {
        if let Some(kvp) = &map.buckets[i as usize] {
            let _ = write!(file, "{}, ", kvp.key);
        }
    }
    let _ = writeln!(file, "]");

    SUN_SUCCESS
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_invalid_capacity() {
        let mut map: Option<SUNHashMap<i32>> = None;
        assert_eq!(SUNHashMap_New(0, &mut map), SUN_ERR_ARG_OUTOFRANGE);
        assert_eq!(SUNHashMap_New(-3, &mut map), SUN_ERR_ARG_OUTOFRANGE);
        assert!(map.is_none());
    }

    #[test]
    fn test_insert_get_remove() {
        let mut map: Option<SUNHashMap<i32>> = None;
        assert_eq!(SUNHashMap_New(8, &mut map), SUN_SUCCESS);
        let map = map.as_mut().unwrap();
        assert_eq!(SUNHashMap_Capacity(map), 8);

        assert_eq!(SUNHashMap_Insert(map, "alpha", 1), 0);
        assert_eq!(SUNHashMap_Insert(map, "beta", 2), 0);
        assert_eq!(SUNHashMap_Insert(map, "gamma", 3), 0);

        /* duplicate keys are rejected */
        assert_eq!(SUNHashMap_Insert(map, "beta", 99), SUNHASHMAP_DUPLICATE);

        let mut v: Option<&i32> = None;
        assert_eq!(SUNHashMap_GetValue(map, "alpha", &mut v), 0);
        assert_eq!(v, Some(&1));
        let mut v: Option<&i32> = None;
        assert_eq!(SUNHashMap_GetValue(map, "beta", &mut v), 0);
        assert_eq!(v, Some(&2));
        let mut v: Option<&i32> = None;
        assert_eq!(
            SUNHashMap_GetValue(map, "delta", &mut v),
            SUNHASHMAP_KEYNOTFOUND
        );
        assert!(v.is_none());

        /* mutable access (Rust adaptation used by the logger) */
        *map.get_mut("gamma").unwrap() = 30;
        let mut v: Option<&i32> = None;
        assert_eq!(SUNHashMap_GetValue(map, "gamma", &mut v), 0);
        assert_eq!(v, Some(&30));

        /* remove hands the value back and leaves the key unfindable */
        let mut removed: Option<i32> = None;
        assert_eq!(SUNHashMap_Remove(map, "beta", &mut removed), 0);
        assert_eq!(removed, Some(2));
        let mut v: Option<&i32> = None;
        assert_eq!(
            SUNHashMap_GetValue(map, "beta", &mut v),
            SUNHASHMAP_KEYNOTFOUND
        );
        let mut removed: Option<i32> = None;
        assert_eq!(
            SUNHashMap_Remove(map, "beta", &mut removed),
            SUNHASHMAP_KEYNOTFOUND
        );
        assert!(removed.is_none());
    }

    #[test]
    fn test_resize_preserves_entries() {
        let mut map: Option<SUNHashMap<usize>> = None;
        assert_eq!(SUNHashMap_New(4, &mut map), SUN_SUCCESS);
        let map = map.as_mut().unwrap();

        let n = 50usize;
        for i in 0..n {
            let key = format!("key_number_{i}");
            assert_eq!(SUNHashMap_Insert(map, &key, i), 0, "insert {key}");
        }

        /* the map had to grow past its initial capacity of 4 */
        assert!(SUNHashMap_Capacity(map) >= n as i64);

        for i in 0..n {
            let key = format!("key_number_{i}");
            let mut v: Option<&usize> = None;
            assert_eq!(SUNHashMap_GetValue(map, &key, &mut v), 0, "get {key}");
            assert_eq!(v, Some(&i));
        }

        let mut values: Vec<&usize> = Vec::new();
        assert_eq!(SUNHashMap_Values(map, &mut values), SUN_SUCCESS);
        assert_eq!(values.len(), n);
        let mut copied: Vec<usize> = values.iter().map(|v| **v).collect();
        copied.sort_unstable();
        assert_eq!(copied, (0..n).collect::<Vec<_>>());
    }

    #[test]
    fn test_collisions_capacity_one() {
        /* capacity 1 forces every key to index 0 and exercises the
        probe + grow path repeatedly */
        let mut map: Option<SUNHashMap<i32>> = None;
        assert_eq!(SUNHashMap_New(1, &mut map), SUN_SUCCESS);
        let map = map.as_mut().unwrap();

        for i in 0..10 {
            let key = format!("k{i}");
            assert_eq!(SUNHashMap_Insert(map, &key, i), 0);
        }
        for i in 0..10 {
            let key = format!("k{i}");
            let mut v: Option<&i32> = None;
            assert_eq!(SUNHashMap_GetValue(map, &key, &mut v), 0, "get {key}");
            assert_eq!(v, Some(&i));
        }
    }

    #[test]
    fn test_sort_and_print_keys() {
        let mut map: Option<SUNHashMap<i32>> = None;
        assert_eq!(SUNHashMap_New(4, &mut map), SUN_SUCCESS);
        let map = map.as_mut().unwrap();
        assert_eq!(SUNHashMap_Insert(map, "b", 2), 0);
        assert_eq!(SUNHashMap_Insert(map, "a", 1), 0);
        assert_eq!(SUNHashMap_Insert(map, "c", 3), 0);

        fn by_key(
            a: &Option<&SUNHashMapKeyValue<i32>>,
            b: &Option<&SUNHashMapKeyValue<i32>>,
        ) -> Ordering {
            /* occupied buckets first, ordered by key; empties last */
            match (a, b) {
                (Some(x), Some(y)) => x.key.cmp(&y.key),
                (Some(_), None) => Ordering::Less,
                (None, Some(_)) => Ordering::Greater,
                (None, None) => Ordering::Equal,
            }
        }

        let mut sorted: Vec<Option<&SUNHashMapKeyValue<i32>>> = Vec::new();
        assert_eq!(SUNHashMap_Sort(map, &mut sorted, by_key), SUN_SUCCESS);
        assert_eq!(sorted.len(), SUNHashMap_Capacity(map) as usize);
        let keys: Vec<&str> = sorted
            .iter()
            .flatten()
            .map(|kvp| kvp.key.as_str())
            .collect();
        assert_eq!(keys, ["a", "b", "c"]);

        let mut out: Vec<u8> = Vec::new();
        assert_eq!(SUNHashMap_PrintKeys(map, &mut out), SUN_SUCCESS);
        let text = String::from_utf8(out).unwrap();
        assert!(text.starts_with('['));
        assert!(text.ends_with("]\n"));
        assert!(text.contains("a, "));
        assert!(text.contains("b, "));
        assert!(text.contains("c, "));
    }

    #[test]
    fn test_destroy() {
        let mut map: Option<SUNHashMap<String>> = None;
        assert_eq!(SUNHashMap_New(2, &mut map), SUN_SUCCESS);
        assert_eq!(
            SUNHashMap_Insert(map.as_mut().unwrap(), "x", "y".to_string()),
            0
        );
        assert_eq!(SUNHashMap_Destroy(&mut map), SUN_SUCCESS);
        assert!(map.is_none());
        /* destroying an empty handle is fine, like C's NULL check */
        assert_eq!(SUNHashMap_Destroy(&mut map), SUN_SUCCESS);
    }
}
