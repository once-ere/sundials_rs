/* -----------------------------------------------------------------
 * Translated from src/kinsol/kinsol_orth.c (KINSOL 7.7.0).
 * Orthogonalization utilities for Anderson acceleration: workspace
 * (re)allocation and selection of the SUNQRAdd kernel.
 *
 * Adaptations (workspace conventions):
 *  - malloc + failure unwinding become infallible Vec/Box
 *    construction (RAII); the KIN_MEM_FAIL paths are unreachable.
 *  - In C the SUNQRData members are *aliases* into KINSOL memory
 *    (vtemp -> kin_vtemp2, vtemp2 -> kin_vtemp3, temp_array ->
 *    kin_T_aa or kin_cv). The Rust SUNQRData owns its workspace
 *    (sundials_iterative.rs), so equally-sized owned buffers are
 *    created here. The one alias with cross-call state — the ICWY
 *    T matrix, updated both by SUNQRAdd_ICWY and by
 *    AndersonAccQRDelete — keeps kin_T_aa canonical: kinsol.rs
 *    swaps it into qr_data.temp_array around each kin_qr_func call.
 *  - The serial N_Vector has no nvdotprodmultiallreduce op, so the
 *    C dot_prod_sb probe is always SUNFALSE and the _SB QRAdd
 *    variants are never selected (they are not ported; see
 *    sundials_iterative.rs).
 * -----------------------------------------------------------------*/
use crate::kinsol_impl::*;
use crate::nvector_serial::{NVector, N_VClone};
use crate::sundials_iterative::{
    SUNQRAdd_CGS2, SUNQRAdd_DCGS2, SUNQRAdd_ICWY, SUNQRAdd_MGS, SUNQRData,
};
use crate::sundials_types::SUNFALSE;

/* KINInitOrth
 *
 * C: int KINInitOrth(KINMem kin_mem) */
pub fn KINInitOrth(kin_mem: &mut KINMem) -> i32 {
    /* Do we need to (re)allocate the orthogonalization workspace? */
    let allocate = kin_mem.kin_m_aa > kin_mem.kin_orth_aa_alloc;

    if allocate {
        /* Free any existing workspace allocations */
        KINFreeOrth(kin_mem);

        /* Template vector for creating clones: kin_mem->kin_unew */

        /* Update the AA workspace size */
        kin_mem.kin_orth_aa_alloc = kin_mem.kin_m_aa;

        let m_aa = kin_mem.kin_m_aa as usize;

        /* Structure of orthogonalization data for QR solve */
        kin_mem.kin_qr_data = Some(SUNQRData::default());

        if kin_mem.kin_orth_aa != KIN_ORTH_MGS {
            kin_mem.kin_vtemp3 = N_VClone(&kin_mem.kin_unew); /* Orth owns */
        }

        if kin_mem.kin_orth_aa == KIN_ORTH_ICWY {
            /* T matrix for ICWY */
            kin_mem.kin_T_aa = vec![0.0; m_aa * m_aa];
        }
    }

    /* Does the vector support dot product with single buffer reductions?
       The serial N_Vector implements nvdotprodlocal/nvdotprodmultilocal
       but not nvdotprodmultiallreduce, so this is always SUNFALSE. */
    kin_mem.kin_dot_prod_sb = SUNFALSE;

    /* Initialize the QRData and set the QRAdd function.
       (C re-points the SUNQRData aliases on every call; here the owned
       buffers are refreshed to the current problem size.) */
    let mut qr_data = kin_mem.kin_qr_data.take().unwrap_or_default();
    if kin_mem.kin_orth_aa == KIN_ORTH_MGS {
        kin_mem.kin_qr_func = Some(SUNQRAdd_MGS);
        qr_data.vtemp = N_VClone(&kin_mem.kin_vtemp2); /* KINSOL owns */
    } else if kin_mem.kin_orth_aa == KIN_ORTH_ICWY {
        /* dot_prod_sb == SUNFALSE always: SUNQRAdd_ICWY_SB unreachable */
        kin_mem.kin_qr_func = Some(SUNQRAdd_ICWY);
        qr_data.vtemp = N_VClone(&kin_mem.kin_vtemp2); /* KINSOL owns */
        qr_data.vtemp2 = N_VClone(&kin_mem.kin_vtemp3); /* Orth owns */
        /* temp_array -> kin_T_aa: kin_T_aa stays canonical and is
           swapped in around each kin_qr_func call (kinsol.rs) */
        qr_data.temp_array = Vec::new();
    } else if kin_mem.kin_orth_aa == KIN_ORTH_CGS2 {
        kin_mem.kin_qr_func = Some(SUNQRAdd_CGS2);
        qr_data.vtemp = N_VClone(&kin_mem.kin_vtemp2); /* KINSOL owns */
        qr_data.vtemp2 = N_VClone(&kin_mem.kin_vtemp3); /* Orth owns */
        /* temp_array -> kin_cv scratch (2*(m_aa+1) sunrealtype) */
        qr_data.temp_array = vec![0.0; 2 * (kin_mem.kin_m_aa as usize + 1)];
    } else if kin_mem.kin_orth_aa == KIN_ORTH_DCGS2 {
        /* dot_prod_sb == SUNFALSE always: SUNQRAdd_DCGS2_SB unreachable */
        kin_mem.kin_qr_func = Some(SUNQRAdd_DCGS2);
        qr_data.vtemp = N_VClone(&kin_mem.kin_vtemp2); /* KINSOL owns */
        qr_data.vtemp2 = N_VClone(&kin_mem.kin_vtemp3); /* Orth owns */
        /* temp_array -> kin_cv scratch (2*(m_aa+1) sunrealtype) */
        qr_data.temp_array = vec![0.0; 2 * (kin_mem.kin_m_aa as usize + 1)];
    }
    kin_mem.kin_qr_data = Some(qr_data);

    KIN_SUCCESS
}

/* KINFreeOrth
 *
 * C: void KINFreeOrth(KINMem kin_mem) */
pub fn KINFreeOrth(kin_mem: &mut KINMem) {
    kin_mem.kin_qr_data = None;
    kin_mem.kin_vtemp3 = NVector::default();
    kin_mem.kin_T_aa = Vec::new();

    /* Reset AA workspace size */
    kin_mem.kin_orth_aa_alloc = 0;
}
