# How to Resume the SUNDIALS → Rust Port

A plain-language guide to picking this work back up in a future session.
It assumes you know how to open a terminal and start Claude Code, but
nothing about the project itself.

---

## 1. What "resuming" actually means here

This is a long software project: translating a big C library (SUNDIALS)
into the Rust programming language, one file at a time. It is far too big
to finish in one sitting, so it is done across many separate Claude Code
sessions.

Each session leaves behind two kinds of "notes to the next session" so the
work can continue without re-reading everything:

- **Files inside the project** — `PROGRESS.md` (a checklist of what is done
  and what is next) and `CLAUDE.md` (the rules). These live *in the code
  folder* and travel with the project.
- **Claude's memory** — a file called `sundials-rs-port.md` plus an index
  called `MEMORY.md`. These live *outside* the project, in Claude's private
  memory area, and are loaded automatically at the start of every new
  session.

"Resuming" means: you start a fresh Claude Code session, point it at the
right folder, and tell it to continue. It reads those notes and keeps going.

---

## 2. The one thing that will trip you up: there are TWO folders

On this Mac there are two folders that look almost identical:

| Folder | What it is |
|---|---|
| `/Users/nsh/Developer/code/rust/cvode/` | An **older / parallel copy**. Do not work here. |
| `/Users/nsh/Developer/code/rust/sundials/` | The **active one**. All real work happens here. |

They both contain a sub-folder called `sundials_rs/` (the actual code),
which is why it is easy to confuse them. **The live, up-to-date work is only
in the `sundials/` one.** If a future session starts editing files under
`cvode/`, it is in the wrong place and its work will be lost or
disconnected.

The memory note spells this out precisely so the next session does not get
confused and start editing the wrong copy.

---

## 3. Step-by-step: how to resume

**Step 1 — Open a terminal and go to the active folder:**

```
cd /Users/nsh/Developer/code/rust/sundials
```

**Step 2 — Start Claude Code** (however you normally launch it — the
`claude` command in the terminal, or the desktop/IDE app opened on that
folder).

**Step 3 — Give it a short instruction.** You do not need to explain the
project; the notes do that. Something like:

> Resume the SUNDIALS→Rust port in sundials_rs/. Read CLAUDE.md,
> PROGRESS.md, and the recent git log first, then continue from the resume
> point.

That is enough, because:

- Claude's memory auto-loads and says "work happens in `rust/sundials/`, the
  IDA library is done, here is what is next."
- `CLAUDE.md` states the rules (no shortcuts, match the C behavior exactly,
  keep it warning-free).
- `PROGRESS.md` has the literal checklist with the next task marked.

---

## 4. What the next session will (correctly) do next

If it reads the notes properly, it will know the current state without your
explaining anything:

- The **IDA library is finished** (that was completed in the last session).
- The **next task** is writing and verifying the 8 IDA *example programs*
  (small demo programs that prove the library produces byte-for-byte
  identical results to the original C). `PROGRESS.md` lists them by name:
  `idaRoberts_dns`, `idaAnalytic_mels`, `idaFoodWeb_bnd`, `idaFoodWeb_kry`,
  `idaHeat2D_bnd`, `idaHeat2D_kry`, `idaKrylovDemo_ls`, `idaSlCrank_dns`.
- After that comes the next large component (called **idas**), then
  **arkode**, then documentation.

---

## 5. How to check the state for yourself

You do not have to take the next session's word for it. From inside
`/Users/nsh/Developer/code/rust/sundials/sundials_rs`, you (or it) can run
these and read the output.

**See the recent history of finished work:**

```
git -C /Users/nsh/Developer/code/rust/sundials/sundials_rs log --oneline -8
```

You should see the last session's commits at the top — their messages
mention `ida_io.c`, `ida_ic.c`, `ida_cli.c`, and `ida_bbdpre.c`, followed by
a `PROGRESS:` update.

**Prove the IDA code still works** (this should report `17 passed`):

```
cd /Users/nsh/Developer/code/rust/sundials/sundials_rs
cargo test -p ida_rs
```

**Read the checklist of what is done vs. to-do:**

```
open /Users/nsh/Developer/code/rust/sundials/sundials_rs/PROGRESS.md
```

---

## 6. One "gotcha" the next session must not panic about

If anyone runs a command that tries to build the *whole* project at once:

```
cargo build --workspace
```

…it **will show errors.** That is *expected and normal right now.* Those
errors come from a different, still-unfinished component called `cvodes_rs`
(a part of the project left half-done in an earlier phase). They are **not**
a regression from recent work.

The correct way to build or test while the project is mid-construction is one
component at a time, using the `-p` flag — for example `cargo test -p ida_rs`
or `cargo test -p sundials_core`. Both the memory note and `PROGRESS.md`
record this so the next session does not waste time "fixing" a problem that
is not a regression.

---

## 7. In one sentence

To resume: **open a terminal, run `cd /Users/nsh/Developer/code/rust/sundials`,
start Claude Code there, and say "resume the SUNDIALS→Rust port, read
CLAUDE.md and PROGRESS.md first"** — the notes left behind will carry it
forward from exactly where the last session stopped (IDA library done; IDA
example programs are next), and just remember that the `sundials/` folder is
the real one, not the near-identical `cvode/` folder.
