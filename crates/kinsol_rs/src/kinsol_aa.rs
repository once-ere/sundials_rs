/* -----------------------------------------------------------------
 * Translated from src/kinsol/kinsol_aa.c (KINSOL 7.7.0).
 * Anderson acceleration utilities: (re)allocation and release of the
 * acceleration workspace attached to KINMem.
 *
 * Adaptations (workspace conventions):
 *  - malloc/N_VCloneVectorArray + the C malloc-failure unwinding
 *    become infallible Vec/NVector construction (RAII); the
 *    KIN_MEM_FAIL error paths are structurally unreachable.
 *  - The kin_cv / kin_Xv fused-operation scratch arrays are not
 *    stored (see kinsol_impl.rs): the serial fused kernels are
 *    reproduced inline in kinsol.rs::AndersonAcc.
 * -----------------------------------------------------------------*/
use crate::kinsol_impl::*;
use crate::nvector_serial::{NVector, N_VClone};

/* KINInitAA
 *
 * C: int KINInitAA(KINMem kin_mem) */
pub fn KINInitAA(kin_mem: &mut KINMem) -> i32 {
    /* Limit the acceleration space size */
    if kin_mem.kin_m_aa >= kin_mem.kin_mxiter {
        kin_mem.kin_m_aa = kin_mem.kin_mxiter - 1;
    }

    /* Initialize the current depth */
    kin_mem.kin_current_depth = 0;

    /* Do we need to (re)allocate the AA workspace? */
    let allocate = kin_mem.kin_m_aa > kin_mem.kin_m_aa_alloc;

    if allocate {
        /* Free any existing workspace allocations */
        KINFreeAA(kin_mem);

        /* Template vector for creating clones: kin_mem->kin_unew */

        /* Update the AA workspace size */
        kin_mem.kin_m_aa_alloc = kin_mem.kin_m_aa;

        let m_aa = kin_mem.kin_m_aa as usize;

        /* Array of acceleration weights */
        kin_mem.kin_gamma_aa = vec![0.0; m_aa];

        /* R matrix for QR factorization */
        kin_mem.kin_R_aa = vec![0.0; m_aa * m_aa];

        /* Q matrix for QR factorization */
        kin_mem.kin_q_aa = (0..m_aa).map(|_| N_VClone(&kin_mem.kin_unew)).collect();

        /* History of residual vector differences */
        kin_mem.kin_df_aa = (0..m_aa).map(|_| N_VClone(&kin_mem.kin_unew)).collect();

        /* History of fixed point function vector differences */
        kin_mem.kin_dg_aa = (0..m_aa).map(|_| N_VClone(&kin_mem.kin_unew)).collect();

        /* Previous residual vector, F(u_{i-1}) = G(u_{i-1}) - u_{i-1} */
        kin_mem.kin_fold_aa = N_VClone(&kin_mem.kin_unew);

        /* Previous fixed point function vector, G(u_{i-1}) */
        kin_mem.kin_gold_aa = N_VClone(&kin_mem.kin_unew);

        /* (C also mallocs the kin_cv / kin_Xv workspace arrays of size
           2*(m_aa+1) for the fused N_VLinearCombination update; those
           are built inline in kinsol.rs::AndersonAcc.) */
    }

    KIN_SUCCESS
}

/* KINFreeAA
 *
 * C: void KINFreeAA(KINMem kin_mem) */
pub fn KINFreeAA(kin_mem: &mut KINMem) {
    kin_mem.kin_gamma_aa = Vec::new();
    kin_mem.kin_R_aa = Vec::new();
    kin_mem.kin_q_aa = Vec::new();
    kin_mem.kin_df_aa = Vec::new();
    kin_mem.kin_dg_aa = Vec::new();
    kin_mem.kin_fold_aa = NVector::default();
    kin_mem.kin_gold_aa = NVector::default();

    /* Reset AA workspace size */
    kin_mem.kin_m_aa_alloc = 0;
}
