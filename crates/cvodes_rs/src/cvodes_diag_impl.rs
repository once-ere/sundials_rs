/* -----------------------------------------------------------------
 * Translated from src/cvodes/cvodes_diag_impl.h (CVODES 7.7.0).
 * CVDIAG diagonal linear solver memory structure.
 * -----------------------------------------------------------------*/
use crate::nvector_serial::NVector;

/* CVDIAG return values (cvodes_diag.h) */
pub const CVDIAG_SUCCESS: i32 = 0;
pub const CVDIAG_MEM_NULL: i32 = -1;
pub const CVDIAG_LMEM_NULL: i32 = -2;
pub const CVDIAG_ILL_INPUT: i32 = -3;
pub const CVDIAG_MEM_FAIL: i32 = -4;

/* Additional last_flag values */
pub const CVDIAG_INV_FAIL: i32 = -5;
pub const CVDIAG_RHSFUNC_UNRECVR: i32 = -6;
pub const CVDIAG_RHSFUNC_RECVR: i32 = -7;

/* Return values for adjoint module */
pub const CVDIAG_NO_ADJ: i32 = -101;

/* -----------------------------------------------------------------
   Types: CVDiagMemRec, CVDiagMem
   -----------------------------------------------------------------*/
pub struct CVDiagMem {
    pub di_gammasv: f64, /* gammasv = gamma at the last call to setup or solve */
    pub di_M: NVector,   /* M = (I - gamma J)^{-1} , gamma = h / l1 */
    pub di_bit: NVector, /* temporary storage vector */
    pub di_bitcomp: NVector, /* temporary storage vector */
    pub di_nfeDI: i64,   /* no. of f calls for difference-quotient diagonal Jacobian */
    pub di_last_flag: i64, /* last error return flag */
}
