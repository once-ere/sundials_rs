/* -----------------------------------------------------------------
 * Translated from src/sunlinsol/pcg/sunlinsol_pcg.c
 * (SUNDIALS 7.7.0): preconditioned conjugate gradient (symmetric
 * systems).
 *
 * The C content struct becomes PcgLS; the ops-table entry points
 * become methods dispatched from sundials_linearsolver.rs. ATimes /
 * PSolve arrive as closures at solve time; a required-but-missing
 * psolve yields SUNLS_PSOLVE_NULL. PCG applies symmetric
 * preconditioning: pretype is stored as given, and any of
 * SUN_PREC_LEFT/RIGHT/BOTH enables the preconditioner (always
 * invoked with lr = SUN_PREC_LEFT). Only the first scaling vector
 * (s1, the C content->s) is used; s2 is ignored, matching
 * SUNLinSolSetScalingVectors_PCG.
 * -----------------------------------------------------------------*/
use crate::nvector_serial::{NVector, N_VClone, N_VDotProd, N_VProd, N_VScale};
use crate::sundials_context::SUNContext;
use crate::sundials_errors::{SUNErrCode, SUN_ERR_ARG_OUTOFRANGE, SUN_SUCCESS};
use crate::sundials_linearsolver::{
    ATimesFn, LinearSolver, PSolveFn, SUNLS_ATIMES_FAIL_REC, SUNLS_ATIMES_FAIL_UNREC,
    SUNLS_CONV_FAIL, SUNLS_PSOLVE_FAIL_REC, SUNLS_PSOLVE_FAIL_UNREC, SUNLS_PSOLVE_NULL,
    SUNLS_RES_REDUCED, SUN_PREC_BOTH, SUN_PREC_LEFT, SUN_PREC_NONE, SUN_PREC_RIGHT,
};
use crate::sundials_math::SUNRsqrt;

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

/* Default PCG parameters (sunlinsol_pcg.h) */
pub const SUNPCG_MAXL_DEFAULT: i32 = 5;

/// SUNLinearSolverContent_PCG
pub struct PcgLS {
    pub numiters: i32,
    pub resnorm: f64,
    pub last_flag: i64,
    pub zeroguess: bool,
    maxl: i32,
    pretype: i32,
    /* workspace vectors, all allocated by the constructor as in C */
    r: NVector,
    p: NVector,
    z: NVector,
    ap: NVector,
}

/* ----------------------------------------------------------------------------
 * Function to create a new PCG linear solver (SUNLinSol_PCG)
 */
pub fn SUNLinSol_PCG(y: &NVector, pretype: i32, maxl: i32, _sunctx: &SUNContext) -> LinearSolver {
    /* check for legal pretype and maxl values; if illegal use defaults */
    let pretype = if (pretype != SUN_PREC_NONE)
        && (pretype != SUN_PREC_LEFT)
        && (pretype != SUN_PREC_RIGHT)
        && (pretype != SUN_PREC_BOTH)
    {
        SUN_PREC_NONE
    } else {
        pretype
    };
    let maxl = if maxl <= 0 { SUNPCG_MAXL_DEFAULT } else { maxl };

    LinearSolver::Pcg(PcgLS {
        numiters: 0,
        resnorm: ZERO,
        last_flag: 0,
        zeroguess: false,
        maxl,
        pretype,
        r: N_VClone(y),
        p: N_VClone(y),
        z: N_VClone(y),
        ap: N_VClone(y),
    })
}

impl PcgLS {
    /// SUNLinSolSpace_PCG: (lenrwLS, leniwLS)
    pub fn space(&self) -> (i64, i64) {
        let (lrw1, liw1) = crate::nvector_serial::N_VSpace(&self.r);
        (1 + lrw1 * 4, 4 + liw1 * 4)
    }

    /// SUNLinSolInitialize_PCG
    pub fn initialize(&mut self) -> SUNErrCode {
        if self.maxl <= 0 {
            self.maxl = SUNPCG_MAXL_DEFAULT;
        }

        if (self.pretype != SUN_PREC_LEFT)
            && (self.pretype != SUN_PREC_RIGHT)
            && (self.pretype != SUN_PREC_BOTH)
        {
            self.pretype = SUN_PREC_NONE;
        }

        /* no additional memory to allocate */
        SUN_SUCCESS
    }

    /// SUNLinSol_PCGSetPrecType
    pub fn set_prec_type(&mut self, pretype: i32) -> SUNErrCode {
        /* Check for legal pretype */
        if (pretype != SUN_PREC_NONE)
            && (pretype != SUN_PREC_LEFT)
            && (pretype != SUN_PREC_RIGHT)
            && (pretype != SUN_PREC_BOTH)
        {
            return SUN_ERR_ARG_OUTOFRANGE;
        }
        self.pretype = pretype;
        SUN_SUCCESS
    }

    /// SUNLinSol_PCGSetMaxl
    pub fn set_maxl(&mut self, maxl: i32) -> SUNErrCode {
        /* Check for legal number of iters */
        self.maxl = if maxl <= 0 { SUNPCG_MAXL_DEFAULT } else { maxl };
        SUN_SUCCESS
    }

    /// SUNLinSolSolve_PCG
    pub fn solve(
        &mut self,
        x: &mut NVector,
        b: &NVector,
        delta: f64,
        atimes: &mut ATimesFn,
        mut psolve: Option<&mut PSolveFn>,
        s1: Option<&NVector>,
        _s2: Option<&NVector>,
    ) -> i32 {
        /* Make local shortcuts to solver variables. */
        let l_max = self.maxl;
        let pretype = self.pretype;
        let PcgLS {
            numiters,
            resnorm,
            last_flag,
            zeroguess,
            r,
            p,
            z,
            ap,
            ..
        } = self;
        let w = s1; /* C: content->s */

        /* Initialize counters and convergence flag */
        *numiters = 0;
        let mut converged = false;

        /* set flags for internal solver options */
        let use_prec = (pretype == SUN_PREC_BOTH)
            || (pretype == SUN_PREC_LEFT)
            || (pretype == SUN_PREC_RIGHT);
        let use_scaling = w.is_some();

        /* If preconditioning, check if psolve has been set */
        if use_prec && psolve.is_none() {
            *last_flag = SUNLS_PSOLVE_NULL as i64;
            return SUNLS_PSOLVE_NULL;
        }

        /* Set r to initial residual r_0 = b - A*x_0 */
        if *zeroguess {
            N_VScale(ONE, b, r);
        } else {
            let status = atimes(x, r);
            if status != 0 {
                *zeroguess = false;
                *last_flag = if status < 0 {
                    SUNLS_ATIMES_FAIL_UNREC
                } else {
                    SUNLS_ATIMES_FAIL_REC
                } as i64;
                return *last_flag as i32;
            }
            /* r = b - r */
            r.linear_sum_with(-ONE, ONE, b);
        }

        /* Set rho to scaled L2 norm of r, and return if small */
        if use_scaling {
            N_VProd(r, w.unwrap(), ap);
        } else {
            N_VScale(ONE, r, ap);
        }
        let mut rho = N_VDotProd(ap, ap);
        rho = SUNRsqrt(rho);
        let r0_norm = rho;
        *resnorm = rho;

        if rho <= delta {
            *zeroguess = false;
            *last_flag = SUN_SUCCESS as i64;
            return SUN_SUCCESS;
        }

        /* Apply preconditioner and b-scaling to r = r_0 */
        if use_prec {
            /* z = P^{-1}r */
            let status = (psolve.as_mut().unwrap())(r, z, delta, SUN_PREC_LEFT);
            if status != 0 {
                *zeroguess = false;
                *last_flag = if status < 0 {
                    SUNLS_PSOLVE_FAIL_UNREC
                } else {
                    SUNLS_PSOLVE_FAIL_REC
                } as i64;
                return *last_flag as i32;
            }
        } else {
            N_VScale(ONE, r, z);
        }

        /* Initialize rz to <r,z> */
        let mut rz = N_VDotProd(r, z);

        /* Copy z to p */
        N_VScale(ONE, z, p);

        /* Begin main iteration loop */
        for l in 0..l_max {
            /* increment counter */
            *numiters += 1;

            /* Generate Ap = A*p */
            let status = atimes(p, ap);
            if status != 0 {
                *zeroguess = false;
                *last_flag = if status < 0 {
                    SUNLS_ATIMES_FAIL_UNREC
                } else {
                    SUNLS_ATIMES_FAIL_REC
                } as i64;
                return *last_flag as i32;
            }

            /* Calculate alpha = <r,z> / <Ap,p> */
            let mut alpha = N_VDotProd(ap, p);
            alpha = rz / alpha;

            /* Update x = x + alpha*p */
            if l == 0 && *zeroguess {
                N_VScale(alpha, p, x);
            } else {
                /* N_VLinearSum(ONE, x, alpha, p, x) */
                x.linear_sum_with(ONE, alpha, p);
            }

            /* Update r = r - alpha*Ap */
            r.linear_sum_with(ONE, -alpha, ap);

            /* Set rho and check convergence */
            if use_scaling {
                N_VProd(r, w.unwrap(), ap);
            } else {
                N_VScale(ONE, r, ap);
            }
            rho = N_VDotProd(ap, ap);
            rho = SUNRsqrt(rho);
            *resnorm = rho;

            if rho <= delta {
                converged = true;
                break;
            }

            /* Exit early on last iteration */
            if l == l_max - 1 {
                break;
            }

            /* Apply preconditioner:  z = P^{-1}*r */
            if use_prec {
                let status = (psolve.as_mut().unwrap())(r, z, delta, SUN_PREC_LEFT);
                if status != 0 {
                    *zeroguess = false;
                    *last_flag = if status < 0 {
                        SUNLS_PSOLVE_FAIL_UNREC
                    } else {
                        SUNLS_PSOLVE_FAIL_REC
                    } as i64;
                    return *last_flag as i32;
                }
            } else {
                N_VScale(ONE, r, z);
            }

            /* update rz */
            let rz_old = rz;
            rz = N_VDotProd(r, z);

            /* Calculate beta = <r,z> / <r_old,z_old> */
            let beta = rz / rz_old;

            /* Update p = z + beta*p */
            p.linear_sum_with(beta, ONE, z);
        }

        /* Main loop finished, return with result */
        *zeroguess = false;
        if converged {
            *last_flag = SUN_SUCCESS as i64;
        } else if rho < r0_norm {
            *last_flag = SUNLS_RES_REDUCED as i64;
        } else {
            *last_flag = SUNLS_CONV_FAIL as i64;
        }
        *last_flag as i32
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /* symmetric positive definite tridiagonal test matrix */
    fn matvec(v: &NVector, av: &mut NVector) {
        let n = v.len();
        for i in 0..n {
            let mut s = 4.0 * v.data[i];
            if i > 0 {
                s -= v.data[i - 1];
            }
            if i + 1 < n {
                s -= v.data[i + 1];
            }
            av.data[i] = s;
        }
    }

    fn residual_norm(x: &NVector, b: &NVector) -> f64 {
        let mut ax = NVector::new(x.len());
        matvec(x, &mut ax);
        let mut sum = 0.0;
        for i in 0..x.len() {
            let d = b.data[i] - ax.data[i];
            sum += d * d;
        }
        sum.sqrt()
    }

    fn rhs(n: usize) -> NVector {
        NVector::from_slice(&(0..n).map(|i| 1.0 + 0.1 * i as f64).collect::<Vec<_>>())
    }

    fn unwrap_pcg(ls: LinearSolver) -> PcgLS {
        match ls {
            LinearSolver::Pcg(s) => s,
            _ => unreachable!(),
        }
    }

    #[test]
    fn pcg_solves_spd_unpreconditioned() {
        let n = 12usize;
        let ctx = SUNContext::default();
        let y = NVector::new(n);
        let mut s = unwrap_pcg(SUNLinSol_PCG(&y, SUN_PREC_NONE, 200, &ctx));
        assert_eq!(s.initialize(), SUN_SUCCESS);

        let b = rhs(n);
        let mut x = NVector::new(n);
        s.zeroguess = true;
        let delta = 1e-10;
        let mut atimes = |v: &NVector, av: &mut NVector| -> i32 {
            matvec(v, av);
            0
        };
        let flag = s.solve(&mut x, &b, delta, &mut atimes, None, None, None);
        assert_eq!(flag, SUN_SUCCESS);
        assert!(s.numiters > 0);
        assert!(s.resnorm <= delta);
        assert!(residual_norm(&x, &b) <= 1e-8);
    }

    #[test]
    fn pcg_jacobi_preconditioner() {
        let n = 12usize;
        let ctx = SUNContext::default();
        let y = NVector::new(n);
        let mut s = unwrap_pcg(SUNLinSol_PCG(&y, SUN_PREC_LEFT, 200, &ctx));
        assert_eq!(s.initialize(), SUN_SUCCESS);

        let b = rhs(n);
        let mut x = NVector::new(n);
        s.zeroguess = true;
        let delta = 1e-10;
        let mut atimes = |v: &NVector, av: &mut NVector| -> i32 {
            matvec(v, av);
            0
        };
        /* Jacobi (diagonal) preconditioner: z = r / 4 */
        let mut psolve = |r: &NVector, z: &mut NVector, _tol: f64, _lr: i32| -> i32 {
            for (zi, ri) in z.data.iter_mut().zip(&r.data) {
                *zi = *ri / 4.0;
            }
            0
        };
        let flag = s.solve(
            &mut x,
            &b,
            delta,
            &mut atimes,
            Some(&mut psolve),
            None,
            None,
        );
        assert_eq!(flag, SUN_SUCCESS);
        assert!(residual_norm(&x, &b) <= 1e-8);
    }
}
