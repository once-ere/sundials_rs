/* -----------------------------------------------------------------
 * Translated from src/sunlinsol/spbcgs/sunlinsol_spbcgs.c
 * (SUNDIALS 7.7.0): scaled preconditioned Bi-CGStab.
 *
 * The C content struct becomes SpbcgsLS; the ops-table entry points
 * become methods dispatched from sundials_linearsolver.rs. ATimes /
 * PSolve arrive as closures at solve time; a required-but-missing
 * psolve yields SUNLS_PSOLVE_NULL. As in C, right preconditioning
 * with a nonzero initial guess is an unsupported configuration and
 * returns SUN_ERR_ARG_INCOMPATIBLE. The fused N_VLinearCombination
 * kernels are reproduced inline with the exact floating-point
 * operation order of the serial implementation.
 * -----------------------------------------------------------------*/
use crate::nvector_serial::{NVector, N_VClone, N_VDiv, N_VDotProd, N_VLinearSum, N_VProd, N_VScale};
use crate::sundials_context::SUNContext;
use crate::sundials_errors::{SUNErrCode, SUN_ERR_ARG_CORRUPT, SUN_ERR_ARG_INCOMPATIBLE, SUN_SUCCESS};
use crate::sundials_linearsolver::{
    ATimesFn, LinearSolver, PSolveFn, SUNLS_ATIMES_FAIL_REC, SUNLS_ATIMES_FAIL_UNREC,
    SUNLS_CONV_FAIL, SUNLS_PSOLVE_FAIL_REC, SUNLS_PSOLVE_FAIL_UNREC, SUNLS_PSOLVE_NULL,
    SUNLS_RES_REDUCED, SUN_PREC_BOTH, SUN_PREC_LEFT, SUN_PREC_NONE, SUN_PREC_RIGHT,
};
use crate::sundials_math::SUNRsqrt;

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

/* Default SPBCGS parameters (sunlinsol_spbcgs.h) */
pub const SUNSPBCGS_MAXL_DEFAULT: i32 = 5;

/// SUNLinearSolverContent_SPBCGS
pub struct SpbcgsLS {
    pub numiters: i32,
    pub resnorm: f64,
    pub last_flag: i64,
    pub zeroguess: bool,
    maxl: i32,
    pretype: i32,
    /* workspace vectors, all allocated by the constructor as in C */
    r_star: NVector,
    r: NVector,
    p: NVector,
    q: NVector,
    u: NVector,
    ap: NVector,
    vtemp: NVector,
}

/* ----------------------------------------------------------------------------
 * Function to create a new SPBCGS linear solver (SUNLinSol_SPBCGS)
 */
pub fn SUNLinSol_SPBCGS(y: &NVector, pretype: i32, maxl: i32, _sunctx: &SUNContext) -> LinearSolver {
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
    let maxl = if maxl <= 0 { SUNSPBCGS_MAXL_DEFAULT } else { maxl };

    LinearSolver::Spbcgs(SpbcgsLS {
        numiters: 0,
        resnorm: ZERO,
        last_flag: 0,
        zeroguess: false,
        maxl,
        pretype,
        r_star: N_VClone(y),
        r: N_VClone(y),
        p: N_VClone(y),
        q: N_VClone(y),
        u: N_VClone(y),
        ap: N_VClone(y),
        vtemp: N_VClone(y),
    })
}

impl SpbcgsLS {
    /// SUNLinSolSpace_SPBCGS: (lenrwLS, leniwLS)
    pub fn space(&self) -> (i64, i64) {
        let (lrw1, liw1) = crate::nvector_serial::N_VSpace(&self.r_star);
        (lrw1 * 9, liw1 * 9)
    }

    /// SUNLinSolInitialize_SPBCGS
    pub fn initialize(&mut self) -> SUNErrCode {
        if self.maxl <= 0 {
            self.maxl = SUNSPBCGS_MAXL_DEFAULT;
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

    /// SUNLinSol_SPBCGSSetPrecType
    pub fn set_prec_type(&mut self, pretype: i32) -> SUNErrCode {
        /* Check for legal pretype */
        if (pretype != SUN_PREC_NONE)
            && (pretype != SUN_PREC_LEFT)
            && (pretype != SUN_PREC_RIGHT)
            && (pretype != SUN_PREC_BOTH)
        {
            return SUN_ERR_ARG_CORRUPT;
        }
        self.pretype = pretype;
        SUN_SUCCESS
    }

    /// SUNLinSol_SPBCGSSetMaxl
    pub fn set_maxl(&mut self, maxl: i32) -> SUNErrCode {
        self.maxl = if maxl <= 0 { SUNSPBCGS_MAXL_DEFAULT } else { maxl };
        SUN_SUCCESS
    }

    /// SUNLinSolSolve_SPBCGS
    pub fn solve(
        &mut self,
        x: &mut NVector,
        b: &NVector,
        delta: f64,
        atimes: &mut ATimesFn,
        mut psolve: Option<&mut PSolveFn>,
        s1: Option<&NVector>,
        s2: Option<&NVector>,
    ) -> i32 {
        /* Make local shortcuts to solver variables. */
        let l_max = self.maxl;
        let pretype = self.pretype;
        let SpbcgsLS {
            numiters,
            resnorm,
            last_flag,
            zeroguess,
            r_star,
            r,
            p,
            q,
            u,
            ap,
            vtemp,
            ..
        } = self;
        let sb = s1;
        let sx = s2;

        /* Initialize counters and convergence flag */
        *numiters = 0;
        let mut converged = false;

        /* set flags for internal solver options */
        let pre_on_left = (pretype == SUN_PREC_LEFT) || (pretype == SUN_PREC_BOTH);
        let pre_on_right = (pretype == SUN_PREC_RIGHT) || (pretype == SUN_PREC_BOTH);
        let scale_x = sx.is_some();
        let scale_b = sb.is_some();

        /* Check for unsupported use case */
        if pre_on_right && !*zeroguess {
            *zeroguess = false;
            *last_flag = SUN_ERR_ARG_INCOMPATIBLE as i64;
            return SUN_ERR_ARG_INCOMPATIBLE;
        }

        /* If preconditioning, check if psolve has been set */
        if (pre_on_left || pre_on_right) && psolve.is_none() {
            *last_flag = SUNLS_PSOLVE_NULL as i64;
            return SUNLS_PSOLVE_NULL;
        }

        /* Set r_star to initial (unscaled) residual r_0 = b - A*x_0 */
        if *zeroguess {
            N_VScale(ONE, b, r_star);
        } else {
            let status = atimes(x, r_star);
            if status != 0 {
                *zeroguess = false;
                *last_flag = if status < 0 {
                    SUNLS_ATIMES_FAIL_UNREC
                } else {
                    SUNLS_ATIMES_FAIL_REC
                } as i64;
                return *last_flag as i32;
            }
            /* r_star = b - r_star */
            r_star.linear_sum_with(-ONE, ONE, b);
        }

        /* Apply left preconditioner and b-scaling to r_star = r_0 */
        if pre_on_left {
            let status = (psolve.as_mut().unwrap())(r_star, r, delta, SUN_PREC_LEFT);
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
            N_VScale(ONE, r_star, r);
        }

        if scale_b {
            N_VProd(sb.unwrap(), r, r_star);
        } else {
            N_VScale(ONE, r, r_star);
        }

        /* Initialize beta_denom to the dot product of r0 with r0 */
        let mut beta_denom = N_VDotProd(r_star, r_star);

        /* Set r_norm to L2 norm of r_star = sb P1_inv r_0, and
           return if small */
        let r_norm = SUNRsqrt(beta_denom);
        let mut rho = r_norm;
        *resnorm = r_norm;

        if r_norm <= delta {
            *zeroguess = false;
            *last_flag = SUN_SUCCESS as i64;
            return SUN_SUCCESS;
        }

        /* Copy r_star to r and p */
        N_VScale(ONE, r_star, r);
        N_VScale(ONE, r_star, p);

        /* Set x = sx x if non-zero guess */
        if scale_x && !*zeroguess {
            /* N_VProd(sx, x, x) */
            x.prod_with(sx.unwrap());
        }

        /* Begin main iteration loop */
        for l in 0..l_max {
            *numiters += 1;

            /* Generate Ap = A-tilde p, where A-tilde = sb P1_inv A P2_inv sx_inv */

            /*   Apply x-scaling: vtemp = sx_inv p */
            if scale_x {
                N_VDiv(p, sx.unwrap(), vtemp);
            } else {
                N_VScale(ONE, p, vtemp);
            }

            /*   Apply right preconditioner: vtemp = P2_inv sx_inv p */
            if pre_on_right {
                N_VScale(ONE, vtemp, ap);
                let status = (psolve.as_mut().unwrap())(ap, vtemp, delta, SUN_PREC_RIGHT);
                if status != 0 {
                    *zeroguess = false;
                    *last_flag = if status < 0 {
                        SUNLS_PSOLVE_FAIL_UNREC
                    } else {
                        SUNLS_PSOLVE_FAIL_REC
                    } as i64;
                    return *last_flag as i32;
                }
            }

            /*   Apply A: Ap = A P2_inv sx_inv p */
            let status = atimes(vtemp, ap);
            if status != 0 {
                *zeroguess = false;
                *last_flag = if status < 0 {
                    SUNLS_ATIMES_FAIL_UNREC
                } else {
                    SUNLS_ATIMES_FAIL_REC
                } as i64;
                return *last_flag as i32;
            }

            /*   Apply left preconditioner: vtemp = P1_inv A P2_inv sx_inv p */
            if pre_on_left {
                let status = (psolve.as_mut().unwrap())(ap, vtemp, delta, SUN_PREC_LEFT);
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
                N_VScale(ONE, ap, vtemp);
            }

            /*   Apply b-scaling: Ap = sb P1_inv A P2_inv sx_inv p */
            if scale_b {
                N_VProd(sb.unwrap(), vtemp, ap);
            } else {
                N_VScale(ONE, vtemp, ap);
            }

            /* Calculate alpha = <r,r_star>/<Ap,r_star> */
            let mut alpha = N_VDotProd(ap, r_star);
            alpha = beta_denom / alpha;

            /* Update q = r - alpha*Ap = r - alpha*(sb P1_inv A P2_inv sx_inv p) */
            N_VLinearSum(ONE, r, -alpha, ap, q);

            /* Generate u = A-tilde q */

            /*   Apply x-scaling: vtemp = sx_inv q */
            if scale_x {
                N_VDiv(q, sx.unwrap(), vtemp);
            } else {
                N_VScale(ONE, q, vtemp);
            }

            /*   Apply right preconditioner: vtemp = P2_inv sx_inv q */
            if pre_on_right {
                N_VScale(ONE, vtemp, u);
                let status = (psolve.as_mut().unwrap())(u, vtemp, delta, SUN_PREC_RIGHT);
                if status != 0 {
                    *zeroguess = false;
                    *last_flag = if status < 0 {
                        SUNLS_PSOLVE_FAIL_UNREC
                    } else {
                        SUNLS_PSOLVE_FAIL_REC
                    } as i64;
                    return *last_flag as i32;
                }
            }

            /*   Apply A: u = A P2_inv sx_inv u */
            let status = atimes(vtemp, u);
            if status != 0 {
                *zeroguess = false;
                *last_flag = if status < 0 {
                    SUNLS_ATIMES_FAIL_UNREC
                } else {
                    SUNLS_ATIMES_FAIL_REC
                } as i64;
                return *last_flag as i32;
            }

            /*   Apply left preconditioner: vtemp = P1_inv A P2_inv sx_inv p */
            if pre_on_left {
                let status = (psolve.as_mut().unwrap())(u, vtemp, delta, SUN_PREC_LEFT);
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
                N_VScale(ONE, u, vtemp);
            }

            /*   Apply b-scaling: u = sb P1_inv A P2_inv sx_inv u */
            if scale_b {
                N_VProd(sb.unwrap(), vtemp, u);
            } else {
                N_VScale(ONE, vtemp, u);
            }

            /* Calculate omega = <u,q>/<u,u> */
            let mut omega_denom = N_VDotProd(u, u);
            if omega_denom == ZERO {
                omega_denom = ONE;
            }
            let mut omega = N_VDotProd(u, q);
            omega /= omega_denom;

            /* Update x = x + alpha*p + omega*q */
            if l == 0 && *zeroguess {
                N_VLinearSum(alpha, p, omega, q, x);
            } else {
                /* N_VLinearCombination(3, {1, alpha, omega}, {x, p, q}, x) */
                for (zi, xi) in x.data.iter_mut().zip(&p.data) {
                    *zi += alpha * *xi;
                }
                for (zi, xi) in x.data.iter_mut().zip(&q.data) {
                    *zi += omega * *xi;
                }
            }

            /* Update the residual r = q - omega*u */
            N_VLinearSum(ONE, q, -omega, u, r);

            /* Set rho = norm(r) and check convergence */
            rho = SUNRsqrt(N_VDotProd(r, r));
            *resnorm = rho;

            if rho <= delta {
                converged = true;
                break;
            }

            /* Not yet converged, continue iteration */
            /* Update beta = <rnew,r_star> / <rold,r_start> * alpha / omega */
            let beta_num = N_VDotProd(r, r_star);
            let beta = (beta_num / beta_denom) * (alpha / omega);

            /* Update p = r + beta*(p - omega*Ap) = beta*p - beta*omega*Ap + r
               (N_VLinearCombination(3, {beta, -alpha*(beta_num/beta_denom), 1},
                                        {p, Ap, r}, p)) */
            let c1 = -alpha * (beta_num / beta_denom);
            for zi in p.data.iter_mut() {
                *zi *= beta;
            }
            for (zi, xi) in p.data.iter_mut().zip(&ap.data) {
                *zi += c1 * *xi;
            }
            for (zi, xi) in p.data.iter_mut().zip(&r.data) {
                *zi += ONE * *xi;
            }

            /* update beta_denom for next iteration */
            beta_denom = beta_num;
        }

        /* Main loop finished */
        if converged || rho < r_norm {
            /* Apply the x-scaling and right preconditioner: x = P2_inv sx_inv x */
            if scale_x {
                /* N_VDiv(x, sx, x) */
                x.div_with(sx.unwrap());
            }
            if pre_on_right {
                let status = (psolve.as_mut().unwrap())(x, vtemp, delta, SUN_PREC_RIGHT);
                if status != 0 {
                    *zeroguess = false;
                    *last_flag = if status < 0 {
                        SUNLS_PSOLVE_FAIL_UNREC
                    } else {
                        SUNLS_PSOLVE_FAIL_REC
                    } as i64;
                    return *last_flag as i32;
                }
                N_VScale(ONE, vtemp, x);
            }

            *zeroguess = false;
            *last_flag = if converged {
                SUN_SUCCESS as i64
            } else {
                SUNLS_RES_REDUCED as i64
            };
            *last_flag as i32
        } else {
            *zeroguess = false;
            *last_flag = SUNLS_CONV_FAIL as i64;
            *last_flag as i32
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /* nonsymmetric, diagonally dominant tridiagonal test matrix */
    fn matvec(v: &NVector, av: &mut NVector) {
        let n = v.len();
        for i in 0..n {
            let mut s = 4.0 * v.data[i];
            if i > 0 {
                s -= v.data[i - 1];
            }
            if i + 1 < n {
                s += 0.5 * v.data[i + 1];
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
        NVector::from_slice(&(0..n).map(|i| 1.0 + (i % 3) as f64).collect::<Vec<_>>())
    }

    fn unwrap_spbcgs(ls: LinearSolver) -> SpbcgsLS {
        match ls {
            LinearSolver::Spbcgs(s) => s,
            _ => unreachable!(),
        }
    }

    #[test]
    fn spbcgs_solves_unpreconditioned() {
        let n = 12usize;
        let ctx = SUNContext::default();
        let y = NVector::new(n);
        let mut s = unwrap_spbcgs(SUNLinSol_SPBCGS(&y, SUN_PREC_NONE, 200, &ctx));
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
    fn spbcgs_identity_preconditioner_left() {
        let n = 12usize;
        let ctx = SUNContext::default();
        let y = NVector::new(n);
        let mut s = unwrap_spbcgs(SUNLinSol_SPBCGS(&y, SUN_PREC_LEFT, 200, &ctx));
        assert_eq!(s.initialize(), SUN_SUCCESS);

        let b = rhs(n);
        let mut x = NVector::new(n);
        s.zeroguess = true;
        let delta = 1e-10;
        let mut atimes = |v: &NVector, av: &mut NVector| -> i32 {
            matvec(v, av);
            0
        };
        let mut psolve = |r: &NVector, z: &mut NVector, _tol: f64, _lr: i32| -> i32 {
            z.data.copy_from_slice(&r.data);
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
