/* -----------------------------------------------------------------
 * Translated from src/sunlinsol/sptfqmr/sunlinsol_sptfqmr.c
 * (SUNDIALS 7.7.0): scaled preconditioned Transpose-Free QMR.
 *
 * The C content struct becomes SptfqmrLS; the ops-table entry points
 * become methods dispatched from sundials_linearsolver.rs. ATimes /
 * PSolve arrive as closures at solve time; a required-but-missing
 * psolve yields SUNLS_PSOLVE_NULL. As in C, right preconditioning
 * with a nonzero initial guess is an unsupported configuration and
 * returns SUN_ERR_ARG_INCOMPATIBLE. The C two-entry work array
 * `r[2]` becomes the fields r_0 / r_1 (the C code only ever indexes
 * it with the constants 0 and 1). The fused N_VLinearCombination
 * kernel is reproduced inline with the exact floating-point
 * operation order of the serial implementation.
 * -----------------------------------------------------------------*/
use crate::nvector_serial::{
    NVector, N_VClone, N_VConst, N_VDiv, N_VDotProd, N_VLinearSum, N_VProd, N_VScale,
};
use crate::sundials_context::SUNContext;
use crate::sundials_errors::{SUNErrCode, SUN_ERR_ARG_INCOMPATIBLE, SUN_ERR_ARG_OUTOFRANGE, SUN_SUCCESS};
use crate::sundials_linearsolver::{
    ATimesFn, LinearSolver, PSolveFn, SUNLS_ATIMES_FAIL_REC, SUNLS_ATIMES_FAIL_UNREC,
    SUNLS_CONV_FAIL, SUNLS_PSOLVE_FAIL_REC, SUNLS_PSOLVE_FAIL_UNREC, SUNLS_PSOLVE_NULL,
    SUNLS_RES_REDUCED, SUN_PREC_BOTH, SUN_PREC_LEFT, SUN_PREC_NONE, SUN_PREC_RIGHT,
};
use crate::sundials_math::{SUNRsqrt, SUNSQR};

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

/* Default SPTFQMR parameters (sunlinsol_sptfqmr.h) */
pub const SUNSPTFQMR_MAXL_DEFAULT: i32 = 5;

/// SUNLinearSolverContent_SPTFQMR
pub struct SptfqmrLS {
    pub numiters: i32,
    pub resnorm: f64,
    pub last_flag: i64,
    pub zeroguess: bool,
    maxl: i32,
    pretype: i32,
    /* workspace vectors, all allocated by the constructor as in C */
    r_star: NVector,
    q: NVector,
    d: NVector,
    v: NVector,
    p: NVector,
    r_0: NVector, /* C: r[0] */
    r_1: NVector, /* C: r[1] */
    u: NVector,
    pub vtemp1: NVector,
    vtemp2: NVector,
    vtemp3: NVector,
}

/* ----------------------------------------------------------------------------
 * Function to create a new SPTFQMR linear solver (SUNLinSol_SPTFQMR)
 */
pub fn SUNLinSol_SPTFQMR(
    y: &NVector,
    pretype: i32,
    maxl: i32,
    _sunctx: &SUNContext,
) -> LinearSolver {
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
    let maxl = if maxl <= 0 { SUNSPTFQMR_MAXL_DEFAULT } else { maxl };

    LinearSolver::Sptfqmr(SptfqmrLS {
        numiters: 0,
        resnorm: ZERO,
        last_flag: 0,
        zeroguess: false,
        maxl,
        pretype,
        r_star: N_VClone(y),
        q: N_VClone(y),
        d: N_VClone(y),
        v: N_VClone(y),
        p: N_VClone(y),
        r_0: N_VClone(y),
        r_1: N_VClone(y),
        u: N_VClone(y),
        vtemp1: N_VClone(y),
        vtemp2: N_VClone(y),
        vtemp3: N_VClone(y),
    })
}

impl SptfqmrLS {
    /// SUNLinSolSpace_SPTFQMR: (lenrwLS, leniwLS)
    pub fn space(&self) -> (i64, i64) {
        let (lrw1, liw1) = crate::nvector_serial::N_VSpace(&self.r_star);
        (lrw1 * 11, liw1 * 11)
    }

    /// SUNLinSolInitialize_SPTFQMR
    pub fn initialize(&mut self) -> SUNErrCode {
        /* ensure valid options */
        if self.maxl <= 0 {
            self.maxl = SUNSPTFQMR_MAXL_DEFAULT;
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

    /// SUNLinSol_SPTFQMRSetPrecType
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

    /// SUNLinSol_SPTFQMRSetMaxl
    pub fn set_maxl(&mut self, maxl: i32) -> SUNErrCode {
        self.maxl = if maxl <= 0 { SUNSPTFQMR_MAXL_DEFAULT } else { maxl };
        SUN_SUCCESS
    }

    /// SUNLinSolSolve_SPTFQMR
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
        let SptfqmrLS {
            numiters,
            resnorm,
            last_flag,
            zeroguess,
            r_star,
            q,
            d,
            v,
            p,
            r_0,
            r_1,
            u,
            vtemp1,
            vtemp2,
            vtemp3,
            ..
        } = self;
        let sb = s1;
        let sx = s2;

        /* Initialize counters and convergence flag */
        let mut temp_val = -ONE;
        let mut r_curr_norm = -ONE;
        *numiters = 0;
        let mut converged = false;
        let mut b_ok = false;

        let mut rho = [ZERO; 2];

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

        /* Set r_star to initial (unscaled) residual r_star = r_0 = b - A*x_0 */
        /* NOTE: if x == 0 then just set residual to b and continue */
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

        /* Apply left preconditioner and b-scaling to r_star (or really just r_0) */
        if pre_on_left {
            let status = (psolve.as_mut().unwrap())(r_star, vtemp1, delta, SUN_PREC_LEFT);
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
            N_VScale(ONE, r_star, vtemp1);
        }

        if scale_b {
            N_VProd(sb.unwrap(), vtemp1, r_star);
        } else {
            N_VScale(ONE, vtemp1, r_star);
        }

        /* Initialize rho[0] */
        /* NOTE: initialized here to reduce number of computations - avoid need
                 to compute r_star^T*r_star twice, and avoid needlessly squaring
                 values */
        rho[0] = N_VDotProd(r_star, r_star);

        /* Compute norm of initial residual (r_0) to see if we really need
           to do anything */
        let r_init_norm = SUNRsqrt(rho[0]);
        *resnorm = r_init_norm;

        if r_init_norm <= delta {
            *zeroguess = false;
            *last_flag = SUN_SUCCESS as i64;
            return SUN_SUCCESS;
        }

        /* Set v = A*r_0 (preconditioned and scaled) */
        if scale_x {
            N_VDiv(r_star, sx.unwrap(), vtemp1);
        } else {
            N_VScale(ONE, r_star, vtemp1);
        }

        if pre_on_right {
            N_VScale(ONE, vtemp1, v);
            let status = (psolve.as_mut().unwrap())(v, vtemp1, delta, SUN_PREC_RIGHT);
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

        let status = atimes(vtemp1, v);
        if status != 0 {
            *zeroguess = false;
            *last_flag = if status < 0 {
                SUNLS_ATIMES_FAIL_UNREC
            } else {
                SUNLS_ATIMES_FAIL_REC
            } as i64;
            return *last_flag as i32;
        }

        if pre_on_left {
            let status = (psolve.as_mut().unwrap())(v, vtemp1, delta, SUN_PREC_LEFT);
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
            N_VScale(ONE, v, vtemp1);
        }

        if scale_b {
            N_VProd(sb.unwrap(), vtemp1, v);
        } else {
            N_VScale(ONE, vtemp1, v);
        }

        /* Initialize remaining variables */
        N_VScale(ONE, r_star, r_0);
        N_VScale(ONE, r_star, u);
        N_VScale(ONE, r_star, p);
        N_VConst(ZERO, d);

        /* Set x = sx x if non-zero guess */
        if scale_x && !*zeroguess {
            /* N_VProd(sx, x, x) */
            x.prod_with(sx.unwrap());
        }

        let mut tau = r_init_norm;
        let mut v_bar = ZERO;
        let mut eta = ZERO;

        /* START outer loop */
        for n in 0..l_max {
            /* Increment linear iteration counter */
            *numiters += 1;

            /* sigma = r_star^T*v */
            let sigma = N_VDotProd(r_star, v);

            /* alpha = rho[0]/sigma */
            let alpha = rho[0] / sigma;

            /* q = u-alpha*v */
            N_VLinearSum(ONE, u, -alpha, v, q);

            /* r[1] = r[0]-alpha*A*(u+q) */
            N_VLinearSum(ONE, u, ONE, q, r_1);
            if scale_x {
                /* N_VDiv(r[1], sx, r[1]) */
                r_1.div_with(sx.unwrap());
            }

            if pre_on_right {
                N_VScale(ONE, r_1, vtemp1);
                let status = (psolve.as_mut().unwrap())(vtemp1, r_1, delta, SUN_PREC_RIGHT);
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

            let status = atimes(r_1, vtemp1);
            if status != 0 {
                *zeroguess = false;
                *last_flag = if status < 0 {
                    SUNLS_ATIMES_FAIL_UNREC
                } else {
                    SUNLS_ATIMES_FAIL_REC
                } as i64;
                return *last_flag as i32;
            }

            if pre_on_left {
                let status = (psolve.as_mut().unwrap())(vtemp1, r_1, delta, SUN_PREC_LEFT);
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
                N_VScale(ONE, vtemp1, r_1);
            }

            if scale_b {
                N_VProd(sb.unwrap(), r_1, vtemp1);
            } else {
                N_VScale(ONE, r_1, vtemp1);
            }
            N_VLinearSum(ONE, r_0, -alpha, vtemp1, r_1);

            /* START inner loop */
            for m in 0..2 {
                /* d = [*]+(v_bar^2*eta/alpha)*d */
                /* NOTES:
                 *   (1) [*] = u if m == 0, and q if m == 1
                 *   (2) using temp_val reduces the number of required computations
                 *       if the inner loop is executed twice
                 */
                let omega;
                if m == 0 {
                    temp_val = SUNRsqrt(N_VDotProd(r_1, r_1));
                    let omega_dot = N_VDotProd(r_0, r_0);
                    omega = SUNRsqrt(SUNRsqrt(omega_dot) * temp_val);
                    /* N_VLinearSum(ONE, u, v_bar^2*eta/alpha, d, d) */
                    d.linear_sum_with(SUNSQR(v_bar) * eta / alpha, ONE, u);
                } else {
                    omega = temp_val;
                    /* N_VLinearSum(ONE, q, v_bar^2*eta/alpha, d, d) */
                    d.linear_sum_with(SUNSQR(v_bar) * eta / alpha, ONE, q);
                }

                /* v_bar = omega/tau */
                v_bar = omega / tau;

                /* c = (1+v_bar^2)^(-1/2) */
                let c = ONE / SUNRsqrt(ONE + SUNSQR(v_bar));

                /* tau = tau*v_bar*c */
                tau = tau * v_bar * c;

                /* eta = c^2*alpha */
                eta = SUNSQR(c) * alpha;

                /* x = x+eta*d */
                if n == 0 && m == 0 && *zeroguess {
                    N_VScale(eta, d, x);
                } else {
                    x.linear_sum_with(ONE, eta, d);
                }

                /* Check for convergence... */
                /* NOTE: just use approximation to norm of residual, if possible */
                r_curr_norm = tau * SUNRsqrt((m + 1) as f64);
                *resnorm = r_curr_norm;

                /* Exit inner loop if iteration has converged based upon approximation
                   to norm of current residual */
                if r_curr_norm <= delta {
                    converged = true;
                    break;
                }

                /* Decide if actual norm of residual vector should be computed */
                /* NOTES:
                 *   (1) if r_curr_norm > delta, then check if actual residual norm
                 *       is OK (recall we first compute an approximation)
                 *   (2) if r_curr_norm >= r_init_norm and m == 1 and n == l_max, then
                 *       compute actual residual norm to see if the iteration can be
                 *       saved
                 *   (3) the scaled and preconditioned right-hand side of the given
                 *       linear system (denoted by b) is only computed once, and the
                 *       result is stored in vtemp3 so it can be reused - reduces the
                 *       number of psolves if using left preconditioning
                 */
                if (r_curr_norm > delta)
                    || (r_curr_norm >= r_init_norm && m == 1 && n == l_max)
                {
                    /* Compute norm of residual ||b-A*x||_2 (preconditioned and scaled) */
                    if scale_x {
                        N_VDiv(x, sx.unwrap(), vtemp1);
                    } else {
                        N_VScale(ONE, x, vtemp1);
                    }

                    if pre_on_right {
                        let status =
                            (psolve.as_mut().unwrap())(vtemp1, vtemp2, delta, SUN_PREC_RIGHT);
                        if status != 0 {
                            *zeroguess = false;
                            *last_flag = if status < 0 {
                                SUNLS_PSOLVE_FAIL_UNREC
                            } else {
                                SUNLS_PSOLVE_FAIL_REC
                            } as i64;
                            return *last_flag as i32;
                        }
                        N_VScale(ONE, vtemp2, vtemp1);
                    }

                    let status = atimes(vtemp1, vtemp2);
                    if status != 0 {
                        *zeroguess = false;
                        *last_flag = if status < 0 {
                            SUNLS_ATIMES_FAIL_UNREC
                        } else {
                            SUNLS_ATIMES_FAIL_REC
                        } as i64;
                        return *last_flag as i32;
                    }

                    if pre_on_left {
                        let status =
                            (psolve.as_mut().unwrap())(vtemp2, vtemp1, delta, SUN_PREC_LEFT);
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
                        N_VScale(ONE, vtemp2, vtemp1);
                    }

                    if scale_b {
                        N_VProd(sb.unwrap(), vtemp1, vtemp2);
                    } else {
                        N_VScale(ONE, vtemp1, vtemp2);
                    }

                    /* Only precondition and scale b once (result saved for reuse) */
                    if !b_ok {
                        b_ok = true;
                        if pre_on_left {
                            let status =
                                (psolve.as_mut().unwrap())(b, vtemp3, delta, SUN_PREC_LEFT);
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
                            N_VScale(ONE, b, vtemp3);
                        }

                        if scale_b {
                            /* N_VProd(sb, vtemp3, vtemp3) */
                            vtemp3.prod_with(sb.unwrap());
                        }
                    }
                    N_VLinearSum(ONE, vtemp3, -ONE, vtemp2, vtemp1);
                    r_curr_norm = N_VDotProd(vtemp1, vtemp1);
                    r_curr_norm = SUNRsqrt(r_curr_norm);
                    *resnorm = r_curr_norm;

                    /* Exit inner loop if inequality condition is satisfied
                       (meaning exit if we have converged) */
                    if r_curr_norm <= delta {
                        converged = true;
                        break;
                    }
                }
            } /* END inner loop */

            /* If converged, then exit outer loop as well */
            if converged {
                break;
            }

            /* rho[1] = r_star^T*r_[1] */
            rho[1] = N_VDotProd(r_star, r_1);

            /* beta = rho[1]/rho[0] */
            let beta = rho[1] / rho[0];

            /* u = r[1]+beta*q */
            N_VLinearSum(ONE, r_1, beta, q, u);

            /* p = u+beta*(q+beta*p) = beta*beta*p + beta*q + u
               (N_VLinearCombination(3, {beta^2, beta, 1}, {p, q, u}, p)) */
            let c0 = SUNSQR(beta);
            for zi in p.data.iter_mut() {
                *zi *= c0;
            }
            for (zi, xi) in p.data.iter_mut().zip(&q.data) {
                *zi += beta * *xi;
            }
            for (zi, xi) in p.data.iter_mut().zip(&u.data) {
                *zi += ONE * *xi;
            }

            /* v = A*p */
            if scale_x {
                N_VDiv(p, sx.unwrap(), vtemp1);
            } else {
                N_VScale(ONE, p, vtemp1);
            }

            if pre_on_right {
                N_VScale(ONE, vtemp1, v);
                let status = (psolve.as_mut().unwrap())(v, vtemp1, delta, SUN_PREC_RIGHT);
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

            let status = atimes(vtemp1, v);
            if status != 0 {
                *zeroguess = false;
                *last_flag = if status < 0 {
                    SUNLS_ATIMES_FAIL_UNREC
                } else {
                    SUNLS_ATIMES_FAIL_REC
                } as i64;
                return *last_flag as i32;
            }

            if pre_on_left {
                let status = (psolve.as_mut().unwrap())(v, vtemp1, delta, SUN_PREC_LEFT);
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
                N_VScale(ONE, v, vtemp1);
            }

            if scale_b {
                N_VProd(sb.unwrap(), vtemp1, v);
            } else {
                N_VScale(ONE, vtemp1, v);
            }

            /* Shift variable values */
            /* NOTE: reduces storage requirements */
            N_VScale(ONE, r_1, r_0);
            rho[0] = rho[1];
        } /* END outer loop */

        /* Determine return value */
        /* If iteration converged or residual was reduced, then return current
         * iterate (x) */
        if converged || r_curr_norm < r_init_norm {
            if scale_x {
                /* N_VDiv(x, sx, x) */
                x.div_with(sx.unwrap());
            }

            if pre_on_right {
                let status = (psolve.as_mut().unwrap())(x, vtemp1, delta, SUN_PREC_RIGHT);
                if status != 0 {
                    *zeroguess = false;
                    *last_flag = if status < 0 {
                        SUNLS_PSOLVE_FAIL_UNREC
                    } else {
                        SUNLS_PSOLVE_FAIL_REC
                    } as i64;
                    return *last_flag as i32;
                }
                N_VScale(ONE, vtemp1, x);
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
        NVector::from_slice(&(0..n).map(|i| 2.0 - (i % 4) as f64 * 0.25).collect::<Vec<_>>())
    }

    fn unwrap_sptfqmr(ls: LinearSolver) -> SptfqmrLS {
        match ls {
            LinearSolver::Sptfqmr(s) => s,
            _ => unreachable!(),
        }
    }

    #[test]
    fn sptfqmr_solves_unpreconditioned() {
        let n = 12usize;
        let ctx = SUNContext::default();
        let y = NVector::new(n);
        let mut s = unwrap_sptfqmr(SUNLinSol_SPTFQMR(&y, SUN_PREC_NONE, 200, &ctx));
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
        assert!(residual_norm(&x, &b) <= 1e-8);
    }

    #[test]
    fn sptfqmr_identity_preconditioner_left() {
        let n = 12usize;
        let ctx = SUNContext::default();
        let y = NVector::new(n);
        let mut s = unwrap_sptfqmr(SUNLinSol_SPTFQMR(&y, SUN_PREC_LEFT, 200, &ctx));
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
