/* -----------------------------------------------------------------
 * Translated from src/sunlinsol/spgmr/sunlinsol_spgmr.c
 * (SUNDIALS 7.7.0): scaled preconditioned GMRES.
 *
 * The C content struct becomes SpgmrLS; the ops-table entry points
 * become methods dispatched from sundials_linearsolver.rs. ATimes /
 * PSolve arrive as closures at solve time instead of being stored,
 * so the C SetATimes/SetPreconditioner/SetScalingVectors plumbing
 * (and the corresponding NULL-pointer asserts in Initialize) is
 * carried by the solve() arguments; a required-but-missing psolve
 * yields SUNLS_PSOLVE_NULL. The command-line option plumbing
 * (SUNLinSolSetOptions_SPGMR) is not applicable. The fused
 * N_VLinearCombination kernels are reproduced inline with the exact
 * floating-point operation order of the serial implementation.
 * -----------------------------------------------------------------*/
use crate::nvector_serial::{NVector, N_VClone, N_VConst, N_VDiv, N_VDotProd, N_VProd, N_VScale};
use crate::sundials_context::SUNContext;
use crate::sundials_errors::{SUNErrCode, SUN_ERR_ARG_OUTOFRANGE, SUN_SUCCESS};
use crate::sundials_iterative::{SUNClassicalGS, SUNModifiedGS, SUNQRfact, SUNQRsol};
use crate::sundials_linearsolver::{
    ATimesFn, LinearSolver, PSolveFn, SUNLS_ATIMES_FAIL_REC, SUNLS_ATIMES_FAIL_UNREC,
    SUNLS_CONV_FAIL, SUNLS_PSOLVE_FAIL_REC, SUNLS_PSOLVE_FAIL_UNREC, SUNLS_PSOLVE_NULL,
    SUNLS_QRFACT_FAIL, SUNLS_QRSOL_FAIL, SUNLS_RES_REDUCED, SUN_CLASSICAL_GS, SUN_MODIFIED_GS,
    SUN_PREC_BOTH, SUN_PREC_LEFT, SUN_PREC_NONE, SUN_PREC_RIGHT,
};
use crate::sundials_math::{SUNRabs, SUNRsqrt};

const ZERO: f64 = 0.0;
const ONE: f64 = 1.0;

/* Default SPGMR parameters (sunlinsol_spgmr.h) */
pub const SUNSPGMR_MAXL_DEFAULT: i32 = 5;
pub const SUNSPGMR_MAXRS_DEFAULT: i32 = 0;
pub const SUNSPGMR_GSTYPE_DEFAULT: i32 = SUN_MODIFIED_GS;

/// SUNLinearSolverContent_SPGMR
pub struct SpgmrLS {
    pub numiters: i32,
    pub resnorm: f64,
    pub last_flag: i64,
    pub zeroguess: bool,
    maxl: i32,
    pretype: i32,
    gstype: i32,
    max_restarts: i32,
    /* workspace: xcor and vtemp are allocated by the constructor,
       V/Hes/givens/yg/cv by initialize(), mirroring where the C code
       allocates them (SUNLinSol_SPGMR vs SUNLinSolInitialize_SPGMR). */
    v: Vec<NVector>,      /* Krylov basis V[0..maxl] */
    hes: Vec<Vec<f64>>,   /* (maxl+1) x maxl Hessenberg matrix */
    givens: Vec<f64>,     /* 2*maxl Givens rotation components */
    yg: Vec<f64>,         /* length maxl+1 */
    cv: Vec<f64>,         /* length maxl+1 scalar workspace (classical GS) */
    xcor: NVector,
    vtemp: NVector,
}

/* ----------------------------------------------------------------------------
 * Function to create a new SPGMR linear solver (SUNLinSol_SPGMR)
 */
pub fn SUNLinSol_SPGMR(y: &NVector, pretype: i32, maxl: i32, _sunctx: &SUNContext) -> LinearSolver {
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
    let maxl = if maxl <= 0 { SUNSPGMR_MAXL_DEFAULT } else { maxl };

    LinearSolver::Spgmr(SpgmrLS {
        numiters: 0,
        resnorm: ZERO,
        last_flag: 0,
        zeroguess: false,
        maxl,
        pretype,
        gstype: SUNSPGMR_GSTYPE_DEFAULT,
        max_restarts: SUNSPGMR_MAXRS_DEFAULT,
        v: Vec::new(),
        hes: Vec::new(),
        givens: Vec::new(),
        yg: Vec::new(),
        cv: Vec::new(),
        xcor: N_VClone(y),
        vtemp: N_VClone(y),
    })
}

impl SpgmrLS {
    /// SUNLinSolSpace_SPGMR: (lenrwLS, leniwLS)
    pub fn space(&self) -> (i64, i64) {
        let maxl = self.maxl as i64;
        let (lrw1, liw1) = crate::nvector_serial::N_VSpace(&self.vtemp);
        (lrw1 * (maxl + 5) + maxl * (maxl + 5) + 2, liw1 * (maxl + 5))
    }

    /// SUNLinSolInitialize_SPGMR
    pub fn initialize(&mut self) -> SUNErrCode {
        /* ensure valid options */
        if self.max_restarts < 0 {
            self.max_restarts = SUNSPGMR_MAXRS_DEFAULT;
        }

        if (self.pretype != SUN_PREC_LEFT)
            && (self.pretype != SUN_PREC_RIGHT)
            && (self.pretype != SUN_PREC_BOTH)
        {
            self.pretype = SUN_PREC_NONE;
        }

        /* allocate solver-specific memory (where the size depends on the
           choice of maxl) here */
        let maxl = self.maxl as usize;
        if self.v.is_empty() {
            self.v = (0..=maxl).map(|_| N_VClone(&self.vtemp)).collect();
        }
        if self.hes.is_empty() {
            self.hes = vec![vec![0.0; maxl]; maxl + 1];
        }
        if self.givens.is_empty() {
            self.givens = vec![0.0; 2 * maxl];
        }
        if self.yg.is_empty() {
            self.yg = vec![0.0; maxl + 1];
        }
        if self.cv.is_empty() {
            self.cv = vec![0.0; maxl + 1];
        }

        SUN_SUCCESS
    }

    /// SUNLinSol_SPGMRSetPrecType
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

    /// SUNLinSol_SPGMRSetGSType
    pub fn set_gs_type(&mut self, gstype: i32) -> SUNErrCode {
        /* Check for legal gstype */
        if gstype != SUN_MODIFIED_GS && gstype != SUN_CLASSICAL_GS {
            return SUN_ERR_ARG_OUTOFRANGE;
        }
        self.gstype = gstype;
        SUN_SUCCESS
    }

    /// SUNLinSol_SPGMRSetMaxRestarts
    pub fn set_max_restarts(&mut self, maxrs: i32) -> SUNErrCode {
        /* Illegal maxrs implies use of default value */
        self.max_restarts = if maxrs < 0 { SUNSPGMR_MAXRS_DEFAULT } else { maxrs };
        SUN_SUCCESS
    }

    /// SUNLinSolSolve_SPGMR
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
        /* Initialize some variables */
        let mut krydim: usize = 0;

        /* Make local shortcuts to solver variables. */
        let l_max = self.maxl as usize;
        let max_restarts = self.max_restarts;
        let gstype = self.gstype;
        let pretype = self.pretype;
        let SpgmrLS {
            numiters,
            resnorm,
            last_flag,
            zeroguess,
            v,
            hes,
            givens,
            yg,
            cv,
            xcor,
            vtemp,
            ..
        } = self;

        /* Initialize counters and convergence flag */
        *numiters = 0;
        let mut converged = false;

        /* Set flags for internal solver options */
        let pre_on_left = (pretype == SUN_PREC_LEFT) || (pretype == SUN_PREC_BOTH);
        let pre_on_right = (pretype == SUN_PREC_RIGHT) || (pretype == SUN_PREC_BOTH);
        let scale1 = s1.is_some();
        let scale2 = s2.is_some();

        /* If preconditioning, check if psolve has been set */
        if (pre_on_left || pre_on_right) && psolve.is_none() {
            *last_flag = SUNLS_PSOLVE_NULL as i64;
            return SUNLS_PSOLVE_NULL;
        }

        /* Set vtemp and V[0] to initial (unscaled) residual r_0 = b - A*x_0 */
        if *zeroguess {
            N_VScale(ONE, b, vtemp);
        } else {
            let status = atimes(x, vtemp);
            if status != 0 {
                *zeroguess = false;
                *last_flag = if status < 0 {
                    SUNLS_ATIMES_FAIL_UNREC
                } else {
                    SUNLS_ATIMES_FAIL_REC
                } as i64;
                return *last_flag as i32;
            }
            /* vtemp = b - vtemp */
            vtemp.linear_sum_with(-ONE, ONE, b);
        }
        N_VScale(ONE, vtemp, &mut v[0]);

        /* Apply left preconditioner and left scaling to V[0] = r_0 */
        if pre_on_left {
            let status = (psolve.as_mut().unwrap())(&v[0], vtemp, delta, SUN_PREC_LEFT);
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
            N_VScale(ONE, &v[0], vtemp);
        }

        if scale1 {
            N_VProd(s1.unwrap(), vtemp, &mut v[0]);
        } else {
            N_VScale(ONE, vtemp, &mut v[0]);
        }

        /* Set r_norm = beta to L2 norm of V[0] = s1 P1_inv r_0, and
           return if small */
        let mut r_norm = N_VDotProd(&v[0], &v[0]);
        r_norm = SUNRsqrt(r_norm);
        let beta = r_norm;
        *resnorm = r_norm;

        if r_norm <= delta {
            *zeroguess = false;
            *last_flag = SUN_SUCCESS as i64;
            return SUN_SUCCESS;
        }

        /* Initialize rho to avoid compiler warning message */
        let mut rho = beta;

        /* Set xcor = 0 */
        N_VConst(ZERO, xcor);

        /* Begin outer iterations: up to (max_restarts + 1) attempts */
        for ntries in 0..=max_restarts {
            /* Initialize the Hessenberg matrix Hes and Givens rotation
               product.  Normalize the initial vector V[0] */
            for i in 0..=l_max {
                for j in 0..l_max {
                    hes[i][j] = ZERO;
                }
            }

            let mut rotation_product = ONE;
            v[0].scale_inplace(ONE / r_norm);

            /* Inner loop: generate Krylov sequence and Arnoldi basis */
            for l in 0..l_max {
                *numiters += 1;
                let l_plus_1 = l + 1;
                krydim = l_plus_1;

                /* Generate A-tilde V[l], where A-tilde = s1 P1_inv A P2_inv s2_inv */

                /*   Apply right scaling: vtemp = s2_inv V[l] */
                if scale2 {
                    N_VDiv(&v[l], s2.unwrap(), vtemp);
                } else {
                    N_VScale(ONE, &v[l], vtemp);
                }

                /*   Apply right preconditioner: vtemp = P2_inv s2_inv V[l] */
                if pre_on_right {
                    N_VScale(ONE, vtemp, &mut v[l_plus_1]);
                    let status =
                        (psolve.as_mut().unwrap())(&v[l_plus_1], vtemp, delta, SUN_PREC_RIGHT);
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

                /* Apply A: V[l+1] = A P2_inv s2_inv V[l] */
                let status = atimes(vtemp, &mut v[l_plus_1]);
                if status != 0 {
                    *zeroguess = false;
                    *last_flag = if status < 0 {
                        SUNLS_ATIMES_FAIL_UNREC
                    } else {
                        SUNLS_ATIMES_FAIL_REC
                    } as i64;
                    return *last_flag as i32;
                }

                /* Apply left preconditioning: vtemp = P1_inv A P2_inv s2_inv V[l] */
                if pre_on_left {
                    let status =
                        (psolve.as_mut().unwrap())(&v[l_plus_1], vtemp, delta, SUN_PREC_LEFT);
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
                    N_VScale(ONE, &v[l_plus_1], vtemp);
                }

                /* Apply left scaling: V[l+1] = s1 P1_inv A P2_inv s2_inv V[l] */
                if scale1 {
                    N_VProd(s1.unwrap(), vtemp, &mut v[l_plus_1]);
                } else {
                    N_VScale(ONE, vtemp, &mut v[l_plus_1]);
                }

                /*  Orthogonalize V[l+1] against previous V[i]: V[l+1] = w_tilde */
                let mut new_vk_norm = ZERO;
                if gstype == SUN_CLASSICAL_GS {
                    let _ =
                        SUNClassicalGS(v, hes, l_plus_1 as i32, l_max as i32, &mut new_vk_norm, cv);
                } else {
                    let _ = SUNModifiedGS(v, hes, l_plus_1 as i32, l_max as i32, &mut new_vk_norm);
                }
                hes[l_plus_1][l] = new_vk_norm;

                /*  Update the QR factorization of Hes */
                if SUNQRfact(krydim as i32, hes, givens, l as i32) != 0 {
                    *zeroguess = false;
                    *last_flag = SUNLS_QRFACT_FAIL as i64;
                    return *last_flag as i32;
                }

                /*  Update residual norm estimate; break if convergence test passes */
                rotation_product *= givens[2 * l + 1];
                rho = SUNRabs(rotation_product * r_norm);
                *resnorm = rho;

                if rho <= delta {
                    converged = true;
                    break;
                }

                /* Normalize V[l+1] with norm value from the Gram-Schmidt routine */
                v[l_plus_1].scale_inplace(ONE / hes[l_plus_1][l]);
            }

            /* Inner loop is done.  Compute the new correction vector xcor */

            /*   Construct g, then solve for y */
            yg[0] = r_norm;
            for i in 1..=krydim {
                yg[i] = ZERO;
            }
            if SUNQRsol(krydim as i32, hes, givens, yg) != 0 {
                *zeroguess = false;
                *last_flag = SUNLS_QRSOL_FAIL as i64;
                return *last_flag as i32;
            }

            /*   Add correction vector V_l y to xcor
               (N_VLinearCombination with cv[0] = 1 and Xv[0] = xcor) */
            for k in 0..krydim {
                let c = yg[k];
                for (zi, xi) in xcor.data.iter_mut().zip(&v[k].data) {
                    *zi += c * *xi;
                }
            }

            /* If converged, construct the final solution vector x and return */
            if converged {
                /* Apply right scaling and right precond.: vtemp = P2_inv s2_inv xcor */
                if scale2 {
                    xcor.div_with(s2.unwrap());
                }

                if pre_on_right {
                    let status = (psolve.as_mut().unwrap())(xcor, vtemp, delta, SUN_PREC_RIGHT);
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
                    N_VScale(ONE, xcor, vtemp);
                }

                /* Add vtemp to initial x to get final solution x, and return */
                if *zeroguess {
                    N_VScale(ONE, vtemp, x);
                } else {
                    x.linear_sum_with(ONE, ONE, vtemp);
                }

                *zeroguess = false;
                *last_flag = SUN_SUCCESS as i64;
                return SUN_SUCCESS;
            }

            /* Not yet converged; if allowed, prepare for restart */
            if ntries == max_restarts {
                break;
            }

            /* Construct last column of Q in yg */
            let mut s_product = ONE;
            for i in (1..=krydim).rev() {
                yg[i] = s_product * givens[2 * i - 2];
                s_product *= givens[2 * i - 1];
            }
            yg[0] = s_product;

            /* Scale r_norm and yg */
            r_norm *= s_product;
            for i in 0..=krydim {
                yg[i] *= r_norm;
            }
            r_norm = SUNRabs(r_norm);

            /* Multiply yg by V_(krydim+1) to get last residual vector; restart
               (N_VLinearCombination with Xv[0] = V[0] and cv[0] = yg[0]) */
            {
                let (v0, vrest) = v.split_at_mut(1);
                let v0 = &mut v0[0];
                for zi in v0.data.iter_mut() {
                    *zi *= yg[0];
                }
                for k in 1..=krydim {
                    let c = yg[k];
                    for (zi, xi) in v0.data.iter_mut().zip(&vrest[k - 1].data) {
                        *zi += c * *xi;
                    }
                }
            }
        }

        /* Failed to converge, even after allowed restarts.
           If the residual norm was reduced below its initial value, compute
           and return x anyway.  Otherwise return failure flag. */
        if rho < beta {
            /* Apply right scaling and right precond.: vtemp = P2_inv s2_inv xcor */
            if scale2 {
                xcor.div_with(s2.unwrap());
            }

            if pre_on_right {
                let status = (psolve.as_mut().unwrap())(xcor, vtemp, delta, SUN_PREC_RIGHT);
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
                N_VScale(ONE, xcor, vtemp);
            }

            /* Add vtemp to initial x to get final solution x, and return */
            if *zeroguess {
                N_VScale(ONE, vtemp, x);
            } else {
                x.linear_sum_with(ONE, ONE, vtemp);
            }

            *zeroguess = false;
            *last_flag = SUNLS_RES_REDUCED as i64;
            return *last_flag as i32;
        }

        *zeroguess = false;
        *last_flag = SUNLS_CONV_FAIL as i64;
        *last_flag as i32
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
        NVector::from_slice(&(0..n).map(|i| (i as f64).sin() + 1.5).collect::<Vec<_>>())
    }

    fn unwrap_spgmr(ls: LinearSolver) -> SpgmrLS {
        match ls {
            LinearSolver::Spgmr(s) => s,
            _ => unreachable!(),
        }
    }

    #[test]
    fn spgmr_solves_unpreconditioned() {
        let n = 12usize;
        let ctx = SUNContext::default();
        let y = NVector::new(n);
        let mut s = unwrap_spgmr(SUNLinSol_SPGMR(&y, SUN_PREC_NONE, n as i32, &ctx));
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
        assert_eq!(s.last_flag, SUN_SUCCESS as i64);
        assert!(s.numiters > 0);
        assert!(s.resnorm <= delta);
        assert!(residual_norm(&x, &b) <= 1e-8);
    }

    #[test]
    fn spgmr_identity_preconditioner_left() {
        let n = 12usize;
        let ctx = SUNContext::default();
        let y = NVector::new(n);
        let mut s = unwrap_spgmr(SUNLinSol_SPGMR(&y, SUN_PREC_LEFT, n as i32, &ctx));
        assert_eq!(s.set_gs_type(SUN_CLASSICAL_GS), SUN_SUCCESS);
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
