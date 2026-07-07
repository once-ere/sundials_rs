/* -----------------------------------------------------------------
 * Translated from src/sunlinsol/spfgmr/sunlinsol_spfgmr.c
 * (SUNDIALS 7.7.0): scaled preconditioned Flexible GMRES.
 *
 * The C content struct becomes SpfgmrLS; the ops-table entry points
 * become methods dispatched from sundials_linearsolver.rs. ATimes /
 * PSolve arrive as closures at solve time; a required-but-missing
 * psolve yields SUNLS_PSOLVE_NULL. SPFGMR supports only right
 * preconditioning: any of SUN_PREC_LEFT/RIGHT/BOTH enables it. The
 * fused N_VLinearCombination kernels are reproduced inline with the
 * exact floating-point operation order of the serial implementation.
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

/* Default SPFGMR parameters (sunlinsol_spfgmr.h) */
pub const SUNSPFGMR_MAXL_DEFAULT: i32 = 5;
pub const SUNSPFGMR_MAXRS_DEFAULT: i32 = 0;
pub const SUNSPFGMR_GSTYPE_DEFAULT: i32 = SUN_MODIFIED_GS;

/// SUNLinearSolverContent_SPFGMR
pub struct SpfgmrLS {
    pub numiters: i32,
    pub resnorm: f64,
    pub last_flag: i64,
    pub zeroguess: bool,
    maxl: i32,
    pretype: i32,
    gstype: i32,
    max_restarts: i32,
    /* workspace: xcor and vtemp are allocated by the constructor,
       V/Z/Hes/givens/yg/cv by initialize(), mirroring where the C code
       allocates them (SUNLinSol_SPFGMR vs SUNLinSolInitialize_SPFGMR). */
    v: Vec<NVector>,      /* Krylov basis V[0..maxl] */
    z: Vec<NVector>,      /* Preconditioned basis Z[0..maxl] */
    hes: Vec<Vec<f64>>,   /* (maxl+1) x maxl Hessenberg matrix */
    givens: Vec<f64>,     /* 2*maxl Givens rotation components */
    yg: Vec<f64>,         /* length maxl+1 */
    cv: Vec<f64>,         /* length maxl+1 scalar workspace (classical GS) */
    xcor: NVector,
    pub vtemp: NVector,
}

/* ----------------------------------------------------------------------------
 * Function to create a new SPFGMR linear solver (SUNLinSol_SPFGMR)
 */
pub fn SUNLinSol_SPFGMR(y: &NVector, pretype: i32, maxl: i32, _sunctx: &SUNContext) -> LinearSolver {
    /* set preconditioning flag (enabling any preconditioner implies right
       preconditioning, since SPFGMR does not support left preconditioning) */
    let pretype = if (pretype == SUN_PREC_LEFT)
        || (pretype == SUN_PREC_RIGHT)
        || (pretype == SUN_PREC_BOTH)
    {
        SUN_PREC_RIGHT
    } else {
        SUN_PREC_NONE
    };

    /* if maxl input is illegal, set to default */
    let maxl = if maxl <= 0 { SUNSPFGMR_MAXL_DEFAULT } else { maxl };

    LinearSolver::Spfgmr(SpfgmrLS {
        numiters: 0,
        resnorm: ZERO,
        last_flag: 0,
        zeroguess: false,
        maxl,
        pretype,
        gstype: SUNSPFGMR_GSTYPE_DEFAULT,
        max_restarts: SUNSPFGMR_MAXRS_DEFAULT,
        v: Vec::new(),
        z: Vec::new(),
        hes: Vec::new(),
        givens: Vec::new(),
        yg: Vec::new(),
        cv: Vec::new(),
        xcor: N_VClone(y),
        vtemp: N_VClone(y),
    })
}

impl SpfgmrLS {
    /// SUNLinSolSpace_SPFGMR: (lenrwLS, leniwLS)
    pub fn space(&self) -> (i64, i64) {
        let maxl = self.maxl as i64;
        let (lrw1, liw1) = crate::nvector_serial::N_VSpace(&self.vtemp);
        (
            lrw1 * (2 * maxl + 4) + maxl * (maxl + 5) + 2,
            liw1 * (2 * maxl + 4),
        )
    }

    /// SUNLinSolInitialize_SPFGMR
    pub fn initialize(&mut self) -> SUNErrCode {
        /* ensure valid options */
        if self.max_restarts < 0 {
            self.max_restarts = SUNSPFGMR_MAXRS_DEFAULT;
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
        if self.z.is_empty() {
            self.z = (0..=maxl).map(|_| N_VClone(&self.vtemp)).collect();
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

    /// SUNLinSol_SPFGMRSetPrecType — toggles preconditioning on/off; turns
    /// on (as right preconditioning) if pretype is any one of SUN_PREC_LEFT,
    /// SUN_PREC_RIGHT or SUN_PREC_BOTH; otherwise turns off.
    pub fn set_prec_type(&mut self, pretype: i32) -> SUNErrCode {
        self.pretype = if (pretype == SUN_PREC_LEFT)
            || (pretype == SUN_PREC_RIGHT)
            || (pretype == SUN_PREC_BOTH)
        {
            SUN_PREC_RIGHT
        } else {
            SUN_PREC_NONE
        };
        SUN_SUCCESS
    }

    /// SUNLinSol_SPFGMRSetGSType
    pub fn set_gs_type(&mut self, gstype: i32) -> SUNErrCode {
        /* Check for legal gstype */
        if gstype != SUN_MODIFIED_GS && gstype != SUN_CLASSICAL_GS {
            return SUN_ERR_ARG_OUTOFRANGE;
        }
        self.gstype = gstype;
        SUN_SUCCESS
    }

    /// SUNLinSol_SPFGMRSetMaxRestarts
    pub fn set_max_restarts(&mut self, maxrs: i32) -> SUNErrCode {
        /* Illegal maxrs implies use of default value */
        self.max_restarts = if maxrs < 0 { SUNSPFGMR_MAXRS_DEFAULT } else { maxrs };
        SUN_SUCCESS
    }

    /// SUNLinSolSolve_SPFGMR
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
        let SpfgmrLS {
            numiters,
            resnorm,
            last_flag,
            zeroguess,
            v,
            z,
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
        let pre_on_right = (pretype == SUN_PREC_LEFT)
            || (pretype == SUN_PREC_RIGHT)
            || (pretype == SUN_PREC_BOTH);
        let scale1 = s1.is_some();
        let scale2 = s2.is_some();

        /* If preconditioning, check if psolve has been set */
        if pre_on_right && psolve.is_none() {
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

        /* Apply left scaling to vtemp = r_0 to fill V[0]. */
        if scale1 {
            N_VProd(s1.unwrap(), vtemp, &mut v[0]);
        } else {
            N_VScale(ONE, vtemp, &mut v[0]);
        }

        /* Set r_norm = beta to L2 norm of V[0] = s1 r_0, and return if small */
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

        /* Set xcor = 0. */
        N_VConst(ZERO, xcor);

        /* Begin outer iterations: up to (max_restarts + 1) attempts. */
        for ntries in 0..=max_restarts {
            /* Initialize the Hessenberg matrix Hes and Givens rotation
               product.  Normalize the initial vector V[0].             */
            for i in 0..=l_max {
                for j in 0..l_max {
                    hes[i][j] = ZERO;
                }
            }
            let mut rotation_product = ONE;
            v[0].scale_inplace(ONE / r_norm);

            /* Inner loop: generate Krylov sequence and Arnoldi basis. */
            for l in 0..l_max {
                *numiters += 1;

                krydim = l + 1;

                /* Generate A-tilde V[l], where A-tilde = s1 A P_inv s2_inv. */

                /*   Apply right scaling: vtemp = s2_inv V[l]. */
                if scale2 {
                    N_VDiv(&v[l], s2.unwrap(), vtemp);
                } else {
                    N_VScale(ONE, &v[l], vtemp);
                }

                /*   Apply right preconditioner: vtemp = Z[l] = P_inv s2_inv V[l]. */
                if pre_on_right {
                    N_VScale(ONE, vtemp, &mut v[l + 1]);
                    let status = (psolve.as_mut().unwrap())(&v[l + 1], vtemp, delta, SUN_PREC_RIGHT);
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
                N_VScale(ONE, vtemp, &mut z[l]);

                /*   Apply A: V[l+1] = A P_inv s2_inv V[l]. */
                let status = atimes(vtemp, &mut v[l + 1]);
                if status != 0 {
                    *zeroguess = false;
                    *last_flag = if status < 0 {
                        SUNLS_ATIMES_FAIL_UNREC
                    } else {
                        SUNLS_ATIMES_FAIL_REC
                    } as i64;
                    return *last_flag as i32;
                }

                /*   Apply left scaling: V[l+1] = s1 A P_inv s2_inv V[l]. */
                if scale1 {
                    /* N_VProd(s1, V[l+1], V[l+1]) */
                    v[l + 1].prod_with(s1.unwrap());
                }

                /* Orthogonalize V[l+1] against previous V[i]: V[l+1] = w_tilde. */
                let mut new_vk_norm = ZERO;
                if gstype == SUN_CLASSICAL_GS {
                    let _ =
                        SUNClassicalGS(v, hes, (l + 1) as i32, l_max as i32, &mut new_vk_norm, cv);
                } else {
                    let _ = SUNModifiedGS(v, hes, (l + 1) as i32, l_max as i32, &mut new_vk_norm);
                }
                hes[l + 1][l] = new_vk_norm;

                /* Update the QR factorization of Hes. */
                if SUNQRfact(krydim as i32, hes, givens, l as i32) != 0 {
                    *zeroguess = false;
                    *last_flag = SUNLS_QRFACT_FAIL as i64;
                    return *last_flag as i32;
                }

                /* Update residual norm estimate; break if convergence test passes. */
                rotation_product *= givens[2 * l + 1];
                rho = SUNRabs(rotation_product * r_norm);
                *resnorm = rho;

                if rho <= delta {
                    converged = true;
                    break;
                }

                /* Normalize V[l+1] with norm value from the Gram-Schmidt routine. */
                v[l + 1].scale_inplace(ONE / hes[l + 1][l]);
            }

            /* Inner loop is done.  Compute the new correction vector xcor. */

            /*   Construct g, then solve for y. */
            yg[0] = r_norm;
            for i in 1..=krydim {
                yg[i] = ZERO;
            }
            if SUNQRsol(krydim as i32, hes, givens, yg) != 0 {
                *zeroguess = false;
                *last_flag = SUNLS_QRSOL_FAIL as i64;
                return *last_flag as i32;
            }

            /*   Add correction vector Z_l y to xcor
               (N_VLinearCombination with cv[0] = 1 and Xv[0] = xcor). */
            for k in 0..krydim {
                let c = yg[k];
                for (zi, xi) in xcor.data.iter_mut().zip(&z[k].data) {
                    *zi += c * *xi;
                }
            }

            /* If converged, construct the final solution vector x and return. */
            if converged {
                if *zeroguess {
                    N_VScale(ONE, xcor, x);
                } else {
                    x.linear_sum_with(ONE, ONE, xcor);
                }
                *zeroguess = false;
                *last_flag = SUN_SUCCESS as i64;
                return SUN_SUCCESS;
            }

            /* Not yet converged; if allowed, prepare for restart. */
            if ntries == max_restarts {
                break;
            }

            /* Construct last column of Q in yg. */
            let mut s_product = ONE;
            for i in (1..=krydim).rev() {
                yg[i] = s_product * givens[2 * i - 2];
                s_product *= givens[2 * i - 1];
            }
            yg[0] = s_product;

            /* Scale r_norm and yg. */
            r_norm *= s_product;
            for i in 0..=krydim {
                yg[i] *= r_norm;
            }
            r_norm = SUNRabs(r_norm);

            /* Multiply yg by V_(krydim+1) to get last residual vector; restart
               (N_VLinearCombination with Xv[0] = V[0] and cv[0] = yg[0]). */
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
            if *zeroguess {
                N_VScale(ONE, xcor, x);
            } else {
                x.linear_sum_with(ONE, ONE, xcor);
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
        NVector::from_slice(&(0..n).map(|i| (i as f64).cos() + 2.0).collect::<Vec<_>>())
    }

    fn unwrap_spfgmr(ls: LinearSolver) -> SpfgmrLS {
        match ls {
            LinearSolver::Spfgmr(s) => s,
            _ => unreachable!(),
        }
    }

    #[test]
    fn spfgmr_solves_unpreconditioned() {
        let n = 12usize;
        let ctx = SUNContext::default();
        let y = NVector::new(n);
        let mut s = unwrap_spfgmr(SUNLinSol_SPFGMR(&y, SUN_PREC_NONE, n as i32, &ctx));
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
    fn spfgmr_identity_preconditioner_classical_gs() {
        let n = 12usize;
        let ctx = SUNContext::default();
        let y = NVector::new(n);
        /* SUN_PREC_LEFT is remapped to right preconditioning by SPFGMR */
        let mut s = unwrap_spfgmr(SUNLinSol_SPFGMR(&y, SUN_PREC_LEFT, n as i32, &ctx));
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
