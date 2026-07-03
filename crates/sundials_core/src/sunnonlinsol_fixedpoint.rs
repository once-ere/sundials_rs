/* -----------------------------------------------------------------
 * Translated from src/sunnonlinsol/fixedpoint/sunnonlinsol_fixedpoint.c
 * (SUNDIALS 7.7.0). Accelerated fixed-point (Anderson acceleration)
 * nonlinear solver state + the AndersonAccelerate kernel. The outer
 * solve loop (SUNNonlinSolSolve_FixedPoint) is driven from
 * cvode_nls.rs with the CVODE Sys/CTest callbacks inlined.
 * -----------------------------------------------------------------*/
use crate::nvector_serial::*;
use crate::sundials_context::SUNContext;
use crate::sundials_math::SUNRsqrt;
use crate::sundials_nonlinearsolver::NonlinearSolver;

pub struct FixedPointSolver {
    /// Anderson acceleration depth (m == 0: basic fixed point)
    pub m: i32,
    pub damping: bool,
    pub beta: f64,
    pub curiter: i32,
    pub maxiters: i32,
    pub niters: i64,
    pub nconvfails: i64,
    /* workspace (AllocateContent) */
    pub yprev: NVector,
    pub gy: NVector,
    pub delta: NVector, // fv workspace in AndersonAccelerate
    pub fold: NVector,
    pub gold: NVector,
    pub imap: Vec<usize>,
    pub R: Vec<f64>,     // m x m, row-major R[i*m + j] as in C
    pub gamma: Vec<f64>, // least-squares coefficients
    pub df: Vec<NVector>,
    pub dg: Vec<NVector>,
    pub q: Vec<NVector>,
}

/// SUNNonlinSol_FixedPoint(y, m)
pub fn SUNNonlinSol_FixedPoint(y: &NVector, m: i32, _sunctx: &SUNContext) -> NonlinearSolver {
    let n = y.len();
    let mm = if m > 0 { m as usize } else { 0 };
    NonlinearSolver::FixedPoint(FixedPointSolver {
        m: if m > 0 { m } else { 0 },
        damping: false,
        beta: 1.0,
        curiter: 0,
        maxiters: 3,
        niters: 0,
        nconvfails: 0,
        yprev: NVector::new(n),
        gy: NVector::new(n),
        delta: NVector::new(n),
        fold: NVector::new(if mm > 0 { n } else { 0 }),
        gold: NVector::new(if mm > 0 { n } else { 0 }),
        imap: vec![0; mm],
        R: vec![0.0; mm * mm],
        gamma: vec![0.0; mm],
        df: (0..mm).map(|_| NVector::new(n)).collect(),
        dg: (0..mm).map(|_| NVector::new(n)).collect(),
        q: (0..mm).map(|_| NVector::new(n)).collect(),
    })
}

impl FixedPointSolver {
    /// SUNNonlinSolSetDamping_FixedPoint
    pub fn set_damping(&mut self, beta: f64) -> i32 {
        if beta <= 0.0 {
            return crate::sundials_errors::SUN_ERR_ARG_OUTOFRANGE;
        }
        if beta < 1.0 {
            self.beta = beta;
            self.damping = true;
        } else {
            self.beta = 1.0;
            self.damping = false;
        }
        0
    }

    /// AndersonAccelerate: computes the Anderson-accelerated fixed-point
    /// iterate. On entry `self.gy` holds g(x_prev) (gval), `self.yprev`
    /// holds the previous iterate (xold); the result is written to `x`.
    /// `x` doubles as the vtemp workspace exactly as in the C code.
    pub fn anderson_accelerate(&mut self, x: &mut NVector, iter: i32) {
        let maa = self.m as usize;
        let beta = self.beta;

        /* reset ipt_map, i_pt */
        for e in self.imap.iter_mut() {
            *e = 0;
        }
        let it = iter as usize;
        let i_pt = if iter > 0 { (it - 1) % maa } else { 0 };

        /* update dg[i_pt], df[i_pt], fv, gold and fold */
        // fv = gval - xold  (fv lives in self.delta)
        N_VLinearSum(1.0, &self.gy, -1.0, &self.yprev, &mut self.delta);
        if iter > 0 {
            // dg_new = gval - gold ; df_new = fv - fold
            N_VLinearSum(1.0, &self.gy, -1.0, &self.gold, &mut self.dg[i_pt]);
            N_VLinearSum(1.0, &self.delta, -1.0, &self.fold, &mut self.df[i_pt]);
        }
        self.gold.data.copy_from_slice(&self.gy.data);
        self.fold.data.copy_from_slice(&self.delta.data);

        /* on first iteration, just do basic fixed-point update */
        if iter == 0 {
            x.data.copy_from_slice(&self.gy.data);
            return;
        }

        /* update data structures based on current iteration index */
        if iter == 1 {
            /* second iteration */
            let mut r0 = N_VDotProd(&self.df[i_pt], &self.df[i_pt]);
            r0 = SUNRsqrt(r0);
            self.R[0] = r0;
            N_VScale(1.0 / r0, &self.df[i_pt], &mut self.q[i_pt]);
            self.imap[0] = 0;
        } else if it <= maa {
            /* another iteration before we've reached maa */
            x.data.copy_from_slice(&self.df[i_pt].data); // vtemp = df[i_pt]
            for j in 0..(it - 1) {
                self.imap[j] = j;
                let r = N_VDotProd(&self.q[j], x);
                self.R[(it - 1) * maa + j] = r;
                x.linear_sum_with(1.0, -r, &self.q[j]); // vtemp -= R*Q[j]
            }
            let mut r = N_VDotProd(x, x);
            r = SUNRsqrt(r);
            self.R[(it - 1) * maa + (it - 1)] = r;
            if r == 0.0 {
                N_VScale(0.0, x, &mut self.q[i_pt]);
            } else {
                N_VScale(1.0 / r, x, &mut self.q[i_pt]);
            }
            self.imap[it - 1] = it - 1;
        } else {
            /* we've filled the acceleration subspace, so start recycling */
            /* delete left-most column vector from QR factorization */
            for i in 0..(maa - 1) {
                let a = self.R[(i + 1) * maa + i];
                let b = self.R[(i + 1) * maa + i + 1];
                let rtemp = SUNRsqrt(a * a + b * b);
                let c = a / rtemp;
                let s = b / rtemp;
                self.R[(i + 1) * maa + i] = rtemp;
                self.R[(i + 1) * maa + i + 1] = 0.0;
                if i < maa - 1 {
                    for j in (i + 2)..maa {
                        let a = self.R[j * maa + i];
                        let b = self.R[j * maa + i + 1];
                        let rtemp = c * a + s * b;
                        self.R[j * maa + i + 1] = -s * a + c * b;
                        self.R[j * maa + i] = rtemp;
                    }
                }
                // vtemp = c*Q[i] + s*Q[i+1]; Q[i+1] = -s*Q[i] + c*Q[i+1]; Q[i] = vtemp
                {
                    let (left, right) = self.q.split_at_mut(i + 1);
                    let qi = &mut left[i];
                    let qip1 = &mut right[0];
                    N_VLinearSum(c, qi, s, qip1, x);
                    qip1.linear_sum_with(c, -s, qi);
                    qi.data.copy_from_slice(&x.data);
                }
            }

            /* shift R to the left by one */
            for i in 1..maa {
                for j in 0..(maa - 1) {
                    self.R[(i - 1) * maa + j] = self.R[i * maa + j];
                }
            }

            /* add the new df vector */
            x.data.copy_from_slice(&self.df[i_pt].data);
            for j in 0..(maa - 1) {
                let r = N_VDotProd(&self.q[j], x);
                self.R[(maa - 1) * maa + j] = r;
                x.linear_sum_with(1.0, -r, &self.q[j]);
            }
            let mut r = N_VDotProd(x, x);
            r = SUNRsqrt(r);
            self.R[(maa - 1) * maa + maa - 1] = r;
            N_VScale(1.0 / r, x, &mut self.q[maa - 1]);

            /* update the iteration map */
            let mut j = 0;
            for i in (i_pt + 1)..maa {
                self.imap[j] = i;
                j += 1;
            }
            for i in 0..(i_pt + 1) {
                self.imap[j] = i;
                j += 1;
            }
        }

        /* solve least squares problem and update solution */
        let laa = if maa < it { maa } else { it };
        // gamma[j] = dot(fv, Q[j])  (N_VDotProdMulti)
        for j in 0..laa {
            self.gamma[j] = N_VDotProd(&self.delta, &self.q[j]);
        }

        /* accumulate the linear combination in the same order as the C
           fused N_VLinearCombination call: x = gval - sum(gamma_i dg_i) [...] */
        x.data.copy_from_slice(&self.gy.data); // cvals[0]=1, Xvecs[0]=gval
        let mut damping_terms: Vec<(f64, usize)> = Vec::new();
        for i in (0..laa).rev() {
            for j in (i + 1)..laa {
                self.gamma[i] -= self.R[j * maa + i] * self.gamma[j];
            }
            if self.gamma[i] == 0.0 {
                self.gamma[i] = 0.0;
            } else {
                self.gamma[i] /= self.R[i * maa + i];
            }
            let c = -self.gamma[i];
            let dgv = &self.dg[self.imap[i]];
            for (xk, dk) in x.data.iter_mut().zip(&dgv.data) {
                *xk += c * *dk;
            }
            if self.damping {
                damping_terms.push(((1.0 - beta) * self.gamma[i], self.imap[i]));
            }
        }

        /* if enabled, apply damping */
        if self.damping {
            let onembeta = 1.0 - beta;
            for (xk, fk) in x.data.iter_mut().zip(&self.delta.data) {
                *xk += -onembeta * *fk;
            }
            for (c, idx) in damping_terms {
                let dfv = &self.df[idx];
                for (xk, dk) in x.data.iter_mut().zip(&dfv.data) {
                    *xk += c * *dk;
                }
            }
        }
    }
}
