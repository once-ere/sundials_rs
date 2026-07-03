/* -----------------------------------------------------------------------------
 * Translated from src/cvodes/cvodes_resize.c (CVODES 7.7.0).
 * Build Nordsieck array from solution history.
 *
 * In C, resizing replaces each internal vector by cloning the (already
 * resized) user history vector y_hist[0]; here the internal vectors are
 * replaced by fresh NVectors of the new length. The C NULL-pointer input
 * checks become slice-length checks. (The C source is identical to the
 * CVODE cvode_resize.c apart from the includes; only the state, not the
 * quadrature/sensitivity extensions, is resized — matching the C.)
 * ---------------------------------------------------------------------------*/
use crate::cvodes_impl::*;
use crate::cvodes_nls::CVodeSetNonlinearSolver;
use crate::nvector_serial::*;
use crate::sundials_types::*;
use crate::sunnonlinsol_newton::SUNNonlinSol_Newton;

const ZERO: f64 = 0.0; /* real 0.0 */
const ONE: f64 = 1.0; /* real 1.0 */

/* -----------------------------------------------------------------------------
 * Build Adams Nordsieck array from f(t,y) history and y(t) value
 * ---------------------------------------------------------------------------*/

fn cvBuildNordsieckArrayAdams(
    t: &[f64],
    y: &NVector,
    f: &[NVector],
    wrk: &mut [NVector],
    order: i32,
    hscale: f64,
    zn: &mut [NVector],
) -> i32 {
    /* Check for valid inputs (C checks the pointers for NULL; here the
       corresponding check is that the slices are long enough) */
    if order < 1 {
        return CV_ILL_INPUT;
    }
    let order = order as usize;

    if t.len() < order || f.len() < order || wrk.len() < order {
        return CV_ILL_INPUT;
    }

    /* Compute Nordsieck array */
    if order > 1 {
        /* Compute Newton polynomial coefficients interpolating f history */
        for i in 0..order {
            let (wrk_i, f_i) = (&mut wrk[i], &f[i]);
            N_VScale(ONE, f_i, wrk_i);
        }

        for i in 1..order {
            for j in (i..order).rev() {
                /* Divided difference */
                let delta_t = ONE / (t[j - i] - t[j]);
                /* N_VLinearSum(delta_t, wrk[j-1], -delta_t, wrk[j], wrk[j]) */
                let (left, right) = wrk.split_at_mut(j);
                right[0].linear_sum_with(-delta_t, delta_t, &left[j - 1]);
            }
        }

        /* Compute derivatives of Newton polynomial of f history */
        N_VScale(ONE, &wrk[order - 1], &mut zn[1]);
        for i in 2..=order {
            N_VConst(ZERO, &mut zn[i]);
        }

        for i in (0..=order - 2).rev() {
            for j in (1..order).rev() {
                /* N_VLinearSum(t[0] - t[i], zn[j+1], j, zn[j], zn[j+1]) */
                let (left, right) = zn.split_at_mut(j + 1);
                right[0].linear_sum_with(t[0] - t[i], j as f64, &left[j]);
            }
            /* N_VLinearSum(t[0] - t[i], zn[1], ONE, wrk[i], zn[1]) */
            zn[1].linear_sum_with(t[0] - t[i], ONE, &wrk[i]);
        }
    }

    /* Overwrite first two columns with input values */
    N_VScale(ONE, y, &mut zn[0]);
    N_VScale(ONE, &f[0], &mut zn[1]);

    /* Scale entries */
    let mut scale = ONE;
    for i in 1..=order {
        scale *= hscale / (i as f64);
        zn[i].scale_inplace(scale);
    }

    CV_SUCCESS
}

/* -----------------------------------------------------------------------------
 * Build BDF Nordsieck array from y(t) history and f(t,y) value
 * ---------------------------------------------------------------------------*/

fn cvBuildNordsieckArrayBDF(
    t: &[f64],
    y: &[NVector],
    f: &NVector,
    wrk: &mut [NVector],
    order: i32,
    hscale: f64,
    zn: &mut [NVector],
) -> i32 {
    /* Check for valid inputs (C checks the pointers for NULL; here the
       corresponding check is that the slices are long enough) */
    if order < 1 {
        return CV_ILL_INPUT;
    }
    let order = order as usize;

    if t.len() < order || y.len() < order || wrk.len() < order + 1 {
        return CV_ILL_INPUT;
    }

    /* Compute Nordsieck array */
    if order > 1 {
        /* Setup extended array of times to incorporate derivative value */
        let mut t_ext = [ZERO; BDF_Q_MAX + 1];

        t_ext[0] = t[0];
        for i in 1..=order {
            t_ext[i] = t[i - 1];
        }

        /* Compute Hermite polynomial coefficients interpolating y history and f */
        N_VScale(ONE, &y[0], &mut wrk[0]);
        for i in 1..=order {
            let (wrk_i, y_im1) = (&mut wrk[i], &y[i - 1]);
            N_VScale(ONE, y_im1, wrk_i);
        }

        for i in 1..=order {
            for j in ((i)..=order).rev() {
                /* j > i - 1  <=>  j >= i */
                if i == 1 && j == 1 {
                    /* Replace with actual derivative value */
                    N_VScale(ONE, f, &mut wrk[j]);
                } else {
                    /* Divided difference */
                    let delta_t = ONE / (t_ext[j - i] - t_ext[j]);
                    /* N_VLinearSum(delta_t, wrk[j-1], -delta_t, wrk[j], wrk[j]) */
                    let (left, right) = wrk.split_at_mut(j);
                    right[0].linear_sum_with(-delta_t, delta_t, &left[j - 1]);
                }
            }
        }

        /* Compute derivatives of Hermite polynomial */
        N_VScale(ONE, &wrk[order], &mut zn[0]);
        for i in 1..=order {
            N_VConst(ZERO, &mut zn[i]);
        }

        for i in (0..=order - 1).rev() {
            for j in (1..=order).rev() {
                /* N_VLinearSum(t_ext[0] - t_ext[i], zn[j], j, zn[j-1], zn[j]) */
                let (left, right) = zn.split_at_mut(j);
                right[0].linear_sum_with(t_ext[0] - t_ext[i], j as f64, &left[j - 1]);
            }
            /* N_VLinearSum(t_ext[0] - t_ext[i], zn[0], ONE, wrk[i], zn[0]) */
            zn[0].linear_sum_with(t_ext[0] - t_ext[i], ONE, &wrk[i]);
        }
    }

    /* Overwrite first two columns with input values */
    N_VScale(ONE, &y[0], &mut zn[0]);
    N_VScale(ONE, f, &mut zn[1]);

    /* Scale entries */
    let mut scale = ONE;
    for i in 1..=order {
        scale *= hscale / (i as f64);
        zn[i].scale_inplace(scale);
    }

    CV_SUCCESS
}

/* -----------------------------------------------------------------------------
 * Compute predicted new state (simplified cvPredict for k = 1 and j = q...1)
 * ---------------------------------------------------------------------------*/

fn cvPredictY(order: i32, zn: &[NVector], ypred: &mut NVector) -> i32 {
    N_VScale(ONE, &zn[0], ypred);
    for j in 1..=order as usize {
        /* N_VLinearSum(ONE, zn[j], ONE, ypred, ypred) */
        ypred.linear_sum_with(ONE, ONE, &zn[j]);
    }

    CV_SUCCESS
}

/* -----------------------------------------------------------------------------
 * Resize CVODE and build new history array
 * ---------------------------------------------------------------------------*/

pub fn CVodeResizeHistory(
    cv_mem: &mut CVodeMem,
    t_hist: &[f64],
    y_hist: &[NVector],
    f_hist: &[NVector],
    num_y_hist: i32,
    num_f_hist: i32,
) -> i32 {
    /* ------------ *
     * Check inputs *
     * ------------ */

    if t_hist.is_empty() {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeResizeHistory", file!(),
                       "Time history array is NULL");
        return CV_ILL_INPUT;
    }

    if y_hist.is_empty() {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeResizeHistory", file!(),
                       "State history array is NULL");
        return CV_ILL_INPUT;
    }

    if f_hist.is_empty() {
        cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeResizeHistory", file!(),
                       "RHS history array is NULL");
        return CV_ILL_INPUT;
    }

    /* Check that the input history is sufficient for the current (next) order */
    let n_hist = std::cmp::min(cv_mem.cv_q + 1, cv_mem.cv_qmax);

    if cv_mem.cv_lmm == CV_ADAMS {
        if num_y_hist < 2 {
            cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeResizeHistory", file!(),
                           "Insufficient solution history");
            return CV_ILL_INPUT;
        }

        for i in 0..n_hist {
            if i as usize >= f_hist.len() || i >= num_f_hist {
                cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeResizeHistory", file!(),
                               "Insufficient right-hand side history");
                return CV_ILL_INPUT;
            }
        }
    } else {
        if num_f_hist < 2 {
            cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeResizeHistory", file!(),
                           "Insufficient right-hand side history");
            return CV_ILL_INPUT;
        }

        for i in 0..n_hist {
            if i as usize >= y_hist.len() || i >= num_y_hist {
                cvProcessError(Some(cv_mem), CV_ILL_INPUT, line!(), "CVodeResizeHistory", file!(),
                               "Insufficient solution history");
                return CV_ILL_INPUT;
            }
        }
    }

    /* -------------- *
     * Resize vectors *
     * -------------- */

    let new_len = y_hist[0].len();

    cv_mem.cv_ewt = NVector::new(new_len);
    cv_mem.cv_acor = NVector::new(new_len);
    cv_mem.cv_tempv = NVector::new(new_len);
    cv_mem.cv_ftemp = NVector::new(new_len);
    cv_mem.cv_vtemp1 = NVector::new(new_len);
    cv_mem.cv_vtemp2 = NVector::new(new_len);
    cv_mem.cv_vtemp3 = NVector::new(new_len);

    /* cv_y holds the current solver state; in C it is not resized here
       because it aliases user-supplied storage, but in this port it is
       owned storage that must match the new problem size */
    cv_mem.cv_y = NVector::new(new_len);

    /* User will need to set a new vector of absolute tolerances */
    if cv_mem.cv_VabstolMallocDone {
        cv_mem.cv_Vabstol = NVector::new(new_len);
    }

    /* User will need to set a new constraints vector
       (CVODES tests constraint presence through the cv_constraints
       vector itself: empty = C NULL) */
    if !cv_mem.cv_constraints.data.is_empty() {
        cv_mem.cv_constraints = NVector::default();
    }

    for j in 0..=cv_mem.cv_qmax_alloc as usize {
        cv_mem.cv_zn[j] = NVector::new(new_len);
    }

    /* ----------------------- *
     * Resize nonlinear solver *
     * ----------------------- */

    if cv_mem.NLS.is_some() && cv_mem.ownNLS {
        /* Destroying the old solver is handled by RAII on overwrite */
        cv_mem.NLS = None;
        cv_mem.ownNLS = SUNFALSE;

        let NLS = SUNNonlinSol_Newton(&y_hist[0], &cv_mem.cv_sunctx);

        let retval = CVodeSetNonlinearSolver(cv_mem, NLS);
        if retval != 0 {
            cvProcessError(Some(cv_mem), CV_MEM_FAIL, line!(), "CVodeResizeHistory", file!(),
                           "Error attaching default Newton solver");
            return CV_MEM_FAIL;
        }
        cv_mem.ownNLS = SUNTRUE;
    }

    /* ----------------------------- *
     * Create workspace for resizing *
     * ----------------------------- */

    let mut wrk_space_size = std::cmp::max(cv_mem.cv_q, cv_mem.cv_qprime);
    if cv_mem.cv_lmm == CV_BDF {
        wrk_space_size += 1;
    }

    let mut resize_wrk: Vec<NVector> = Vec::with_capacity(wrk_space_size as usize);
    for _ in 0..wrk_space_size {
        resize_wrk.push(NVector::new(new_len));
    }

    /* ------------------------------------------------------------------------ *
     * Construct Nordsieck array at the old time but with the new size to
     * compute correction vector at the new state size.
     * ------------------------------------------------------------------------ */

    if cv_mem.cv_q < cv_mem.cv_qmax {
        /* Compute z_{n-1} with new history size */
        let retval = if cv_mem.cv_lmm == CV_ADAMS {
            cvBuildNordsieckArrayAdams(&t_hist[1..], &y_hist[1], &f_hist[1..],
                                       &mut resize_wrk, cv_mem.cv_q, cv_mem.cv_hscale,
                                       &mut cv_mem.cv_zn)
        } else {
            cvBuildNordsieckArrayBDF(&t_hist[1..], &y_hist[1..], &f_hist[1],
                                     &mut resize_wrk, cv_mem.cv_q, cv_mem.cv_hscale,
                                     &mut cv_mem.cv_zn)
        };

        if retval != 0 {
            cvProcessError(Some(cv_mem), retval, line!(), "CVodeResizeHistory", file!(),
                           "Building the Nordsieck array failed");
            return retval;
        }

        /* Get predicted value */
        let retval = {
            let CVodeMem { cv_zn, cv_vtemp1, .. } = cv_mem;
            cvPredictY(cv_mem.cv_q, cv_zn, cv_vtemp1)
        };

        if retval != 0 {
            cvProcessError(Some(cv_mem), retval, line!(), "CVodeResizeHistory", file!(),
                           "Computing the predictor failed");
            return retval;
        }

        /* Resized correction */
        {
            let qmax = cv_mem.cv_qmax as usize;
            let CVodeMem { cv_zn, cv_vtemp1, .. } = cv_mem;
            N_VLinearSum(ONE, &y_hist[0], -ONE, cv_vtemp1, &mut cv_zn[qmax]);
        }
    }

    /* ----------------------------- *
     * Construct new Nordsieck Array *
     * ----------------------------- */

    let retval = if cv_mem.cv_lmm == CV_ADAMS {
        cvBuildNordsieckArrayAdams(t_hist, &y_hist[0], f_hist, &mut resize_wrk,
                                   cv_mem.cv_qprime, cv_mem.cv_hscale, &mut cv_mem.cv_zn)
    } else {
        cvBuildNordsieckArrayBDF(t_hist, y_hist, &f_hist[0], &mut resize_wrk,
                                 cv_mem.cv_qprime, cv_mem.cv_hscale, &mut cv_mem.cv_zn)
    };

    if retval != 0 {
        cvProcessError(Some(cv_mem), retval, line!(), "CVodeResizeHistory", file!(),
                       "Building the Nordsieck array failed");
        return retval;
    }

    /* ------------------- *
     * Update time history *
     * ------------------- */

    /* Ensure internal time and step history match the input history */
    cv_mem.cv_tn = t_hist[0];

    for i in 1..n_hist as usize {
        /* (C indexes t_hist unconditionally; guard keeps this in safe bounds) */
        if i < t_hist.len() {
            cv_mem.cv_tau[i] = t_hist[i - 1] - t_hist[i];
        }
    }

    /* In the next step, perform initialization needed after a resize */
    cv_mem.first_step_after_resize = SUNTRUE;

    /* Workspace for resizing is destroyed when resize_wrk drops */

    CV_SUCCESS
}
