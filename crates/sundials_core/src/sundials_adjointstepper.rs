/* -----------------------------------------------------------------
 * Translation of sundials-7.7.0/src/sundials/sundials_adjointstepper.c
 * (+ include/sundials/sundials_adjointstepper.h and
 *  src/sundials/sundials_adjointstepper_impl.h).
 *
 * SUNAdjointStepper pairs a forward SUNStepper (for checkpoint
 * recomputation) with an adjoint SUNStepper (for the backward
 * integration) plus the shared checkpoint scheme.
 *
 * Ownership adaptations: the C struct holds POINTERS to the two
 * steppers and the checkpoint scheme (the scheme is also referenced
 * by the integrator that created it). Here the fields are owning
 * `Option<...>` slots. Destroy runs SUNStepper_Destroy only for the
 * steppers whose own_* flag is set, exactly as C; a non-owned
 * stepper can be retrieved (`.take()`) before Destroy — otherwise it
 * drops with the struct (in C it would stay alive with the caller).
 * `sunctx` is carried only by the constructor signature.
 *
 * C's SUNAdjointStepper_ReInit discards the SUNStepper_ReInit return
 * values (no SUNCheckCall) — mirrored.
 * -----------------------------------------------------------------*/

use crate::nvector_serial::NVector;
use crate::sundials_adjointcheckpointscheme::{
    SUNAdjointCheckpointScheme, SUNAdjointCheckpointScheme_EnableDense,
};
use crate::sundials_context::SUNContext;
use crate::sundials_errors::{SUNErrCode, SUN_SUCCESS};
use crate::sundials_stepper::{
    SUNStepper, SUNStepper_Destroy, SUNStepper_Evolve, SUNStepper_GetNumSteps,
    SUNStepper_OneStep, SUNStepper_ReInit, SUNStepper_Reset, SUNStepper_ResetCheckpointIndex,
    SUNStepper_SetStopTime,
};
use crate::sundials_types::{
    suncountertype, sunbooleantype, sunrealtype, SUNOutputFormat, UserData,
    SUN_OUTPUTFORMAT_TABLE, SUNFALSE, SUNTRUE,
};

/// struct SUNAdjointStepper_ (impl header)
pub struct SUNAdjointStepper {
    pub nrecompute: suncountertype,
    pub final_step_idx: suncountertype,

    pub adj_sunstepper: Option<SUNStepper>,
    pub fwd_sunstepper: Option<SUNStepper>,
    pub own_adj_sunstepper: sunbooleantype,
    pub own_fwd_sunstepper: sunbooleantype,
    pub checkpoint_scheme: Option<SUNAdjointCheckpointScheme>,

    pub user_data: UserData,
    pub content: UserData,

    pub tf: sunrealtype,
}

#[allow(clippy::too_many_arguments)]
pub fn SUNAdjointStepper_Create(
    fwd_sunstepper: SUNStepper,
    own_fwd: sunbooleantype,
    adj_sunstepper: SUNStepper,
    own_adj: sunbooleantype,
    final_step_idx: suncountertype,
    tf: sunrealtype,
    _sf: &NVector,
    checkpoint_scheme: Option<SUNAdjointCheckpointScheme>,
    _sunctx: &SUNContext,
    adj_stepper_ptr: &mut Option<SUNAdjointStepper>,
) -> SUNErrCode {
    *adj_stepper_ptr = Some(SUNAdjointStepper {
        fwd_sunstepper: Some(fwd_sunstepper),
        own_fwd_sunstepper: own_fwd,
        adj_sunstepper: Some(adj_sunstepper),
        own_adj_sunstepper: own_adj,
        checkpoint_scheme,

        tf,
        final_step_idx,

        nrecompute: 0,

        user_data: None,
        content: None,
    });

    SUN_SUCCESS
}

pub fn SUNAdjointStepper_ReInit(
    self_: &mut SUNAdjointStepper,
    t0: sunrealtype,
    y0: &NVector,
    tf: sunrealtype,
    sf: &NVector,
) -> SUNErrCode {
    self_.tf = tf;
    self_.nrecompute = 0;
    /* C discards both return values here */
    let _ = SUNStepper_ReInit(self_.adj_sunstepper.as_mut().unwrap(), tf, sf);
    let _ = SUNStepper_ReInit(self_.fwd_sunstepper.as_mut().unwrap(), t0, y0);
    SUN_SUCCESS
}

pub fn SUNAdjointStepper_Evolve(
    self_: &mut SUNAdjointStepper,
    tout: sunrealtype,
    sens: &mut NVector,
    tret: &mut sunrealtype,
) -> SUNErrCode {
    let err = SUNStepper_Evolve(self_.adj_sunstepper.as_mut().unwrap(), tout, sens, tret);
    if err != SUN_SUCCESS {
        return err;
    }
    SUN_SUCCESS
}

pub fn SUNAdjointStepper_OneStep(
    self_: &mut SUNAdjointStepper,
    tout: sunrealtype,
    sens: &mut NVector,
    tret: &mut sunrealtype,
) -> SUNErrCode {
    let err = SUNStepper_OneStep(self_.adj_sunstepper.as_mut().unwrap(), tout, sens, tret);
    if err != SUN_SUCCESS {
        return err;
    }
    SUN_SUCCESS
}

pub fn SUNAdjointStepper_RecomputeFwd(
    self_: &mut SUNAdjointStepper,
    start_idx: suncountertype,
    t0: sunrealtype,
    y0: &mut NVector,
    tf: sunrealtype,
) -> SUNErrCode {
    let retcode = SUN_SUCCESS;

    let mut fwd_t = t0;
    let fwd_stepper = self_.fwd_sunstepper.as_mut().unwrap();
    let err = SUNStepper_Reset(fwd_stepper, t0, y0);
    if err != SUN_SUCCESS {
        return err;
    }
    let err = SUNStepper_ResetCheckpointIndex(fwd_stepper, start_idx);
    if err != SUN_SUCCESS {
        return err;
    }

    let err =
        SUNAdjointCheckpointScheme_EnableDense(self_.checkpoint_scheme.as_mut().unwrap(), true);
    if err != SUN_SUCCESS {
        return err;
    }

    let fwd_stepper = self_.fwd_sunstepper.as_mut().unwrap();
    let err = SUNStepper_SetStopTime(fwd_stepper, tf);
    if err != SUN_SUCCESS {
        return err;
    }

    let mut nst_before: suncountertype = 0;
    let mut nst_after: suncountertype = 0;
    let err = SUNStepper_GetNumSteps(fwd_stepper, &mut nst_before);
    if err != SUN_SUCCESS {
        return err;
    }
    let err = SUNStepper_Evolve(fwd_stepper, tf, y0, &mut fwd_t);
    if err != SUN_SUCCESS {
        return err;
    }
    let err = SUNStepper_GetNumSteps(fwd_stepper, &mut nst_after);
    if err != SUN_SUCCESS {
        return err;
    }
    self_.nrecompute += nst_after - nst_before;

    let err =
        SUNAdjointCheckpointScheme_EnableDense(self_.checkpoint_scheme.as_mut().unwrap(), false);
    if err != SUN_SUCCESS {
        return err;
    }

    retcode
}

pub fn SUNAdjointStepper_Destroy(self_ptr: &mut Option<SUNAdjointStepper>) -> SUNErrCode {
    if let Some(self_) = self_ptr.as_mut() {
        if self_.own_fwd_sunstepper {
            SUNStepper_Destroy(&mut self_.fwd_sunstepper);
        }
        if self_.own_adj_sunstepper {
            SUNStepper_Destroy(&mut self_.adj_sunstepper);
        }
        *self_ptr = None;
    }
    SUN_SUCCESS
}

pub fn SUNAdjointStepper_SetUserData(
    self_: &mut SUNAdjointStepper,
    user_data: UserData,
) -> SUNErrCode {
    self_.user_data = user_data;
    SUN_SUCCESS
}

pub fn SUNAdjointStepper_GetNumSteps(
    self_: &mut SUNAdjointStepper,
    num_steps: &mut suncountertype,
) -> SUNErrCode {
    SUNStepper_GetNumSteps(self_.adj_sunstepper.as_mut().unwrap(), num_steps)
}

pub fn SUNAdjointStepper_GetNumRecompute(
    self_: &SUNAdjointStepper,
    num_recompute: &mut suncountertype,
) -> SUNErrCode {
    *num_recompute = self_.nrecompute;
    SUN_SUCCESS
}

/* sunfprintf_long (src/sundials/sundials_utils.h), private per the
   established per-module pattern */
const SUN_TABLE_WIDTH: usize = 29;

fn sunfprintf_long(
    outfile: &mut dyn std::io::Write,
    fmt: SUNOutputFormat,
    start: bool,
    name: &str,
    value: i64,
) {
    if fmt == SUN_OUTPUTFORMAT_TABLE {
        let _ = writeln!(outfile, "{:<width$} = {}", name, value, width = SUN_TABLE_WIDTH);
    } else {
        if !start {
            let _ = write!(outfile, ",");
        }
        let _ = write!(outfile, "{},{}", name, value);
    }
}

pub fn SUNAdjointStepper_PrintAllStats(
    self_: &mut SUNAdjointStepper,
    outfile: &mut dyn std::io::Write,
    fmt: SUNOutputFormat,
) -> SUNErrCode {
    let mut nst: suncountertype = 0;
    let err = SUNStepper_GetNumSteps(self_.adj_sunstepper.as_mut().unwrap(), &mut nst);
    if err != SUN_SUCCESS {
        return err;
    }
    sunfprintf_long(outfile, fmt, SUNTRUE, "Num backwards steps", nst);
    sunfprintf_long(outfile, fmt, SUNFALSE, "Num recompute steps", self_.nrecompute);

    SUN_SUCCESS
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::nvector_serial::{N_VConst, N_VNew_Serial};
    use crate::sunadjointcheckpointscheme_fixed::SUNAdjointCheckpointScheme_Create_Fixed;
    use crate::sundials_context::SUNContext_Create;
    use crate::sundials_stepper::*;
    use crate::sundials_system_memory::SUNMemoryHelper_Sys;
    use crate::sundials_types::{SUNDATAIOMODE_INMEM, SUN_OUTPUTFORMAT_CSV};

    /* dummy forward stepper: content holds the step counter */
    fn fwd_evolve(
        stepper: &mut SUNStepper,
        tout: sunrealtype,
        _vret: &mut NVector,
        tret: &mut sunrealtype,
    ) -> SUNErrCode {
        *tret = tout;
        let nst = stepper
            .content
            .as_mut()
            .unwrap()
            .downcast_mut::<suncountertype>()
            .unwrap();
        *nst += 3; /* pretend the recompute took 3 steps */
        SUN_SUCCESS
    }

    fn fwd_getnumsteps(stepper: &mut SUNStepper, nst: &mut suncountertype) -> SUNErrCode {
        *nst = *stepper
            .content
            .as_ref()
            .unwrap()
            .downcast_ref::<suncountertype>()
            .unwrap();
        SUN_SUCCESS
    }

    fn ok_reset(_s: &mut SUNStepper, _t: sunrealtype, _v: &NVector) -> SUNErrCode {
        SUN_SUCCESS
    }
    fn ok_rci(_s: &mut SUNStepper, _i: suncountertype) -> SUNErrCode {
        SUN_SUCCESS
    }
    fn ok_stop(_s: &mut SUNStepper, _t: sunrealtype) -> SUNErrCode {
        SUN_SUCCESS
    }

    /* dummy adjoint stepper */
    fn adj_evolve(
        _stepper: &mut SUNStepper,
        tout: sunrealtype,
        vret: &mut NVector,
        tret: &mut sunrealtype,
    ) -> SUNErrCode {
        N_VConst(-1.0, vret);
        *tret = tout;
        SUN_SUCCESS
    }

    fn adj_getnumsteps(_stepper: &mut SUNStepper, nst: &mut suncountertype) -> SUNErrCode {
        *nst = 11;
        SUN_SUCCESS
    }

    fn make_adjoint_stepper() -> SUNAdjointStepper {
        let ctx = SUNContext_Create();
        let helper = SUNMemoryHelper_Sys(&ctx);

        let mut fwd = SUNStepper_Create(&ctx);
        SUNStepper_SetContent(&mut fwd, Some(Box::new(0 as suncountertype)));
        SUNStepper_SetEvolveFn(&mut fwd, fwd_evolve);
        SUNStepper_SetGetNumStepsFn(&mut fwd, fwd_getnumsteps);
        SUNStepper_SetResetFn(&mut fwd, ok_reset);
        SUNStepper_SetResetCheckpointIndexFn(&mut fwd, ok_rci);
        SUNStepper_SetStopTimeFn(&mut fwd, ok_stop);

        let mut adj = SUNStepper_Create(&ctx);
        SUNStepper_SetEvolveFn(&mut adj, adj_evolve);
        SUNStepper_SetOneStepFn(&mut adj, adj_evolve);
        SUNStepper_SetGetNumStepsFn(&mut adj, adj_getnumsteps);

        let mut scheme = None;
        assert_eq!(
            SUNAdjointCheckpointScheme_Create_Fixed(
                SUNDATAIOMODE_INMEM,
                &helper,
                1,
                8,
                true,
                &ctx,
                &mut scheme
            ),
            SUN_SUCCESS
        );

        let sf = N_VNew_Serial(2, &ctx);
        let mut adj_stepper = None;
        assert_eq!(
            SUNAdjointStepper_Create(fwd, SUNTRUE, adj, SUNTRUE, 100, 10.0, &sf, scheme,
                                     &ctx, &mut adj_stepper),
            SUN_SUCCESS
        );
        adj_stepper.unwrap()
    }

    #[test]
    fn evolve_delegates_and_stats_print() {
        let ctx = SUNContext_Create();
        let mut s = make_adjoint_stepper();
        assert_eq!(s.tf, 10.0);
        assert_eq!(s.final_step_idx, 100);

        let mut sens = N_VNew_Serial(2, &ctx);
        let mut tret = 0.0;
        assert_eq!(SUNAdjointStepper_Evolve(&mut s, 4.0, &mut sens, &mut tret), SUN_SUCCESS);
        assert_eq!(tret, 4.0);
        assert_eq!(sens.data, vec![-1.0; 2]);

        let mut nst = 0;
        assert_eq!(SUNAdjointStepper_GetNumSteps(&mut s, &mut nst), SUN_SUCCESS);
        assert_eq!(nst, 11);

        let mut out = Vec::new();
        assert_eq!(
            SUNAdjointStepper_PrintAllStats(&mut s, &mut out, SUN_OUTPUTFORMAT_CSV),
            SUN_SUCCESS
        );
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "Num backwards steps,11,Num recompute steps,0"
        );

        let mut slot = Some(s);
        assert_eq!(SUNAdjointStepper_Destroy(&mut slot), SUN_SUCCESS);
        assert!(slot.is_none());
    }

    #[test]
    fn recompute_fwd_counts_steps() {
        let ctx = SUNContext_Create();
        let mut s = make_adjoint_stepper();

        let mut y0 = N_VNew_Serial(2, &ctx);
        assert_eq!(SUNAdjointStepper_RecomputeFwd(&mut s, 0, 0.0, &mut y0, 5.0), SUN_SUCCESS);
        let mut nrec = 0;
        assert_eq!(SUNAdjointStepper_GetNumRecompute(&s, &mut nrec), SUN_SUCCESS);
        assert_eq!(nrec, 3);

        /* second recompute accumulates */
        assert_eq!(SUNAdjointStepper_RecomputeFwd(&mut s, 3, 1.0, &mut y0, 8.0), SUN_SUCCESS);
        assert_eq!(SUNAdjointStepper_GetNumRecompute(&s, &mut nrec), SUN_SUCCESS);
        assert_eq!(nrec, 6);

        /* ReInit resets the recompute counter (return values discarded as in C) */
        let sf = N_VNew_Serial(2, &ctx);
        assert_eq!(SUNAdjointStepper_ReInit(&mut s, 0.0, &y0, 9.0, &sf), SUN_SUCCESS);
        assert_eq!(s.tf, 9.0);
        assert_eq!(SUNAdjointStepper_GetNumRecompute(&s, &mut nrec), SUN_SUCCESS);
        assert_eq!(nrec, 0);
    }
}
