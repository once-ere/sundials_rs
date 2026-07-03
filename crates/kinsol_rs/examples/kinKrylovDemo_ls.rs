/* -----------------------------------------------------------------
 * Translated from examples/kinsol/serial/kinKrylovDemo_ls.c (KINSOL 7.7.0)
 *
 * This example loops through the available iterative linear solvers:
 * SPGMR, SPBCGS, SPTFQMR, and SPFGMR.
 *
 * Example (serial):
 *
 * This example solves a nonlinear system that arises from a system
 * of partial differential equations. The PDE system is a food web
 * population model, with predator-prey interaction and diffusion
 * on the unit square in two dimensions. The dependent variable
 * vector is the following:
 *
 *       1   2         ns
 * c = (c , c ,  ..., c  )     (denoted by the variable cc)
 *
 * and the PDE's are as follows:
 *
 *                    i       i
 *         0 = d(i)*(c     + c    )  +  f  (x,y,c)   (i=1,...,ns)
 *                    xx      yy         i
 *
 *   where
 *
 *                   i             ns         j
 *   f  (x,y,c)  =  c  * (b(i)  + sum a(i,j)*c )
 *    i                           j=1
 *
 * The number of species is ns = 2 * np, with the first np being
 * prey and the last np being predators. The number np is both the
 * number of prey and predator species. The coefficients a(i,j),
 * b(i), d(i) are:
 *
 *   a(i,i) = -AA   (all i)
 *   a(i,j) = -GG   (i <= np , j >  np)
 *   a(i,j) =  EE   (i >  np,  j <= np)
 *   b(i) = BB * (1 + alpha * x * y)   (i <= np)
 *   b(i) =-BB * (1 + alpha * x * y)   (i >  np)
 *   d(i) = DPREY   (i <= np)
 *   d(i) = DPRED   ( i > np)
 *
 * The various scalar parameters are set using define's or in
 * routine InitUserData.
 *
 * The boundary conditions are: normal derivative = 0, and the
 * initial guess is constant in x and y, but the final solution
 * is not.
 *
 * The PDEs are discretized by central differencing on an MX by
 * MY mesh.
 *
 * The nonlinear system is solved by KINSOL using the method
 * specified in local variable globalstrat.
 *
 * The preconditioner matrix is a block-diagonal matrix based on
 * the partial derivatives of the interaction terms f only.
 *
 * Constraints are imposed to make all components of the solution
 * positive.
 * -----------------------------------------------------------------
 * References:
 *
 * 1. Peter N. Brown and Youcef Saad,
 *    Hybrid Krylov Methods for Nonlinear Systems of Equations
 *    LLNL report UCRL-97645, November 1987.
 *
 * 2. Peter N. Brown and Alan C. Hindmarsh,
 *    Reduced Storage Matrix Methods in Stiff ODE systems,
 *    Lawrence Livermore National Laboratory Report  UCRL-95088,
 *    Rev. 1, June 1987, and  Journal of Applied Mathematics and
 *    Computation, Vol. 31 (May 1989), pp. 40-91. (Presents a
 *    description of the time-dependent version of this test
 *    problem.)
 * -----------------------------------------------------------------*/
#![allow(non_snake_case)]
#![allow(non_upper_case_globals)]

use kinsol_rs::sundials_dense::{SUNDlsMat_denseGETRF, SUNDlsMat_denseGETRS};
use kinsol_rs::sundials_utils::fmt_g;
use kinsol_rs::*;

/* Problem Constants */

/* must equal 2*(number of prey or predators)
number of prey = number of predators       */
const NUM_SPECIES: i64 = 6;

const MX: i64 = 5; /* MX = number of x mesh points */
const MY: i64 = 5; /* MY = number of y mesh points */
const NSMX: i64 = NUM_SPECIES * MX;
const NEQ: i64 = NSMX * MY; /* number of equations in the system */
const AA: f64 = 1.0; /* value of coefficient AA in above eqns */
const EE: f64 = 10000.; /* value of coefficient EE in above eqns */
const GG: f64 = 0.5e-6; /* value of coefficient GG in above eqns */
const BB: f64 = 1.0; /* value of coefficient BB in above eqns */
const DPREY: f64 = 1.0; /* value of coefficient dprey above */
const DPRED: f64 = 0.5; /* value of coefficient dpred above */
const ALPHA: f64 = 1.0; /* value of coefficient alpha above */
const AX: f64 = 1.0; /* total range of x variable */
const AY: f64 = 1.0; /* total range of y variable */
const FTOL: f64 = 1.0e-7; /* ftol tolerance */
const STOL: f64 = 1.0e-13; /* stol tolerance */
const THOUSAND: f64 = 1000.0; /* one thousand */
const ZERO: f64 = 0.; /* 0. */
const ONE: f64 = 1.0; /* 1. */
const TWO: f64 = 2.0; /* 2. */
const PREYIN: f64 = 1.0; /* initial guess for prey concentrations. */
const PREDIN: f64 = 30000.0; /* initial guess for predator concs.      */

/* Linear Solver Loop Constants */

const USE_SPGMR: i32 = 0;
const USE_SPBCGS: i32 = 1;
const USE_SPTFQMR: i32 = 2;
const USE_SPFGMR: i32 = 3;

const NS: usize = NUM_SPECIES as usize;

/* User-defined vector access macro: IJ_V */

/* IJ_V is defined in order to translate from the underlying 3D structure
   of the dependent variable vector to the 1D storage scheme for an N-vector.
   IJ_V(i,j) returns the base index in a vector data array corresponding to
   indices is = 0, jx = i, jy = j (C's IJ_Vptr). */

macro_rules! IJ_V {
    ($i:expr, $j:expr) => {
        (($i) * NUM_SPECIES + ($j) * NSMX) as usize
    };
}

/* Type : UserDataStruct
   contains preconditioner blocks, pivot arrays, and problem constants */

struct UserDataStruct {
    P: Vec<Vec<DenseMatrix>>,   /* [MX][MY] NUM_SPECIES x NUM_SPECIES blocks */
    pivot: Vec<Vec<[i64; NS]>>, /* [MX][MY] pivot arrays */
    acoef: [[f64; NS]; NS],     /* acoef[i] = C's i-th column pointer acoef[i] */
    bcoef: [f64; NS],
    rates: NVector,
    cox: [f64; NS],
    coy: [f64; NS],
    ax: f64,
    ay: f64,
    dx: f64,
    dy: f64,
    uround: f64,
    sqruround: f64,
    #[allow(dead_code)]
    mx: i64,
    #[allow(dead_code)]
    my: i64,
    #[allow(dead_code)]
    ns: i64,
    np: i64,
}

/*
 *--------------------------------------------------------------------
 * MAIN PROGRAM
 *--------------------------------------------------------------------
 */

fn main() {
    /* Create the SUNDIALS context object for this simulation. */
    let sunctx = SUNContext_Create();

    /* Allocate memory, and set problem data, initial values, tolerances */
    let globalstrategy = KIN_NONE;

    /* (data->rates is allocated here rather than after cc/sc as in the
       C main; order is inconsequential) */
    let mut data = AllocUserData(&sunctx);
    InitUserData(&mut data);

    /* Create serial vectors of length NEQ */
    let mut cc = N_VNew_Serial(NEQ, &sunctx);
    let mut sc = N_VNew_Serial(NEQ, &sunctx);

    let mut constraints = N_VNew_Serial(NEQ, &sunctx);
    N_VConst(TWO, &mut constraints);

    let fnormtol = FTOL;
    let scsteptol = STOL;

    /* The single C data block is shared by all loop passes: hand it to
       KINSOL each pass and take it back after the solve */
    let mut data_holder: UserData = Some(Box::new(data));

    /* START: Loop through SPGMR, SPBCGS, SPTFQMR and SPFGMR linear solver modules */
    for linsolver in 0..4 {
        /* (Re-)Initialize user data */
        SetInitialProfiles(&mut cc, &mut sc);

        /* Call KINCreate/KINInit to initialize KINSOL:
        A pointer to KINSOL problem memory is returned and stored in kmem. */
        let mut kmem = KINCreate(&sunctx);

        /* Vector cc passed as template vector. */
        let mut flag = KINInit(&mut kmem, func, &cc);
        if check_flag(flag, "KINInit") {
            std::process::exit(1);
        }

        flag = KINSetUserData(&mut kmem, data_holder.take());
        if check_flag(flag, "KINSetUserData") {
            std::process::exit(1);
        }
        flag = KINSetConstraints(&mut kmem, Some(&constraints));
        if check_flag(flag, "KINSetConstraints") {
            std::process::exit(1);
        }
        flag = KINSetFuncNormTol(&mut kmem, fnormtol);
        if check_flag(flag, "KINSetFuncNormTol") {
            std::process::exit(1);
        }
        flag = KINSetScaledStepTol(&mut kmem, scsteptol);
        if check_flag(flag, "KINSetScaledStepTol") {
            std::process::exit(1);
        }

        /* Attach a linear solver module */
        let mut maxl = 0;
        let mut maxlrst = 0;
        match linsolver {
            /* (a) SPGMR */
            USE_SPGMR => {
                /* Print header */
                print!(" -------");
                print!(" \n| SPGMR |\n");
                print!(" -------\n");

                /* Create SUNLinSol_SPGMR object with right preconditioning and the
                maximum Krylov dimension maxl */
                maxl = 15;
                let mut LS = SUNLinSol_SPGMR(&cc, SUN_PREC_RIGHT, maxl, &sunctx);

                /* Set the maximum number of restarts (C calls
                   SUNLinSol_SPGMRSetMaxRestarts after KINSetLinearSolver;
                   the LS object moves into kmem there, so set it just
                   before attaching — the value is only read at solve time) */
                maxlrst = 2;
                if let LinearSolver::Spgmr(s) = &mut LS {
                    flag = s.set_max_restarts(maxlrst);
                    if check_flag(flag, "SUNLinSol_SPGMRSetMaxRestarts") {
                        std::process::exit(1);
                    }
                }

                /* Attach the linear solver to KINSOL */
                flag = KINSetLinearSolver(&mut kmem, LS, None);
                if check_flag(flag, "KINSetLinearSolver") {
                    std::process::exit(1);
                }
            }

            /* (b) SPBCGS */
            USE_SPBCGS => {
                /* Print header */
                print!(" --------");
                print!(" \n| SPBCGS |\n");
                print!(" --------\n");

                /* Create SUNLinSol_SPBCGS object with right preconditioning and the
                maximum Krylov dimension maxl */
                maxl = 15;
                let LS = SUNLinSol_SPBCGS(&cc, SUN_PREC_RIGHT, maxl, &sunctx);

                /* Attach the linear solver to KINSOL */
                flag = KINSetLinearSolver(&mut kmem, LS, None);
                if check_flag(flag, "KINSetLinearSolver") {
                    std::process::exit(1);
                }
            }

            /* (c) SPTFQMR */
            USE_SPTFQMR => {
                /* Print header */
                print!(" ---------");
                print!(" \n| SPTFQMR |\n");
                print!(" ---------\n");

                /* Create SUNLinSol_SPTFQMR object with right preconditioning and the
                maximum Krylov dimension maxl */
                maxl = 25;
                let LS = SUNLinSol_SPTFQMR(&cc, SUN_PREC_RIGHT, maxl, &sunctx);

                /* Attach the linear solver to KINSOL */
                flag = KINSetLinearSolver(&mut kmem, LS, None);
                if check_flag(flag, "KINSetLinearSolver") {
                    std::process::exit(1);
                }
            }

            /* (d) SPFGMR */
            USE_SPFGMR => {
                /* Print header */
                print!(" -------");
                print!(" \n| SPFGMR |\n");
                print!(" -------\n");

                /* Create SUNLinSol_SPFGMR object with right preconditioning and the
                maximum Krylov dimension maxl */
                maxl = 15;
                let mut LS = SUNLinSol_SPFGMR(&cc, SUN_PREC_RIGHT, maxl, &sunctx);

                /* Set the maximum number of restarts (C calls
                   SUNLinSol_SPGMRSetMaxRestarts on the SPFGMR object,
                   which writes the SPFGMR content's max_restarts field —
                   the first four content fields coincide) */
                maxlrst = 2;
                if let LinearSolver::Spfgmr(s) = &mut LS {
                    flag = s.set_max_restarts(maxlrst);
                    if check_flag(flag, "SUNLinSol_SPGMRSetMaxRestarts") {
                        std::process::exit(1);
                    }
                }

                /* Attach the linear solver to KINSOL */
                flag = KINSetLinearSolver(&mut kmem, LS, None);
                if check_flag(flag, "KINSetLinearSolver") {
                    std::process::exit(1);
                }
            }

            _ => {}
        }

        /* Set preconditioner functions */
        flag = KINSetPreconditioner(&mut kmem, Some(PrecSetupBD), Some(PrecSolveBD));
        if check_flag(flag, "KINSetPreconditioner") {
            std::process::exit(1);
        }

        /* Print out the problem size, solution parameters, initial guess. */
        PrintHeader(globalstrategy, maxl, maxlrst, fnormtol, scsteptol, linsolver);

        /* Call KINSol and print output concentration profile */
        flag = KINSol(
            &mut kmem,      /* KINSol memory block */
            &mut cc,        /* initial guess on input; solution vector */
            globalstrategy, /* global strategy choice */
            &sc,            /* scaling vector, for the variable cc */
            &sc,            /* scaling vector for function values fval */
        );
        if check_flag(flag, "KINSol") {
            std::process::exit(1);
        }

        print!("\n\nComputed equilibrium species concentrations:\n");
        PrintOutput(&cc);

        /* Print final statistics and free memory */
        PrintFinalStats(&mut kmem, linsolver);

        /* Recover the shared user data block before releasing kmem */
        data_holder = KINGetUserData(&mut kmem).take();

        KINFree(kmem);
    } /* END: Loop through SPGMR, SPBCGS, SPTFQMR, and SPFGMR linear solver modules */
}

/*
 *--------------------------------------------------------------------
 * FUNCTIONS CALLED BY KINSOL
 *--------------------------------------------------------------------
 */

/*
 * System function for predator-prey system
 */

fn func(cc: &NVector, fval: &mut NVector, user_data: &mut UserData) -> i32 {
    let data = user_data
        .as_mut()
        .unwrap()
        .downcast_mut::<UserDataStruct>()
        .unwrap();
    let delx = data.dx;
    let dely = data.dy;

    let ccdata = &cc.data;
    /* split borrows: rates written while acoef/bcoef/cox/coy read */
    let UserDataStruct {
        rates,
        acoef,
        bcoef,
        cox,
        coy,
        ..
    } = data;
    let rdata = &mut rates.data;
    let fdata = &mut fval.data;

    /* Loop over all mesh points, evaluating rate array at each point*/
    for jy in 0..MY {
        let yy = dely * jy as f64;

        /* Set lower/upper index shifts, special at boundaries. */
        let idyl = if jy != 0 { NSMX } else { -NSMX };
        let idyu = if jy != MY - 1 { NSMX } else { -NSMX };

        for jx in 0..MX {
            let xx = delx * jx as f64;

            /* Set left/right index shifts, special at boundaries. */
            let idxl = if jx != 0 { NUM_SPECIES } else { -NUM_SPECIES };
            let idxr = if jx != MX - 1 { NUM_SPECIES } else { -NUM_SPECIES };

            let loc = IJ_V!(jx, jy); /* cxy, rxy, fxy base index */

            /* Get species interaction rate array at (xx,yy) */
            WebRate(
                xx,
                yy,
                &ccdata[loc..loc + NS],
                &mut rdata[loc..loc + NS],
                acoef,
                bcoef,
            );

            for is in 0..NS {
                /* Differencing in x direction */
                let dcyli = ccdata[loc + is] - ccdata[(loc as i64 - idyl) as usize + is];
                let dcyui = ccdata[(loc as i64 + idyu) as usize + is] - ccdata[loc + is];

                /* Differencing in y direction */
                let dcxli = ccdata[loc + is] - ccdata[(loc as i64 - idxl) as usize + is];
                let dcxri = ccdata[(loc as i64 + idxr) as usize + is] - ccdata[loc + is];

                /* Compute the total rate value at (xx,yy) */
                fdata[loc + is] =
                    coy[is] * (dcyui - dcyli) + cox[is] * (dcxri - dcxli) + rdata[loc + is];
            } /* end of is loop */
        } /* end of jx loop */
    } /* end of jy loop */

    0
}

/*
 * Preconditioner setup routine. Generate and preprocess P.
 */

fn PrecSetupBD(
    cc: &NVector,
    cscale: &NVector,
    fval: &NVector,
    fscale: &NVector,
    user_data: &mut UserData,
) -> i32 {
    let data = user_data
        .as_mut()
        .unwrap()
        .downcast_mut::<UserDataStruct>()
        .unwrap();
    let delx = data.dx;
    let dely = data.dy;

    let uround = data.uround;
    let sqruround = data.sqruround;
    let mut fac = N_VWL2Norm(fval, fscale);
    let mut r0 = THOUSAND * uround * fac * NEQ as f64;
    if r0 == ZERO {
        r0 = ONE;
    }

    let ccdata = &cc.data;
    let csdata = &cscale.data;
    let UserDataStruct {
        P,
        pivot,
        rates,
        acoef,
        bcoef,
        ..
    } = data;
    let rdata = &rates.data;

    let mut perturb_rates = [0.0f64; NS];

    /* Loop over spatial points; get size NUM_SPECIES Jacobian block at each */
    for jy in 0..MY {
        let yy = jy as f64 * dely;

        for jx in 0..MX {
            let xx = jx as f64 * delx;
            let Pxy = &mut P[jx as usize][jy as usize];
            let loc = IJ_V!(jx, jy); /* cxy, scxy, ratesxy base index */

            /* C perturbs cc in place and restores it; cc is immutable
               here, so perturb a local copy of the species block */
            let mut cxy = [0.0f64; NS];
            cxy.copy_from_slice(&ccdata[loc..loc + NS]);

            /* Compute difference quotients of interaction rate fn. */
            for j in 0..NS {
                let csave = cxy[j]; /* Save the j,jx,jy element of cc */
                let a = sqruround * SUNRabs(csave);
                let b = r0 / csdata[loc + j];
                let r = if a > b { a } else { b }; /* MAX(A, B) */
                cxy[j] += r; /* Perturb the j,jx,jy element of cc */
                fac = ONE / r;

                WebRate(xx, yy, &cxy, &mut perturb_rates, acoef, bcoef);

                /* Restore j,jx,jy element of cc */
                cxy[j] = csave;

                /* Load the j-th column of difference quotients */
                let Pxycol = Pxy.col_mut(j as i64);
                for i in 0..NS {
                    Pxycol[i] = (perturb_rates[i] - rdata[loc + i]) * fac;
                }
            } /* end of j loop */

            /* Do LU decomposition of size NUM_SPECIES preconditioner block */
            let ret = SUNDlsMat_denseGETRF(Pxy, &mut pivot[jx as usize][jy as usize]);
            if ret != 0 {
                return 1;
            }
        } /* end of jx loop */
    } /* end of jy loop */

    0
}

/*
 * Preconditioner solve routine
 */

fn PrecSolveBD(
    _cc: &NVector,
    _cscale: &NVector,
    _fval: &NVector,
    _fscale: &NVector,
    vv: &mut NVector,
    user_data: &mut UserData,
) -> i32 {
    let data = user_data
        .as_mut()
        .unwrap()
        .downcast_mut::<UserDataStruct>()
        .unwrap();

    let vdata = &mut vv.data;

    for jx in 0..MX {
        for jy in 0..MY {
            /* For each (jx,jy), solve a linear system of size NUM_SPECIES.
            vxy is the address of the corresponding portion of the vector vv;
            Pxy is the address of the corresponding block of the matrix P;
            piv is the address of the corresponding block of the array pivot. */
            let loc = IJ_V!(jx, jy);
            let vxy = &mut vdata[loc..loc + NS];
            let Pxy = &data.P[jx as usize][jy as usize];
            let piv = &data.pivot[jx as usize][jy as usize];
            SUNDlsMat_denseGETRS(Pxy, piv, vxy);
        } /* end of jy loop */
    } /* end of jx loop */

    0
}

/*
 * Interaction rate function routine
 */

fn WebRate(
    xx: f64,
    yy: f64,
    cxy: &[f64],
    ratesxy: &mut [f64],
    acoef: &[[f64; NS]; NS],
    bcoef: &[f64; NS],
) {
    for i in 0..NS {
        ratesxy[i] = DotProd(NUM_SPECIES, cxy, &acoef[i]);
    }

    let fac = ONE + ALPHA * xx * yy;

    for i in 0..NS {
        ratesxy[i] = cxy[i] * (bcoef[i] * fac + ratesxy[i]);
    }
}

/*
 * Dot product routine for sunrealtype arrays
 */

fn DotProd(size: i64, x1: &[f64], x2: &[f64]) -> f64 {
    let mut temp = ZERO;
    for i in 0..size as usize {
        temp += x1[i] * x2[i];
    }
    temp
}

/*
 *--------------------------------------------------------------------
 * PRIVATE FUNCTIONS
 *--------------------------------------------------------------------
 */

/*
 * Allocate memory for data structure of type UserDataStruct
 */

fn AllocUserData(sunctx: &SUNContext) -> UserDataStruct {
    let mut P = Vec::with_capacity(MX as usize);
    let mut pivot = Vec::with_capacity(MX as usize);
    for _jx in 0..MX {
        let mut prow = Vec::with_capacity(MY as usize);
        let mut pvrow = Vec::with_capacity(MY as usize);
        for _jy in 0..MY {
            prow.push(DenseMatrix::new(NUM_SPECIES, NUM_SPECIES));
            pvrow.push([0i64; NS]);
        }
        P.push(prow);
        pivot.push(pvrow);
    }

    UserDataStruct {
        P,
        pivot,
        acoef: [[ZERO; NS]; NS],
        bcoef: [ZERO; NS],
        rates: N_VNew_Serial(NEQ, sunctx),
        cox: [ZERO; NS],
        coy: [ZERO; NS],
        ax: ZERO,
        ay: ZERO,
        dx: ZERO,
        dy: ZERO,
        uround: ZERO,
        sqruround: ZERO,
        mx: 0,
        my: 0,
        ns: 0,
        np: 0,
    }
}

/*
 * Load problem constants in data
 */

fn InitUserData(data: &mut UserDataStruct) {
    data.mx = MX;
    data.my = MY;
    data.ns = NUM_SPECIES;
    data.np = NUM_SPECIES / 2;
    data.ax = AX;
    data.ay = AY;
    data.dx = data.ax / (MX - 1) as f64;
    data.dy = data.ay / (MY - 1) as f64;
    data.uround = SUN_UNIT_ROUNDOFF;
    data.sqruround = data.uround.sqrt();

    /* Set up the coefficients a and b plus others found in the equations */
    let np = data.np as usize;

    let dx2 = data.dx * data.dx;
    let dy2 = data.dy * data.dy;

    for i in 0..np {
        /*  Fill in the portion of acoef in the four quadrants, row by row */
        for j in 0..np {
            data.acoef[i][np + j] = -GG; /* a1 */
            data.acoef[i + np][j] = EE; /* a2 */
            data.acoef[i][j] = ZERO; /* a3 */
            data.acoef[i + np][np + j] = ZERO; /* a4 */
        }

        /* and then change the diagonal elements of acoef to -AA */
        data.acoef[i][i] = -AA;
        data.acoef[i + np][i + np] = -AA;

        data.bcoef[i] = BB;
        data.bcoef[i + np] = -BB;

        data.cox[i] = DPREY / dx2;
        data.cox[i + np] = DPRED / dx2;

        data.coy[i] = DPREY / dy2;
        data.coy[i + np] = DPRED / dy2;
    }
}

/*
 * Set initial conditions in cc
 */

fn SetInitialProfiles(cc: &mut NVector, sc: &mut NVector) {
    let mut ctemp = [0.0f64; NS];
    let mut stemp = [0.0f64; NS];

    let ccdata = &mut cc.data;
    let scdata = &mut sc.data;

    /* Initialize arrays ctemp and stemp used in the loading process */
    for i in 0..NS / 2 {
        ctemp[i] = PREYIN;
        stemp[i] = ONE;
    }
    for i in NS / 2..NS {
        ctemp[i] = PREDIN;
        stemp[i] = 0.00001;
    }

    /* Load initial profiles into cc and sc vector from ctemp and stemp. */
    for jy in 0..MY {
        for jx in 0..MX {
            let loc = IJ_V!(jx, jy);
            for i in 0..NS {
                ccdata[loc + i] = ctemp[i];
                scdata[loc + i] = stemp[i];
            }
        }
    }
}

/*
 * Print first lines of output (problem description)
 */

fn PrintHeader(
    globalstrategy: i32,
    maxl: i32,
    maxlrst: i32,
    fnormtol: f64,
    scsteptol: f64,
    linsolver: i32,
) {
    print!("\nPredator-prey test problem --  KINSol (serial version)\n\n");
    println!("Mesh dimensions = {} X {}", MX, MY);
    println!("Number of species = {}", NUM_SPECIES);
    print!("Total system size = {}\n\n", NEQ);
    println!(
        "Flag globalstrategy = {} (0 = None, 1 = Linesearch)",
        globalstrategy
    );

    match linsolver {
        USE_SPGMR => println!(
            "Linear solver is SPGMR with maxl = {}, maxlrst = {}",
            maxl, maxlrst
        ),
        USE_SPBCGS => println!("Linear solver is SPBCGS with maxl = {}", maxl),
        USE_SPTFQMR => println!("Linear solver is SPTFQMR with maxl = {}", maxl),
        USE_SPFGMR => println!(
            "Linear solver is SPFGMR with maxl = {}, maxlrst = {}",
            maxl, maxlrst
        ),
        _ => {}
    }

    println!("Preconditioning uses interaction-only block-diagonal matrix");
    println!("Positivity constraints imposed on all components ");
    println!(
        "Tolerance parameters:  fnormtol = {}   scsteptol = {}",
        fmt_g(fnormtol, 0, 6),
        fmt_g(scsteptol, 0, 6)
    );

    println!("\nInitial profile of concentration");
    println!(
        "At all mesh points:  {} {} {}   {} {} {}",
        fmt_g(PREYIN, 0, 6),
        fmt_g(PREYIN, 0, 6),
        fmt_g(PREYIN, 0, 6),
        fmt_g(PREDIN, 0, 6),
        fmt_g(PREDIN, 0, 6),
        fmt_g(PREDIN, 0, 6)
    );
}

/*
 * Print sampled values of current cc
 */

fn PrintOutput(cc: &NVector) {
    let ccdata = &cc.data;

    let mut jy = 0;
    let mut jx = 0;
    let mut ct = IJ_V!(jx, jy);
    print!("\nAt bottom left:");

    /* Print out lines with up to 6 values per line */
    for is in 0..NS {
        if (is % 6) * 6 == is {
            println!();
        }
        print!(" {}", fmt_g(ccdata[ct + is], 0, 6));
    }

    jy = MY - 1;
    jx = MX - 1;
    ct = IJ_V!(jx, jy);
    print!("\n\nAt top right:");

    /* Print out lines with up to 6 values per line */
    for is in 0..NS {
        if (is % 6) * 6 == is {
            println!();
        }
        print!(" {}", fmt_g(ccdata[ct + is], 0, 6));
    }
    print!("\n\n");
}

/*
 * Print final statistics contained in iopt
 */

fn PrintFinalStats(kmem: &mut KINMem, linsolver: i32) {
    let (mut nni, mut nfe, mut nli, mut npe, mut nps, mut ncfl, mut nfeSG) =
        (0i64, 0i64, 0i64, 0i64, 0i64, 0i64, 0i64);

    let mut flag = KINGetNumNonlinSolvIters(kmem, &mut nni);
    check_flag(flag, "KINGetNumNonlinSolvIters");
    flag = KINGetNumFuncEvals(kmem, &mut nfe);
    check_flag(flag, "KINGetNumFuncEvals");

    flag = KINGetNumLinIters(kmem, &mut nli);
    check_flag(flag, "KINGetNumLinIters");
    flag = KINGetNumPrecEvals(kmem, &mut npe);
    check_flag(flag, "KINGetNumPrecEvals");
    flag = KINGetNumPrecSolves(kmem, &mut nps);
    check_flag(flag, "KINGetNumPrecSolves");
    flag = KINGetNumLinConvFails(kmem, &mut ncfl);
    check_flag(flag, "KINGetNumLinConvFails");
    flag = KINGetNumLinFuncEvals(kmem, &mut nfeSG);
    check_flag(flag, "KINGetNumLinFuncEvals");

    println!("Final Statistics.. ");
    println!("nni    = {:5}    nli   = {:5}", nni, nli);
    println!("nfe    = {:5}    nfeSG = {:5}", nfe, nfeSG);
    println!(
        "nps    = {:5}    npe   = {:5}     ncfl  = {:5}",
        nps, npe, ncfl
    );

    if linsolver < 3 {
        print!("\n=========================================================\n\n");
    }
}

/*
 * Check function return value (opt == 1 case of the C check_flag)
 */

fn check_flag(flag: i32, funcname: &str) -> bool {
    if flag < 0 {
        eprintln!(
            "\nSUNDIALS_ERROR: {}() failed with flag = {}\n",
            funcname, flag
        );
        return true;
    }
    false
}
