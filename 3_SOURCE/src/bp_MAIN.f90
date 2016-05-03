    program bp_main

    ! *********************************************
    ! ***     MAIN BAYES MODULE PEST PROGRAM    ***
    ! ***      Implementation of: Bayesian      ***
    ! *** Geostatistical inverse method in PEST.***
    ! ***                                       ***
    ! ***             a m!ke@usgs joint         ***
    ! ***          Michael N. Fienen, PhD       ***
    ! ***    UNITED STATES GEOLOGICAL SURVEY    ***
    ! ***          mnfienen@usgs.gov            ***
    ! ***                 AND                   ***
    ! ***           Marco D'Oria, PhD           ***
    ! ***    UNIVERSITY OF PARMA, ITALY         ***
    ! ***       marco.doria@unipr.it            ***
    ! ***           January 3, 2012             ***
    ! ***      Modified by M.D. 1/10/2009       ***
    ! ***  Further Modifications MNF/MD 11/2010 ***
    ! ***  Further Modifications   MD   09/2011 ***
    ! ***  Further Modifications MNF/MD 12/2011 ***
    ! ***   Version 1.0 --- November 19, 2012   ***
    ! *********************************************

    !When compiling the program with a new MSVS solution several settings need to be set:
    ! Errors in bso_structural_optimization. Set bgaPEST_PCGA->property->Diagnostics->Check Routine Interfaces to No. This is however not recommend: https://software.intel.com/en-us/forums/intel-visual-fortran-compiler-for-windows/topic/282188.
    ! Compiler is not happy with double declaration of e.g. iso in MIO. Set bgaPEST_PCGA->Property->Fortran->Preprocessor->Preprocess Source file: Yes (/fpp)
    ! Errors related to MKL_intel_thread.dll and MKL_core.dll:
    ! At bgaPEST_PCGA->Properties->Debugging->Environment added PATH=%IFORT_COMPILER16%redist\ia32\mkl;%IFORT_COMPILER16%redist\intel64\mkl
    ! At Computer->Properties->Advanced system settings->Environment variables->System variables->Path->Edit->Variable value added ;%IFORT_COMPILER16%redist\ia32_win\mkl;%IFORT_COMPILER16%redist\intel64_win\mkl. Solution found at https://software.intel.com/en-us/forums/intel-visual-fortran-compiler-for-windows/topic/542583
    ! MKL is activated in library
    ! Avoid rebuilding entire solution when compiling: http://stackoverflow.com/questions/1022485/why-is-visual-studio-2008-always-rebuilding-my-whole-project: Configuration Properties -> C/C++ -> General -> Debug Information Format : none (before: full). This supress the demand of recompiling every f90 file. now only the f90 files which has been changed are recompiled.
    ! Stack overflow problems: Solution properties -> Fortran -> Optimization -> Heap Array  is set to zero, which means that all arrays larger than 0 KB is send to the heap. This prevents stack overflow problems, however could lead to a slower program. This is however to prefer rather than running larger models which suddently crashes due to lack of memory. The value could maybe also be increased to make the program perform faster.
    ! It was expirenced that there was a memory leakage problem related to the read and write function in one of the first versions of the XE2016 Intel Fortran compiler. The solution was to update the compiler. Therefore, it is advised to use the latest version of the Intel Fortran compiler.

    use jupiter_input_data_support
    use bayes_pest_control
    use model_input_output
    use bayes_pest_mio_setup
    use bayes_pest_model
    use bayes_pest_reader
    use bayes_pest_finalize
    use bayes_output_control
    use error_message
    use utilities
    use make_kernels
    use bayes_matrix_operations
    use param_trans
    use linesearch
    use extern_derivs
    use bayes_output_control
    use struct_param_optimization
    use PCGA_objective_function
    use PCGA_normalize_X
    use PCGA_low_rank_app
    use PCGA_solve_linear_equations
    use PCGA_posterior_cov
    use PCGA_param_trans
    use PCGA_struct_param_optimization
    use posterior_cov_operations
    use Compare_derivative
    use mkl_dfti !MKL FFT
    use mkl_vsl  !MKL Random numbers

    implicit none

    !--  Main Data Arrays for OBS and PARS and ALGORITHM
    integer                      :: n1, i, j, imax
    integer                      :: forward_flag_der, ci95_flag
    integer                      :: s_ind, p_ind, b_ind  !Indices for structural parameters, quasi-linear and bga method loops
    type (mio_struc)             :: miostruc
    type (err_failure_struc)     :: errstruc
    integer                      :: ifail, restart, outunit, bprunit, cparunit, cobsunit, finalparunit, postcovunit
    type (cv_algorithmic)        :: cv_A
    type (d_algorithmic)         :: d_A
    type (cv_prior_mean)         :: cv_PM
    type (d_prior_mean)          :: d_PM
    type (cv_struct)             :: cv_S
    type (d_struct)              :: d_S
    type (cv_param)              :: cv_PAR
    type (Q0_compr), pointer     :: Q0_All(:)
    type (d_param)               :: d_PAR
    type (cv_observ)             :: cv_OBS
    type (d_observ)              :: d_OBS
    type (d_comlin)              :: d_MOD
    type (cv_minout)             :: cv_MIO
    type (d_minout)              :: d_MIO
    type (kernel_XQR)            :: d_XQR
    type (d_anisotropy)          :: d_ANI
    character (len=ERRORWIDTH)   :: retmsg
    character (len=100)          :: command_line, curr_par_file, curr_resid_file,post_cov_file
    character (len=FILEWIDTH)    :: ctlfile
    character (len=FILEWIDTH)    :: casename
    character (len=FILEWIDTH)    :: atemp
    character (len=20)           :: inner_iter ! aka p_ind - this is the temporary holder for printing out the inner iteration number
    character (len=20)           :: outer_iter ! aka b_ind - this is the temporary holder for printing out the outer iteration number
    double precision, dimension(1) :: curr_structural_conv, curr_phi_conv       !Current iteration convergence values for
    double precision, dimension(1) :: curr_bga_conv, curr_bga_phi,prev_bga_phi, prev_curr_phi_conv  !structural parameters and quasi linear objective function
    double precision, dimension(1) :: curr_phi !Current value for quasi linear objective function
    double precision, allocatable :: prev_struct(:) !Previous vector of theta and sigma values to be optimized for or previous objective function
    double precision, pointer    :: VV(:,:), V(:) !VV is the posterior covariance matrix, V is only the diagonal of VV
    double precision             :: structural_conv
    double precision             :: huge_val=huge(huge_val) !Largest machine number

    !PCGA
    integer :: status, p, reactivate_LM
    integer :: brng_Gaussian, seed_Gaussian, method_Gaussian
    double precision :: lambda_org, norm_diff, norm_pars, norm_oldpars

    !-----------------------------------------------

    nullify(Q0_All)
    nullify(VV)
    nullify(V)
    ! -- PRINT OUT THE BANNER INFORMATION
    call bpo_write_banner()

    !-- READ AND PARSE THE COMMAND LINE TO OBTAIN THE CONTROL FILE NAME
    call UTL_GET_COMMAND_LINE(COMMAND_LINE)
    call UTL_PARSE_COMMAND_LINE(IFAIL,COMMAND_LINE,CTLFILE,RESTART)

    !-- handle the case where no control file was indicated
    IF (IFAIL.NE.0) THEN
        call bpo_write_nocmd()
        stop
    endif

    ! -- An extension of ".bgp" is added to the bgaPEST control file if necessary.
    i=LEN_TRIM(CTLFILE)
    IF(i.GE.5) THEN
        ATEMP=CTLFILE(I-3:I)
        CALL UTL_CASETRANS(ctlfile,'lo')
        IF(ATEMP.NE.'.bgp') THEN
            CASENAME = CTLFILE
            CTLFILE(I+1:)='.bgp'
        ELSE
            CTLFILE = trim(CTLFILE)
            CASENAME = CTLFILE(1:I-4)
        ENDIF
    ELSE
        CASENAME = trim(CTLFILE)
        CTLFILE  = trim(CTLFILE) // '.bgp'
    ENDIF

    !--  INITIALIZE MIO STRUCTURE SO IT CAN BE PASSED
    if(mio_initialise(errstruc,miostruc).ne.0) then  !MD mio_initialise is an integer
        call utl_bomb_out(errstruc)
        n1=mio_finalise(errstruc,miostruc)
        stop
    endif

    ! open up the main output record file
    bprunit = utl_nextunit()
    call bpc_openfile(bprunit,trim(trim(casename) // '.bpr'),1) ![1] at end indicates open with write access

!--  READ INPUT FILE AND PERFORM ASSOCIATED ALLOCATIONS AND PARSING     
    call bpr_read(errstruc,CTLFILE,cv_A, d_A, cv_PM, d_PM, cv_S, d_S, cv_PAR,Q0_All, d_PAR, &
        cv_OBS, d_OBS, d_MOD, cv_MIO, d_MIO, d_ANI,miostruc)

    !--  SETUP THE MODEL INPUT AND OUTPUT INFORMATION (TPL/INPUT AND INS/OUTPUT PAIRINGS)
    call bpm_setup_mio(errstruc, cv_MIO, d_MIO,  cv_PAR%npargp, &
        cv_OBS%nobsgp, cv_PAR%npar, cv_OBS%nobs, miostruc)

    !--  INITIALIZE THE INVERSE MODEL
    !Make Q0, R0, X0, and InvQbb if necessary
    call bxq_make_X0_Q0_R0_InvQbb(d_PAR,cv_S,d_S,cv_PAR,d_XQR,cv_A,d_OBS,cv_OBS%nobs,d_PM,Q0_All,cv_PM,d_ANI)

    allocate(d_OBS%h(cv_OBS%nobs)) ! Allocate the current model output [y]

    !-- IF STRUCTURAL PARAMETERS WILL BE OPTIMIZED FOR, SET UP REQUIRED INFORMATION
    if ((maxval(cv_S%struct_par_opt).eq.1).or.(d_S%sig_opt.eq.1)) then
        call bxq_theta_cov_calcs(cv_PAR,cv_S,d_S,cv_PM,cv_A)
        if (cv_A%structural_conv.ge.0.) then
            allocate(prev_struct(1))
            prev_struct = UNINIT_REAL
        else
            allocate(prev_struct(cv_S%num_theta_opt))
            prev_struct = d_S%struct_par_opt_vec_0
        endif
    endif

    !-- CALL THE SETUP OF EXTERNAL DERIVATIVES FILES (IF REQUIRED).  THIS HAPPENS ONLY ONCE FOR ALL BUT PARAMETERS FILE
    if (cv_A%PCGA == 0) then
        if ((cv_A%deriv_mode .eq. 0) .or. (cv_A%deriv_mode .eq. 4)) then
            call bxd_write_ext_PEST_files(d_MOD, cv_MIO, d_MIO, cv_OBS, cv_PAR, d_OBS,cv_A)
        endif
    endif

    !-- WRITE THE HEADER INFORMATION TO THE BPR Run Record FILE
    call bpo_write_bpr_header(bprunit,casename,cv_PAR,cv_OBS,d_MOD, cv_A, &
        cv_MIO, d_MIO,Q0_all,cv_PM,d_PM,cv_S,d_S,d_PAR,d_ANI)

    !!! --- Initialize outer loop convergence values
    curr_bga_conv = huge_val ! initialize current outer loop convergence
    prev_bga_phi  = huge_val ! initialize current outer loop convergence
    curr_bga_phi  = huge_val ! initialize current bga objective function

    !-- Run forward model once to obtain the initial modeled observation vector (consistent with the initial set of parameters)
    if (cv_A%eigenvalue_test /= 1) then
        if (cv_A%PCGA == 1) write(6,*) "Performing forward computation of initial parameters"
        call bpf_model_run(errstruc, d_MOD, cv_PAR,d_PAR, cv_OBS, cv_A,  d_OBS, d_A, 0, miostruc)
    endif

    !PCGA: Normalize X -> U_X
    if(cv_A%PCGA == 1) then
        do p = 1, cv_PAR%p
            if (Q0_All(p)%Toep_flag /= 1) then
                stop "PCGA requires Toeplitz covariance. Set Toep_flag = 1"
            endif
        enddo

        write(6,'(1a1,A40,$)') char(13), "Normalize the drift matrix: In progress"
        call Normalize_X(d_XQR, cv_PAR)
        write(6,'(1a1,A40,$)') char(13), "Normalize the drift matrix: Complete   "
        write(6,*) " "

        !Settings for the random sampling in "PCGA_low_rank_app"
        call seed_random_samples(cv_A)

        !Settings for Levenberg-Marquardt
        if (cv_A%LM == 1) then
            reactivate_LM = -1
            lambda_org = cv_A%LM_lambda
        else
            cv_A%LM_lambda = 0.0
        endif

        !Derivative test
        if (cv_A%derivative_test == 1 .or. cv_A%derivative_test == 2 .or. cv_A%derivative_test == 3) then !Analyse the derivative or use PCGA bgaPEST.

            if (maxval(d_PM%Partrans).ge.1) then  !If true, the parameters transformation is required
                call PCGA_sen_par_trans(d_PAR%pars, cv_PAR, d_PAR, d_PM) !Converting sensitivity and parameters in the estimation space
                d_PAR%pars_old = d_PAR%pars
            end if

            call Compare_derivative_MFA(miostruc, errstruc, cv_A, cv_MIO, cv_S, d_S, cv_PAR, d_PAR, cv_OBS, d_MIO, d_OBS, d_MOD, d_XQR, d_PM, Q0_All)
            cv_A%it_max_bga    = 0
            cv_A%it_max_phi    = 0
            cv_A%post_cov_flag = 0
        endif
    endif

    do b_ind = 1, cv_A%it_max_bga  !*********************************************************************** (more external loop)

        !***************************************************************************************************************************
        !****************************** FROM HERE THE QUASI-LINEAR PARAMETER ESTIMATION LOOP ***************************************
        !***************************************************************************************************************************

        curr_phi_conv = huge_val !Initialize current quasi linear objective function convergence
        curr_phi      = huge_val !Initialize current quasi-linear objective function value

        if (cv_A%PCGA == 1) then !PCGA: Construct truncated U_z, C and Ap:
            write(6,'(1a1,A57,$)') char(13), "Construcing the parameter covariance matrix: In progress"
            call Low_rank_app(cv_A, cv_PAR, cv_S, d_ANI, d_PAR, d_S, d_XQR, Q0_All)
            write(6,'(1a1,A57,$)') char(13), "Construcing the parameter covariance matrix: Complete   "
            write(6,*) " "
            write(6,*) " "
            if (cv_A%eigenvalue_test == 1) exit !exit b_ind loop

            !Computing the objective function of the 'starting parameters'
            call PCGA_cal_ob_funcs(d_XQR, d_S, cv_PAR, cv_OBS, d_OBS,  cv_A, d_PAR, b_ind, 0)
            curr_phi = d_PAR%phi_T
        endif

        do p_ind = 1, cv_A%it_max_phi !************************************************************* (first intermediate loop)
            !********** quasi-liner parameter estimation for given structural parameters ***********

            !-- RUN THE FORWARD MODEL (INCLUDES DELETING OLD OUTPUT, WRITING NEW INPUT, RUNNING MODEL, AND READING NEW OUTPUT)
            if (cv_A%PCGA == 0) then !Normal bgaPEST
                select case(cv_A%deriv_mode)
                case (0)
                    forward_flag_der = 1
                case (1)
                    forward_flag_der = 2
                case (4)
                    forward_flag_der = 4
                end select
                call bpf_model_run(errstruc, d_MOD, cv_PAR,d_PAR, cv_OBS, cv_A,  d_OBS, d_A, forward_flag_der, miostruc)
            endif
            !For PCGA the derivative is computed using the matrix-free approach.

            !-- CONVERT OR NOT SENSITIVITY AND PARAMETERS IN THE ESTIMATION SPACE
            if (maxval(d_PM%Partrans).ge.1) then  !If yes, the parameters transformation is required
                select case(cv_A%PCGA)
                case(0) !Normal bgaPEST
                    d_PAR%pars_old = d_PAR%pars   !MD At the beginning pars is the vector of the initial values of the parameters
                    !as read in the parameter file. Then became the best estimate.
                    call sen_par_trans(cv_PAR, cv_OBS, d_PAR, d_A, d_PM) !Converting sensitivity and parameters in the estimation space
                case(1) !PCGA bgaPEST
                    call PCGA_sen_par_trans(d_PAR%pars, cv_PAR, d_PAR, d_PM) !Converting sensitivity and parameters in the estimation space
                    d_PAR%pars_old = d_PAR%pars !The transformed parameters are saved. Later pars will be the new estimate.
                end select
            end if

            !-- SOLVE THE BAYESIAN LINEAR SYSTEM AND CALCULATE THE OBJECTIVE FUNCTIONS
            select case(cv_A%PCGA)
            case(0) !Normal bgaPEST
                call  bmo_form_Qss_Qsy_HQsy(d_XQR, d_S%theta, cv_PAR, cv_OBS, cv_S, cv_A, d_A, d_PAR, Q0_All)
                call  bmo_form_Qyy(d_XQR, d_S%sig, cv_OBS, d_A)
                call  bmo_H_only_operations(d_XQR, d_A,cv_OBS,d_PAR,cv_PAR)
                call  bmo_solve_linear_system(d_XQR, d_S, d_PM, cv_PAR, cv_OBS, d_OBS, d_A, d_PAR,cv_PM)
            case(1) !PCGA bgaPEST
                write(6,*) '*******************************************'
                write(6,*) " "
789             FORMAT(' Iteration No.       :',I3)
                write(6,789) p_ind
                write(6,*) " "
                if (reactivate_LM == b_ind) cv_A%LM_lambda = lambda_org !Reactivating the Levenberg-Marquardt part if structural parameters has been updated.
                call PCGA_solve_linear_equation(miostruc, errstruc, cv_A, cv_MIO, d_S, cv_PAR, d_PAR, cv_OBS, d_MIO, d_OBS, d_MOD, d_XQR, d_PM)
            end select

            !-- PERFORM LINESEARCH IF REQUESTED
            if (cv_A%lns_flag.eq.1) then  !If yes, we perform the linesearch procedure
                write(6,*) " "
                write(6,*) "Performing linesearch"
                d_XQR%F_count = 0
                call lns_proc(d_XQR,d_S,cv_PAR,d_A,d_PAR,d_PM,cv_OBS,d_OBS,cv_PM,d_MOD,cv_A,p_ind,miostruc,errstruc)
                write(6,*) " "
            endif

            !Print of the new and old estimates as well as the difference in estimate
            if (cv_A%PCGA == 1) then
787             FORMAT(" Norm of estimated parameters (estimation space): ",ES11.4)
786             FORMAT(" Norm of previous  parameters (estimation space): ",ES11.4)
785             FORMAT(" Difference in norms of parameters              : ",ES11.4)
                norm_pars = dsqrt(sum((d_PAR%pars)**2))
                norm_oldpars = dsqrt(sum((d_PAR%pars_old)**2))
                norm_diff = dsqrt(sum((d_PAR%pars_old - d_PAR%pars)**2))
                write(6,787) norm_pars
                write(6,786) norm_oldpars
                write(6,785) norm_diff
                write (6,*) " "
            endif

            !-- BACK-TRANSFORM OR NOT PARAMETERS INTO PHYSICAL SPACE
            if (maxval(d_PM%Partrans).ge.1) then  !If yes, we need to back-transform the parameters in the physical space
                select case(cv_A%PCGA)
                case(0)
                    call par_back_trans(cv_PAR, d_PAR, d_PM)
                case(1)
                    call PCGA_par_back_trans(d_PAR%pars, cv_PAR, d_PAR, d_PM) !Converting sensitivity and parameters in the estimation space
                end select
            end if

            !Run the forward model to obtain the current modeled observation vector (consistent with the estimated parameters)
            !Could this one be omitted if linesearch is used? MC
            if (cv_A%PCGA == 1) write(6,*) "Performing forward computation of estimated parameters"
            call bpf_model_run(errstruc, d_MOD, cv_PAR,d_PAR, cv_OBS, cv_A,  d_OBS, d_A, 0, miostruc)

            ! UPDATE THE OBJECTIVE FUNCTION COMPONENTS
            select case(cv_A%PCGA)
            case(0)
                if (cv_A%lns_flag.eq.1) then
                    call cal_ob_funcs(d_XQR, d_S, d_PM, cv_PAR, cv_OBS, d_OBS,  d_A, d_PAR, cv_PM)
                endif
            case(1)
                call PCGA_cal_ob_funcs(d_XQR, d_S, cv_PAR, cv_OBS, d_OBS,  cv_A, d_PAR, b_ind, p_ind)
            end select

            !-- set temporary string version of iteration numbers and phi to write out
            curr_phi_conv = abs(curr_phi - d_PAR%phi_T) !Alternative convergence criterion: abs(norm_diff)/abs(norm_pars)

            !Print of the objective function
            if (cv_A%PCGA == 1) then
784             FORMAT(" Current objective function:   ",ES11.4)
783             FORMAT(" Change in objective function: ",ES11.4)
782             FORMAT(" Js: ",ES11.4)
781             FORMAT(" Jd group name ",a20,": ",ES11.4)
                write (6,*) " "
                write(*,*) "Objective function:"
                write(6,782) d_PAR%phi_R
                do i = 1,cv_OBS%nobsgp
                    write(6,781) cv_OBS%grp_name(i), d_PAR%phi_M_vec(i)
                enddo
                write(6,784) d_PAR%phi_T
                write(6,783) (d_PAR%phi_T-curr_phi)
                write (6,*) " "
            endif

            curr_phi = d_PAR%phi_T

            call UTL_INT2CHAR(p_ind,inner_iter)
            call UTL_INT2CHAR(b_ind,outer_iter)
            curr_par_file = trim(casename) // '.bpp.' // trim(outer_iter) // '_' // trim(inner_iter)
            curr_resid_file = trim(casename) // '.bre.' // trim(outer_iter) // '_' // trim(inner_iter)
            !-- Write intermediate values out to BPR record file
            call bpo_write_bpr_intermed(bprunit,p_ind,b_ind,curr_par_file,curr_resid_file,d_PAR)
            ! --Write the intermediate parameter and residuals files

            cparunit = utl_nextunit()
            call bpc_openfile(cparunit,trim(curr_par_file),1) ![1] at end indicates open with write access
            call bpo_write_allpars(cv_PAR,d_PAR,cparunit)
            close(cparunit)
            cobsunit = utl_nextunit()
            call bpc_openfile(cobsunit,trim(curr_resid_file),1) ![1] at end indicates open with write access
            call bpo_write_residuals(cv_OBS,d_OBS,cobsunit)
            close(cobsunit)

            !Levenberg-Marquardt
            if (cv_A%PCGA == 1) then
                if (cv_A%LM == 1) then !Not nessesary since when Levenberg-Marquad is not used lambda is zero
                    if (curr_phi_conv(1) < prev_curr_phi_conv(1)) then !Improvement in the objective function -> the Levenberg-Marquardt parameter is decreased
                        cv_A%LM_lambda = cv_A%LM_lambda / cv_A%LM_down
                    else !No improvement in the objective function -> the Levenberg-Marquardt parameter is increased
                        cv_A%LM_lambda = cv_A%LM_lambda * cv_A%LM_up
                    endif

                    if (cv_A%LM_lambda < cv_A%LM_min .or. p_ind >= cv_A%LM_max_it ) then !Deactivating Levenberg-Marquardt when the correction is small.
                        cv_A%LM_lambda = 0.d0
                        cv_A%LM = 0
                        reactivate_LM = b_ind+1
                        write(6,*) "Levenberg-Marquardt is deactivated since the parameter is small"
                        write(6,*) " "
                    elseif (cv_A%LM_lambda > cv_A%LM_max) then !If the Levenberg-Marquardt parameter gets to large it is reduced
                        cv_A%LM_lambda = cv_A%LM_max
                        write(6,*) "Levenberg-Marquardt parameter is redueced since it is to large"
                        write(6,*) " "
                    endif
                endif
            endif
            prev_curr_phi_conv(1) = curr_phi_conv(1)

            !-- check for convergence - exit if convergence has been achieved
            if (curr_phi_conv(1) .le. cv_A%phi_conv) then
                exit
            elseif (p_ind .ge. cv_A%it_max_phi) then
                write(retmsg,10) p_ind
10              format(' Warning: Maximum number of iterations exceeded in quasi-linear parameter optimization loop during bgaPEST iteration',i4, &
                    & '. Convergence was not achieved, but iterations will cease.')
                call utl_writmess(6,retmsg)
            endif !- checking for convergence or exceeding maximum iterations
        enddo  !(first intermediate loop) quasi-linear method  --> p_ind
        curr_bga_phi = d_PAR%phi_T  ! Set the current bga outer loop convergence to equal the current value of PHI_T
        !***************************************************************************************************************************
        !************************************* END OF QUASI-LINEAR PARAMETER ESTIMATION LOOP ***************************************
        !***************************************************************************************************************************

        !***************************************************************************************************************************
        !********************** FROM HERE THE STRUCTURAL PARAMETER ESTIMATION LOOP  (ONLY IF REQUIRED) *****************************
        !************************************************************************** *************************************************
        if ((maxval(cv_S%struct_par_opt).eq.1).or.(d_S%sig_opt.eq.1)) then !Enter the structural pars estimation loop only if required
            if (b_ind .ge. cv_A%it_max_bga) then !-- do not re-estimate structural parameters if we have exceeded maximum number of it_max_bga
                write(retmsg,11) b_ind
11              format(' Warning: Maximum number of iterations exceeded in quasi-linear parameter optimization loop during bgaPEST iteration',i4, &
                    & '. Structural parameters will not be re-calculated, so final structural parameters and posterior covariance (if requested)' &
                    & ' are based on the last iteration.')
                call utl_writmess(6,retmsg)
                write (*,*) " "
                cv_S%struct_par_opt = 1
                call bpo_write_bpr_final_structpar(bprunit,cv_S,d_S,cv_PAR,d_PAR)
                cv_S%struct_par_opt = 0
                exit
            else
                if ((cv_A%structural_conv.ge.0.).and.(b_ind.eq.1)) then !In case of objective function monitoring, here the evaluation of the objective function with the initial parameters (only at the first bga loop)
                    select case(cv_A%PCGA)
                    case(0)
                        call beg_str_object_fun(cv_OBS,d_OBS,d_A,cv_S,d_PM,cv_PAR,cv_PM)
                    case(1)
                        call PCGA_beg_str_object_fun(d_XQR,cv_OBS,d_OBS,d_A,cv_S,d_PM,cv_PAR,cv_PM,d_S)
                    end select
                    prev_struct=cv_S%str_obj_fun
                endif

                select case(cv_A%PCGA) !Here d_S%struct_par_opt_vec is the vector of the optimized theta and sigma values
                case(0)
                    call marginal_struct_param_optim(d_XQR,Q0_all,cv_OBS,d_OBS,cv_A,d_A,d_PAR,cv_S,d_S,d_PM,cv_PAR,cv_PM,b_ind,cv_S%num_theta_opt)
                case(1)
                    call PCGA_marginal_struct_param_optim(d_XQR,Q0_all,cv_OBS,d_OBS,cv_A,d_A,d_PAR,cv_S,d_S,d_PM,cv_PAR,cv_PM,errstruc,miostruc,cv_MIO, d_MIO,d_MOD,d_ANI,b_ind,cv_S%num_theta_opt)
                end select

                if (cv_A%structural_conv.ge.0.) then
                    curr_structural_conv=abs(prev_struct - cv_S%str_obj_fun) !Calculate difference between actual and previous objective function values
                    structural_conv = cv_A%structural_conv                   !The structural convergenge value remain the one assigned into .bgp file
                    prev_struct = cv_S%str_obj_fun                           !Assign the current value of the objective function to the previous
                else
                    !curr_structural_conv = sqrt(sum((prev_struct - d_S%struct_par_opt_vec)**2)) !Calculate norm of difference between actual and previous vectors
                    curr_structural_conv = sqrt(sum(((prev_struct - d_S%struct_par_opt_vec)/prev_struct)**2))!calculate the normalized root square difference
                    !between last iteration and current one
                    structural_conv = -cv_A%structural_conv             !The structural convergenge value changes sign respect to the one assigned into .bgp file
                    prev_struct = d_S%struct_par_opt_vec                !Assign the current vector of structural parameters to the previous vector
                endif

                if (curr_structural_conv(1).le.structural_conv) then !If yes, structural parameters have converged.
                    call bpo_write_bpr_intermed_structpar(bprunit,cv_S,d_S,cv_PAR,d_PAR) ! write out the final structural parameters
                    cv_S%struct_par_opt = 0  !Set to zero cv_S%struct_par_opt and d_S%sig_opt so the structural parameters estimation loop is no more entered.
                    d_S%sig_opt = 0          !The optimized struct_par_opt_vec is used to run the quasi-linear loop that is the last one.
                    if (allocated(prev_struct)) deallocate(prev_struct)
                endif
            endif !-- special warning if exceed maximum number of main algorithm iterations (it_max_bga) without convergence
            ! -- write the intermediate files to the BPR file
            call bpo_write_bpr_intermed_structpar(bprunit,cv_S,d_S,cv_PAR,d_PAR)
        else
            ! need to be sure that final values get written out.
            cv_S%struct_par_opt = 1
            call bpo_write_bpr_final_structpar(bprunit,cv_S,d_S,cv_PAR,d_PAR)
            cv_S%struct_par_opt = 0
            exit !If the structural pars optimization is not required or structural pars have converged (run the last quasi_linear), exit the bga_loop
        endif
        !***************************************************************************************************************************
        !*************************** END OF STRUCTURAL PARAMETER ESTIMATION LOOP  (ONLY IF REQUIRED) *******************************
        !***************************************************************************************************************************



        !*********************************************
        ! Evaluate outer bga convergence
        !*********************************************
        curr_bga_conv = abs(prev_bga_phi - curr_bga_phi)

        if (curr_bga_conv(1) .le. cv_A%bga_conv) then
            exit
        else
            prev_bga_phi = curr_bga_phi
        endif

    enddo      !(more external loop) --> b_ind
    ! write out a final residuals file
    cobsunit = utl_nextunit()
    curr_resid_file = trim(casename) // '.bre.fin'
    call bpc_openfile(cobsunit,trim(curr_resid_file),1) ![1] at end indicates open with write access
    call bpo_write_residuals(cv_OBS,d_OBS,cobsunit)
    close(cobsunit)

    !*************************************************************************************************************************
    !******** FROM HERE THE EVALUATION OF THE POSTERIOR COVARIANCE (ONLY IF REQUIRED --> cv_A%post_cov_flag = 1 **************
    !*********** The posterior covariance is the full matrix (matrix VV) in case of no compression of Q, *********************
    !********************* it is only the diagonal (vector V) in case of compression of Q ************************************
    !*************************************************************************************************************************
    if (cv_A%PCGA == 0 .or. cv_A%eigenvalue_test == 0) then !Not making posterior covariance when testing the eigenvalues
        if (cv_A%post_cov_flag.eq.1) then
            write(6,'(1a1,A45,$)') char(13), " Calulating Posterior Covariance: In progress"
            ci95_flag = 1 ! - flag determining whether to write 95% confidence intervals or not

            select case(cv_A%PCGA)
            case(0)
                call form_post_covariance(d_XQR, cv_PAR, cv_OBS, cv_S, cv_A, d_A, d_PAR,Q0_All,cv_PM,d_PM,d_S,VV,V)
                !-- write out the final parameter values and confidence intervals
                if (cv_A%Q_compression_flag .eq. 0) then  !Select if the Q0 matrix is compressed or not
                    allocate(V(cv_PAR%npar))
                    V = 0.D0 ! initialize the vector for V
                    do i = 1,cv_PAR%npar
                        V(i) = VV(i,i)
                    enddo
                endif
            case(1)
                cv_A%Q_compression_flag = cv_A%posterior_cov_compression_flag
                select case(cv_A%Q_compression_flag)
                case(0) !Full posterior covariance matrix
                    allocate(VV(cv_PAR%npar,cv_PAR%npar))
                    call PCGA_posterior_covariance_full(cv_A, cv_S, cv_PAR, Q0_All, cv_OBS, d_XQR, d_S, VV)
                    allocate(V(cv_PAR%npar))
                    do i = 1,cv_PAR%npar
                        V(i) = VV(i,i) !bpp can only contain the diagonal
                    enddo
                case(1) !Diagonal of posterior covariance matrix
                    allocate(V(cv_PAR%npar))
                    do i = 1,cv_PAR%npar
                        call PCGA_posterior_covariance_point(cv_A, cv_S, cv_PAR, Q0_All, cv_OBS, d_XQR, d_S, i, i, V(i))
                    enddo

                    case default
                    write(*,*) "Incorrect Q_compression_flag value. Printed posterior covariance is incorrect."
                end select
            end select

            finalparunit = utl_nextunit()
            curr_par_file = trim(casename) // '.bpp.fin'
            call bpc_openfile(finalparunit,trim(curr_par_file),1) ![1] at end indicates open with write access
            call bpo_write_allpars_95ci(cv_PAR,d_PAR,d_PM,V,finalparunit,ci95_flag)
            close(finalparunit)
            ! --- Also write a separate file with only the posterior covariance values
            postcovunit = utl_nextunit()
            post_cov_file = trim(casename) // '.post.cov'
            call bpc_openfile(postcovunit,trim(post_cov_file),1) ![1] at end indicates open with write access
            if (cv_A%Q_compression_flag.eq.0) then
                if (associated(V)) deallocate(V)
            endif
            call  bpo_write_posterior_covariance(cv_A%Q_compression_flag,cv_PAR,d_PAR,d_PM,V,VV,postcovunit)
            if (associated(VV)) deallocate(VV)
            if (associated(V)) deallocate(V)
            close(postcovunit)
            write(6,'(1a1,A45,$)') char(13), " Calulating Posterior Covariance: Complete   "
            write(6,*) " "
        else ! still need to write the final parameters out
            ci95_flag = 0 ! - flag determining whether to write 95% confidence intervals or not
            allocate(V(1))
            V = 0.
            finalparunit = utl_nextunit()
            curr_par_file = trim(casename) // '.bpp.fin'
            call bpc_openfile(finalparunit,trim(curr_par_file),1) ![1] at end indicates open with write access
            call bpo_write_allpars_95ci(cv_PAR,d_PAR,d_PM,V,finalparunit,ci95_flag)
            close(finalparunit)

            if (associated(V)) deallocate(V)
        endif
    endif
    !*************************************************************************************************************************
    !*********** END OF THE EVALUATION OF THE POSTERIOR COVARIANCE (ONLY IF REQUIRED --> cv_A%post_cov_flag = 1 **************
    !*************************************************************************************************************************

    if (cv_A%PCGA == 1) write(*,*) 'Memory cleaning up'

    !-- FINALIZE and CLEANUP - deallocate all the structures
    call bpd_finalize(errstruc, miostruc, Q0_All, d_PM, d_S, cv_A, cv_S, cv_PAR, d_PAR, cv_OBS, d_OBS, d_MOD, d_MIO, d_XQR)

    if (cv_A%PCGA == 1) write(*,*) 'Parameter estimation is complete!'
    
    end program bp_main