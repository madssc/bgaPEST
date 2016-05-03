    module PCGA_solve_linear_equations

    use bayes_pest_control
    use model_input_output
    use utilities
    use error_message
    use jupiter_input_data_support
    use PCGA_param_trans

    contains
    
    !******************************************************************************
    ! This function solves the cokriging equations and computes the new estimate
    !******************************************************************************
    subroutine PCGA_solve_linear_equation(miostruc, errstruc, cv_A, cv_MIO, d_S, cv_PAR, d_PAR, cv_OBS, d_MIO, d_OBS, d_MOD, d_XQR, d_PM)

    implicit none
    type (mio_struc)         , intent(inout)        :: miostruc
    type (err_failure_struc) , intent(inout)        :: errstruc
    type (cv_algorithmic)    , intent(inout)        :: cv_A
    type (cv_minout)         , intent(in)           :: cv_MIO
    type (d_struct)          , intent(inout)        :: d_S
    type (cv_param)          , intent(inout)        :: cv_PAR
    type (d_param)           , intent(inout)        :: d_PAR
    type (cv_observ)         , intent(inout)        :: cv_OBS
    type (d_minout)          , intent(in)           :: d_MIO
    type (d_observ)          , intent(inout)        :: d_OBS
    type (d_comlin)          , intent(inout)        :: d_MOD
    type (kernel_XQR)        , intent(inout)        :: d_XQR
    type (d_prior_mean)      , intent(in)           :: d_PM
    integer 							            :: i, j, ierr, iii, jjj, F_count, work_Z_size, info_sysv
    integer, parameter                              :: out_unit=20
    integer, allocatable                            :: ipiv(:)
    logical, save                                   :: first_entry = .true.
    double precision, allocatable                   :: work_sysv(:)

    !Nullify
    if (first_entry) then
        nullify(d_XQR%HX)
        nullify(d_XQR%hs0Hs0)
        nullify(d_XQR%B)
        nullify(d_XQR%BC)
        nullify(d_XQR%HQH)
        nullify(d_XQR%LApy)
        nullify(d_XQR%LHS)
        nullify(d_XQR%HQ)
        nullify(d_XQR%HQAp)
        nullify(d_XQR%XA_p)
        nullify(d_XQR%RHS)
        nullify(d_XQR%C_S)
        nullify(d_XQR%LA_p)
        nullify(d_XQR%MA_p)
        first_entry = .false.
    endif

    if (cv_A%LM == 1) then
789     FORMAT(" Current Levenberg-Marquardt parameter:",F7.2)
        write(6,789) cv_A%LM_lambda
    endif

    !Performing the matrix-free computations
    call matrix_free_approach_computations(miostruc,errstruc,cv_A,cv_MIO,cv_PAR,d_PAR,cv_OBS,d_MIO,d_OBS,d_MOD,d_XQR,d_PM,1)

    if (associated(d_XQR%HQ)) deallocate(d_XQR%HQ)
    if (associated(d_XQR%LA_p)) deallocate(d_XQR%LA_p)
    if (associated(d_XQR%MA_p)) deallocate(d_XQR%MA_p)

    !Compute BC
    if (associated(d_XQR%BC)) deallocate(d_XQR%BC)
    allocate(d_XQR%BC(cv_OBS%nobs,cv_A%K))
    !dgemm ( transa='n' ,  transb='n' ,  m=cv_OBS%nobs ,  n=cv_A%K ,  k=cv_A%K ,  alpha=1.d0 ,  a=d_XQR%B ,  lda=cv_OBS%nobs ,  b=d_XQR%tC ,  ldb=cv_A%K ,  beta=0.d0 ,  c=d_XQR%BC ,  ldc=cv_OBS%nobs )
    call dgemm ( 'n' ,  'n' ,  cv_OBS%nobs ,  cv_A%K ,  cv_A%K ,  1.d0 ,  d_XQR%B ,  cv_OBS%nobs ,  d_XQR%tC ,  cv_A%K ,  0.d0 ,  d_XQR%BC ,  cv_OBS%nobs )

    !Compute HQH^T = BC*B^T
    if (associated(d_XQR%HQH)) deallocate(d_XQR%HQH)
    allocate(d_XQR%HQH(cv_OBS%nobs,cv_OBS%nobs))
    !dgemm ( transa='n' ,  transb='t' ,  m=cv_OBS%nobs ,  n=cv_OBS%nobs ,  k=cv_A%K ,  alpha=1.d0 ,  a=d_XQR%BC ,  lda=cv_OBS%nobs ,  b=B ,  ldb=cv_OBS%nobs ,  beta=0.d0 ,  c=d_XQR%HQH ,  ldc=cv_OBS%nobs )
    call dgemm ( 'n' ,  't' ,  cv_OBS%nobs ,  cv_OBS%nobs ,  cv_A%K ,  1.d0 ,  d_XQR%BC ,  cv_OBS%nobs ,  d_XQR%B ,  cv_OBS%nobs ,  0.d0 ,  d_XQR%HQH ,  cv_OBS%nobs )

    !Construct LHS
    call construct_LHS(d_S, cv_A, cv_PAR, cv_OBS, d_XQR)

    !Construct RHS
    call construct_RHS(cv_A, cv_PAR, cv_OBS, d_XQR)

    !write out some of the computations. Should be removed in the final version. However, it is very usefull for debugging and also debugging new models for errors
    open (unit=out_unit,file="results_bga.txt",action="write",status="replace")
    write (out_unit,*) "tU_Z = ["
    do jjj = 1,cv_PAR%npar
        write (out_unit,'(*(D30.20))') ((d_XQR%tU_Z(jjj,:)))
    enddo
    write (out_unit,*) "];"
    write (out_unit,*) "tC = ["
    do jjj = 1,cv_A%K
        write (out_unit,'(*(D30.20))') ((d_XQR%tC(jjj,:)))
    enddo
    write (out_unit,*) "];"
    write (out_unit,*) "B = ["
    do jjj = 1,cv_OBS%nobs
        write (out_unit,'(*(D30.20))') ((d_XQR%B(jjj,:)))
    enddo
    write (out_unit,*) "];"
    write (out_unit,*) "LHS = ["
    do jjj = 1,cv_OBS%nobs+cv_PAR%p
        write (out_unit,'(*(D30.20))') ((d_XQR%LHS(jjj,:)))
    enddo
    write (out_unit,*) "];"

    write (out_unit,*) "RHS = ["
    do jjj = 1,cv_OBS%nobs+cv_PAR%p
        write (out_unit,'(*(D30.20))') ((d_XQR%RHS(jjj,:)))
    enddo
    write (out_unit,*) "];"

    !Normalization of LHS and RHS
    call normalize_LHS_RHS(cv_PAR, cv_OBS, d_XQR)

    write (out_unit,*) "LHS_Normalized = ["
    do jjj = 1,cv_OBS%nobs+cv_PAR%p
        write (out_unit,'(*(D30.20))') ((d_XQR%LHS(jjj,:)))
    enddo
    write (out_unit,*) "];"

    write (out_unit,*) "RHS_Normalized = ["
    do jjj = 1,cv_OBS%nobs+cv_PAR%p
        write (out_unit,'(*(D30.20))') ((d_XQR%RHS(jjj,:)))
    enddo
    write (out_unit,*) "];"

    !Solve linear set of equations
    !Internal function: SLVSSU(N,nrhs,A,X) -- RHS becomes x
    !call SLVSSU(cv_OBS%nobs+cv_PAR%p,cv_A%K+cv_PAR%p,d_XQR%LHS,d_XQR%RHS)
    allocate(ipiv(cv_OBS%nobs+cv_PAR%p))
    allocate(work_sysv(1))
    !call dsysv (  uplo='U' ,  n=cv_OBS%nobs+cv_PAR%p ,  nrhs=cv_PAR%p+cv_A%K ,  a=d_XQR%LHS ,  lda=cv_OBS%nobs+cv_PAR%p ,  ipiv=ipiv ,  b=d_XQR%RHS ,  ldb=cv_OBS%nobs+cv_PAR%p ,  work=work_sysv ,  lwork=-1 ,  info=info_sysv )
    call dsysv (  'U' ,  cv_OBS%nobs+cv_PAR%p ,  cv_PAR%p+cv_A%K ,  d_XQR%LHS ,  cv_OBS%nobs+cv_PAR%p ,  ipiv ,  d_XQR%RHS ,  cv_OBS%nobs+cv_PAR%p ,  work_sysv ,  -1 ,  info_sysv )
    work_Z_size = int(work_sysv(1))
    deallocate(work_sysv)
    allocate(work_sysv(work_Z_size))
    !call dsysv (  uplo='U' ,  n=cv_OBS%nobs+cv_PAR%p ,  nrhs=cv_PAR%p+cv_A%K ,  a=d_XQR%LHS ,  lda=cv_OBS%nobs+cv_PAR%p ,  ipiv=ipiv ,  b=d_XQR%RHS ,  ldb=cv_OBS%nobs+cv_PAR%p ,  work=work_sysv ,  lwork=work_Z_size ,  info=info_sysv )
    call dsysv (  'U' ,  cv_OBS%nobs+cv_PAR%p ,  cv_PAR%p+cv_A%K ,  d_XQR%LHS ,  cv_OBS%nobs+cv_PAR%p ,  ipiv ,  d_XQR%RHS ,  cv_OBS%nobs+cv_PAR%p ,  work_sysv ,  work_Z_size ,  info_sysv )
    deallocate(ipiv)
    deallocate(work_sysv)

    write (out_unit,*) "res_inv = ["
    do jjj = 1,cv_OBS%nobs+cv_PAR%p
        write (out_unit,'(*(D30.20))') ((d_XQR%RHS(jjj,:)))
    enddo
    write (out_unit,*) "];"

    !Rescale solutions
    call rescale_solution(cv_A, cv_PAR, cv_OBS, d_XQR)

    write (out_unit,*) "LA_p = ["
    do jjj = 1,cv_OBS%nobs
        write (out_unit,'(*(D30.20))') ((d_XQR%LA_p(jjj,:)))
    enddo
    write (out_unit,*) "];"

    write (out_unit,*) "MA_p = ["
    do jjj = 1,cv_PAR%p
        write (out_unit,'(*(D30.20))') ((d_XQR%MA_p(jjj,:)))
    enddo
    write (out_unit,*) "];"

    write (out_unit,*) "C_S = ["
    do jjj = 1,cv_OBS%nobs+cv_PAR%p
        write (out_unit,'(*(D30.20))') ((d_XQR%C_S(jjj)))
    enddo
    write (out_unit,*) "];"

    !Compute the estimate of the unknown parameters:
    
    !LApy = LAp*(y-hs0+Hs0)
    if (associated(d_XQR%LApy))   deallocate(d_XQR%LApy,stat=ierr)
    allocate(d_XQR%LApy(cv_A%K+cv_PAR%p))
    d_XQR%LApy = 0.d0
    !dgemv ( trans='t' ,  m=(cv_A%K+cv_PAR%p) ,  n=cv_OBS%nobs ,  alpha=1.d0 ,  a=d_XQR%LA_p ,  lda=(cv_A%K+cv_PAR%p) ,  x=(d_OBS%obs-hs0Hs0) ,  incx=1 ,  beta=0.d0 ,  y=LApy ,  incy=1 )
    call dgemv ( 't' ,  cv_OBS%nobs , (cv_A%K+cv_PAR%p) ,  1.d0 ,  d_XQR%LA_p , cv_OBS%nobs , (d_OBS%obs+d_XQR%hs0Hs0)  ,  1 ,  0.d0 ,  d_XQR%LApy ,  1 )

    !Estimate = Ap*LApy
    !dgemv ( trans='n' ,  m=cv_PAR%npar ,  n=(cv_A%K+cv_PAR%p) ,  alpha=1.d0 ,  a=d_XQR%Ap ,  lda=cv_PAR%npar ,  x=LApy ,  incx=1 ,  beta=0.d0 ,  y=d_XQR%s_hat ,  incy=1 )
    call dgemv ( 'n' ,  cv_PAR%npar ,  (cv_A%K+cv_PAR%p) ,  1.d0 ,  d_XQR%Ap ,  cv_PAR%npar ,  d_XQR%LApy ,  1 ,  0.d0 ,  d_PAR%pars ,  1 )
    
    write (out_unit,*) "Ap = ["
    do jjj = 1,cv_PAR%npar
        write (out_unit,'(*(D30.20))') ((d_XQR%Ap(jjj,:)))
    enddo
    write (out_unit,*) "];"
    write (out_unit,*) "obs = ["
    do jjj = 1,cv_OBS%nobs
        write (out_unit,'(*(D30.20))') ((d_OBS%obs(jjj)))
    enddo
    write (out_unit,*) "];"
    write (out_unit,*) "hs0Hs0 = ["
    do jjj = 1,cv_OBS%nobs
        write (out_unit,'(*(D30.20))') ((d_XQR%hs0Hs0(jjj)))
    enddo
    write (out_unit,*) "];"
    write (out_unit,*) "s_new = ["
    do jjj = 1,cv_PAR%npar
        write (out_unit,'(*(D30.20))') ((d_PAR%pars(jjj)))
    enddo
    write (out_unit,*) "];"
    close (out_unit)

    deallocate(d_XQR%B)
    deallocate(d_XQR%C_S)
    deallocate(d_XQR%HX)
    deallocate(d_XQR%LHS)
    deallocate(d_XQR%RHS)
    if (associated(d_XQR%HQAp)) deallocate(d_XQR%HQAp)
    if (associated(d_XQR%XA_p)) deallocate(d_XQR%XA_p)

    end subroutine PCGA_solve_linear_equation

    !****************************************************************************************************
    ! This functions computes all the Jacobian-vector products using the matrix-free approach function
    !****************************************************************************************************
    subroutine matrix_free_approach_computations(miostruc,errstruc,cv_A,cv_MIO,cv_PAR,d_PAR,cv_OBS,d_MIO,d_OBS,d_MOD,d_XQR,d_PM,flag)

    type (mio_struc)         , intent(inout)        :: miostruc
    type (err_failure_struc) , intent(inout)        :: errstruc
    type (cv_algorithmic)    , intent(inout)        :: cv_A
    type (cv_minout)         , intent(in)           :: cv_MIO
    type (cv_param)          , intent(inout)        :: cv_PAR
    type (d_param)           , intent(inout)        :: d_PAR
    type (cv_observ)         , intent(inout)        :: cv_OBS
    type (d_minout)          , intent(in)           :: d_MIO
    type (d_observ)          , intent(inout)        :: d_OBS
    type (d_comlin)          , intent(inout)        :: d_MOD
    type (kernel_XQR)        , intent(inout)        :: d_XQR
    type (d_prior_mean)      , intent(in)           :: d_PM
    integer                  , intent(in)           :: flag !1 = PCGA_solve_linear_equations, 2 = PCGA_structural_optimization
    integer 							            :: i, F_count, F_count_tot

    if (cv_A%Nthread == 0) then
        if (cv_A%deriv_type == 1 .or. cv_A%deriv_type == 2) then
            F_count_tot = 1+cv_A%K*cv_A%deriv_accuracy
        elseif (cv_A%deriv_type == 3) then
            F_count_tot = 1+2*cv_A%K*cv_A%deriv_accuracy
        endif
        if (flag == 1) F_count_tot = F_count_tot + cv_PAR%p

        !Compute -h(s0) + Hs0
        write(6,788) char(13), 1, F_count_tot
        if (associated(d_XQR%hs0Hs0)) deallocate(d_XQR%hs0Hs0)
        allocate(d_XQR%hs0Hs0(cv_OBS%nobs))
        call matrix_free_approach(errstruc, miostruc, cv_A, cv_MIO, cv_OBS, cv_PAR, d_MIO, d_MOD, d_PAR, d_OBS, d_PM, 2, d_PAR%pars, d_XQR%hs0Hs0)

        !Compute B = H*tU_Z
        allocate(d_XQR%B(cv_OBS%nobs,cv_A%K))
        do i=1,cv_A%K
            if (cv_A%deriv_type < 3) then
                F_count = 1+i
            else
                F_count = 1+2*i
            endif
            write(6,788) char(13), F_count, F_count_tot
            call matrix_free_approach(errstruc, miostruc, cv_A, cv_MIO, cv_OBS, cv_PAR, d_MIO, d_MOD, d_PAR, d_OBS, d_PM, 1, d_XQR%tU_Z(:,i), d_XQR%B(:,i))
        end do

        !Compute HX, only for solve linear system of equations
        if (flag == 1) then
            allocate(d_XQR%HX(cv_OBS%nobs,cv_PAR%p))
            do i=1,cv_PAR%p
                if (cv_A%deriv_type < 3) then
                    F_count = i+1+cv_A%K
                else
                    F_count = i+1+2*cv_A%K
                endif
                write(6,788) char(13), F_count, F_count_tot
                !matrix_free_approach(errstruc, miostruc, d_MOD, d_PAR, d_OBS, flag, a, Ha)
                call matrix_free_approach(errstruc, miostruc, cv_A, cv_MIO,cv_OBS, cv_PAR, d_MIO, d_MOD, d_PAR, d_OBS, d_PM, 1, d_XQR%U_X(:,i), d_XQR%HX(:,i))
            end do
        endif

        write(6,*) " "
    elseif (cv_A%Nthread > 0) then
        write (*,*) "Multi thread is not implemented yet. bgaPEST cannot continue."
        stop
    else
        write (*,*) "Nthread has to be zero or positive. bgaPEST cannot continue."
        stop
    endif

    !Formats
788 FORMAT(1a1,' Performing forward computation No.',I3,' /',I3,$)

    end subroutine matrix_free_approach_computations

    !*************************************************
    ! Computing H*a using the matrix free approach
    !*************************************************
    subroutine matrix_free_approach(errstruc, miostruc, cv_A, cv_MIO, cv_OBS, cv_PAR, d_MIO, d_MOD, d_PAR, d_OBS, d_PM, flag, a, Ha)

    implicit none
    !--  Main Data Arrays for OBS and PARS
    type (mio_struc)                        :: miostruc
    type (err_failure_struc)                :: errstruc
    type (cv_algorithmic), intent(in)       :: cv_A
    type (cv_minout),      intent(in)       :: cv_MIO
    type (cv_observ),      intent(in)       :: cv_OBS
    type (cv_param),       intent(in)       :: cv_PAR
    type (d_minout),       intent(in)       :: d_MIO
    type (d_comlin),       intent(in)       :: d_MOD
    type (d_param),        intent(inout)    :: d_PAR
    type (d_observ),       intent(in)       :: d_OBS
    type (d_prior_mean),   intent(in)       :: d_PM
    integer,               intent(in)       :: flag
    double precision,      intent(in)       :: a(:)
    double precision,      intent(inout)    :: Ha(:)
    real(kind=selected_real_kind(30)), allocatable  :: Ha_tmp(:)
    double precision, allocatable                   :: pars_perturbed(:), h_perturbed(:,:)
    double precision                                :: Norm_a, Norm_s0, sign, forward_diff_coeff(6,7), central_diff_coeff(4,-4:4)
    integer                                         :: i, numinfile, numinfile_index

    !-- MIO delete the last set of output files
    if(mio_delete_model_output_files(errstruc,miostruc).ne.0) then
        call utl_bomb_out(errstruc)
    endif

    !Check if the parameters deriv_type and deriv_accuracy are set correctly
    if (cv_A%deriv_type == 1 .or. cv_A%deriv_type == 2) then !forward/backward
        if (cv_A%deriv_accuracy < 1 .or. cv_A%deriv_accuracy > 6) stop "Accuracy for the derivative is set incorrectly. It should be between 1 and 6 for forward/backward finite difference"
    elseif (cv_A%deriv_type == 3) then
        if (cv_A%deriv_accuracy < 1 .or. cv_A%deriv_accuracy > 4) stop "Accuracy for the derivative is set incorrectly. It should be between 1 and 4 for central finite difference"
    else
        stop "deriv_type is set incorrectly. It should be 1, 2 or 3."
    endif

    !Set the coefficients for forward/backward/centreal finite difference:
    forward_diff_coeff(1,:) = (/         -1.d0, 1.d0,        0.d0,       0.d0,        0.d0,       0.d0,       0.d0 /)
    forward_diff_coeff(2,:) = (/    -3.d0/2.d0, 2.d0,  -1.d0/2.d0,       0.d0,        0.d0,       0.d0,       0.d0 /)
    forward_diff_coeff(3,:) = (/   -11.d0/6.d0, 3.d0,  -3.d0/2.d0,  1.d0/3.d0,        0.d0,       0.d0,       0.d0 /)
    forward_diff_coeff(4,:) = (/  -25.d0/12.d0, 4.d0,       -3.d0,  4.d0/3.d0,  -1.d0/4.d0,       0.d0,       0.d0 /)
    forward_diff_coeff(5,:) = (/ -137.d0/60.d0, 5.d0,       -5.d0, 10.d0/3.d0,  -5.d0/4.d0,  1.d0/5.d0,       0.d0 /)
    forward_diff_coeff(6,:) = (/  -49.d0/20.d0, 6.d0, -15.d0/2.d0, 20.d0/3.d0, -15.d0/4.d0,  6.d0/5.d0, -1.d0/6.d0 /)

    sign = 1.d0
    if (cv_A%deriv_type == 2) then !backward
        forward_diff_coeff = -forward_diff_coeff
        sign = -1.d0
    endif

    central_diff_coeff(1,:) = (/        0.d0,         0.d0,        0.d0, -1.d0/2.d0, 0.d0, 1.d0/2.d0,        0.d0,        0.d0,         0.d0 /)
    central_diff_coeff(2,:) = (/        0.d0,         0.d0,  1.d0/12.d0, -2.d0/3.d0, 0.d0, 2.d0/3.d0, -1.d0/12.d0,        0.d0,         0.d0 /)
    central_diff_coeff(3,:) = (/        0.d0,  -1.d0/60.d0,  3.d0/20.d0, -3.d0/4.d0, 0.d0, 3.d0/4.d0, -3.d0/20.d0,  1.d0/60.d0,         0.d0 /)
    central_diff_coeff(4,:) = (/ 1.d0/280.d0, -4.d0/105.d0,   1.d0/5.d0, -4.d0/5.d0, 0.d0, 4.d0/5.d0,  -1.d0/5.d0, 4.d0/105.d0, -1.d0/280.d0 /)

    !Allocate for perturbed parameters and corresponding forward computations.
    allocate(pars_perturbed(cv_PAR%npar))
    if (cv_A%deriv_type == 1 .or. cv_A%deriv_type == 2) then ! .and. flag == 1
        allocate(h_perturbed(cv_OBS%nobs, cv_A%deriv_accuracy))
    elseif (cv_A%deriv_type == 3) then ! .and. flag == 1
        allocate(h_perturbed(cv_OBS%nobs, 2*cv_A%deriv_accuracy+1)) !One of them are unused, however it is possible easier to read afterwards compared to a fix of the problem. The number of obs are probably not that big (e.g. 1000) wherefore it does not matter to allocate to much.
        !    elseif (flag == 2) then
        !allocate(h_perturbed(cv_OBS%nobs, 1))
    endif

    !Compute the norm
    if (flag == 1) then
        Norm_a  = sqrt(sum(a**2)) !Norm2 a
        !Norm_s0 = sqrt(sum(d_PAR%pars**2)) !Norm2 s0
        Norm_s0 = 0.d0
        do i = 1,cv_PAR%npar
            if (abs(a(i)) > 1.d-14) Norm_s0 = Norm_s0 + d_PAR%pars(i)**2
        enddo
        Norm_s0 = sqrt(norm_s0)
        
    elseif (flag == 2) then
        Norm_a  = 1.0
        Norm_s0 = 1.0
    endif

    numinfile = cv_MIO%ntplfle  !mio_get_numinfile(miostruc)
    do i = -cv_A%deriv_accuracy, cv_A%deriv_accuracy

        !non-used forward computations is skipped
        if (cv_A%deriv_type <= 2 .and. i < 1) then
            cycle
        elseif (cv_A%deriv_type == 3 .and. i == 0) then
            cycle
        endif

        do numinfile_index = 1, numinfile

            !setting the perturbed parameters. hmm gør det omvendt
            pars_perturbed = d_PAR%pars + (sign * i * (Norm_s0 * d_MIO%derinc_PCGA(numinfile_index) / Norm_a)) * a

            !Transforming the values to physical space
            if (maxval(d_PM%Partrans).ge.1) then
                !_matrix_free
                call PCGA_par_back_trans_matrix_free(pars_perturbed, cv_PAR, d_PAR, d_PM) !Converting sensitivity and parameters in the estimation space
            endif

            !-- MIO write the model input files
            if(mio_write_model_input_files_mult_delta(errstruc,miostruc, pars_perturbed, numinfile_index).ne.0) then
                call utl_bomb_out(errstruc)
            endif
        enddo !numinfile_index

        !Not nessesary
        !Transforming the values back to estimation space
        !call PCGA_sen_par_trans(pars_perturbed, cv_PAR, d_PAR, d_PM) !Converting sensitivity and parameters in the estimation space

        !-- Matrix-free Approach
        call system(d_MOD%com)

        !-- MIO read the ouput file results and update
        !if (flag ==1) then
        if (cv_A%deriv_type <= 2) then
            if(mio_read_model_output_files(errstruc,miostruc, h_perturbed(:,i)).ne.0) then
                call utl_bomb_out(errstruc)
            endif
        elseif (cv_A%deriv_type == 3) then
            if(mio_read_model_output_files(errstruc,miostruc, h_perturbed(:,i+cv_A%deriv_accuracy+1)).ne.0) then
                call utl_bomb_out(errstruc)
            endif
        endif

    enddo !i

    allocate(Ha_tmp(cv_OBS%nobs))

    !Adding the forward computations, multiplied by different coefficients
    if (cv_A%deriv_type <= 2) then

        Ha_tmp = forward_diff_coeff(cv_A%deriv_accuracy,1) * d_OBS%h
        do i = 1, cv_A%deriv_accuracy
            Ha_tmp = Ha_tmp + forward_diff_coeff(cv_A%deriv_accuracy,i+1) * h_perturbed(:,i)
        enddo

    elseif (cv_A%deriv_type == 3) then

        Ha_tmp = 0.d0
        do i = -cv_A%deriv_accuracy, cv_A%deriv_accuracy
            if (i /= 0) then
                Ha_tmp = Ha_tmp + central_diff_coeff(cv_A%deriv_accuracy,i) * h_perturbed(:,i+cv_A%deriv_accuracy+1)
            endif
        enddo

    endif

    !Normalization
    Ha = (Norm_a / Norm_s0) * dble(Ha_tmp)

    deallocate(Ha_tmp)

    !Dividing by delta
    do i = 1, cv_OBS%nobs
        Ha(i) = Ha(i) / d_MIO%derinc_PCGA(d_OBS%TemplateFileNo(i))
    enddo

    if (flag == 2) then
        Ha = Ha - d_OBS%h
    endif

    deallocate(h_perturbed)
    deallocate(pars_perturbed)

    end subroutine matrix_free_approach

    !*************************************************
    ! Constructing the left hand side matrix (LHS)
    !*************************************************
    subroutine construct_LHS(d_S, cv_A, cv_PAR, cv_OBS, d_XQR)

    type (d_struct)         , intent(in)    :: d_S
    type (cv_algorithmic)   , intent(in)    :: cv_A
    type (cv_param)         , intent(in)    :: cv_PAR
    type (cv_observ)        , intent(in)    :: cv_OBS
    type (kernel_XQR)       , intent(inout) :: d_XQR
    integer                                 :: i

    allocate(d_XQR%LHS((cv_OBS%nobs+cv_PAR%p),(cv_OBS%nobs+cv_PAR%p)))
    d_XQR%LHS = 0.d0

    !Incert HQH
    d_XQR%LHS(:cv_OBS%nobs,:cv_OBS%nobs) = d_XQR%HQH(:,:)
    
    !Adding R
    do i = 1, cv_OBS%nobs
        d_XQR%LHS(i,i) = d_XQR%LHS(i,i) + ((1.d0 + cv_A%LM_lambda) * d_S%sig * d_XQR%R0(i,i)) !adding R, diagonal.
    end do
    
    !Incerting HX and HX^T
    d_XQR%LHS(:cv_OBS%nobs,(cv_OBS%nobs+1):) = d_XQR%HX(:,:)
    d_XQR%LHS((cv_OBS%nobs+1):,:cv_OBS%nobs) = transpose(d_XQR%HX(:,:))
    
    !Zeros in the right lower corner
    d_XQR%LHS((cv_OBS%nobs+1):,(cv_OBS%nobs+1):) = 0.d0

    end subroutine construct_LHS

    !*************************************************
    ! Constructing the right hand side (RHS)
    !*************************************************
    subroutine construct_RHS(cv_A, cv_PAR, cv_OBS, d_XQR)

    type (cv_algorithmic)   , intent(in)    :: cv_A
    type (cv_param)         , intent(in)    :: cv_PAR
    type (cv_observ)        , intent(in)    :: cv_OBS
    type (kernel_XQR)       , intent(inout) :: d_XQR
    integer                                 :: i

    !Original way of computing the RHS
    !allocate(d_XQR%HQ(cv_OBS%nobs,cv_PAR%npar))
    !allocate(d_XQR%HQAp(cv_OBS%nobs,(cv_PAR%p+cv_A%K)))
    !allocate(d_XQR%XA_p(cv_PAR%p,(cv_PAR%p+cv_A%K)))
    !!--HQ = BC*A^T = BC*tU_Z^T
    !!dgemm ( transa='n' ,  transb='t' ,  m=cv_OBS%nobs ,  n=cv_PAR%npar ,  k=cv_A%K ,  alpha=1.d0 ,  a=d_XQR%BC ,  lda=cv_OBS%nobs ,  b=d_XQR%tU_Z ,  ldb=cv_PAR%npar ,  beta=0.d0 ,  c=d_XQR%HQ ,  ldc=cv_OBS%nobs )
    !call dgemm ( 'n' ,  't' ,  cv_OBS%nobs ,  cv_PAR%npar ,  cv_A%K ,  1.d0 ,  d_XQR%BC ,  cv_OBS%nobs ,  d_XQR%tU_Z ,  cv_PAR%npar ,  0.d0 ,  d_XQR%HQ ,  cv_OBS%nobs )
    !!--Compute HQ*Ap and X^T*Ap
    !!dgemm ( transa='n' ,  transb='n' ,  m=cv_OBS%nobs ,  n=cv_A%K+cv_PAR%p ,  k=cv_PAR%npar ,  alpha=1.d0 ,  a=d_XQR%HQ ,  lda=cv_OBS%nobs ,  b=d_XQR%Ap ,  ldb=cv_PAR%npar ,  beta=0.d0 ,  c=d_XQR%HQAp ,  ldc=cv_OBS%nobs )
    !call dgemm ( 'n' ,  'n' ,  cv_OBS%nobs ,  cv_A%K+cv_PAR%p ,  cv_PAR%npar ,  1.d0 ,  d_XQR%HQ ,  cv_OBS%nobs ,  d_XQR%Ap ,  cv_PAR%npar , 0.d0 ,  d_XQR%HQAp ,  cv_OBS%nobs )
    !!dgemm ( transa='t' ,  transb='n' ,  m=cv_PAR%p ,  n=cv_A%K+cv_PAR%p ,  k=cv_PAR%npar ,  alpha=1.d0 ,  a=d_XQR%U_X ,  lda=cv_PAR%npar ,  b=d_XQR%Ap ,  ldb=cv_PAR%npar ,  beta=0.d0 ,  c=d_XQR%XA_p ,  ldc=cv_PAR%p )
    !call dgemm ( 't' ,  'n' ,  cv_PAR%p ,  cv_A%K+cv_PAR%p ,  cv_PAR%npar ,  1.d0 ,  d_XQR%U_X ,  cv_PAR%npar ,  d_XQR%Ap ,  cv_PAR%npar ,  0.d0 ,  d_XQR%XA_p ,  cv_PAR%p )
    !d_XQR%RHS(:cv_OBS%nobs,:)     = d_XQR%HQAp(:,:)
    !d_XQR%RHS((cv_OBS%nobs+1):,:) = d_XQR%XA_p(:,:)

    !The matrix multiplications can be avoided and the result is instead
    allocate(d_XQR%RHS((cv_OBS%nobs+cv_PAR%p),(cv_PAR%p+cv_A%K)))
    d_XQR%RHS = 0.d0
    d_XQR%RHS(:cv_OBS%nobs,:cv_A%K) = d_XQR%BC
    do i = 1, cv_PAR%p
        d_XQR%RHS(cv_OBS%nobs+i,cv_A%K+i) = 1.d0
    enddo

    end subroutine construct_RHS

    !**********************************************************
    ! Normalizing the left hand side and the right hand side
    !**********************************************************
    subroutine normalize_LHS_RHS(cv_PAR, cv_OBS, d_XQR)

    type (cv_param)     , intent(in)    :: cv_PAR
    type (cv_observ)    , intent(in)    :: cv_OBS
    type (kernel_XQR)   , intent(inout) :: d_XQR
    integer                             :: i

    !Compute one over the square-root of the left hand side
    allocate(d_XQR%C_S(cv_OBS%nobs+cv_PAR%p)) !This is the vector with the diagonal values of the scaling matrix
    d_XQR%C_S = UNINIT_REAL                   !Initialization
    do i=1,cv_OBS%nobs+cv_PAR%p
        if (d_XQR%LHS(i,i).eq.0.) then        !If the value on the LHS diagonal is 0 then the scaling factor will be 1
            d_XQR%C_S(i) = 1.d0
        else
            d_XQR%C_S(i) = 1/sqrt(abs(d_XQR%LHS(i,i)))!**(-0.5)
        endif
    enddo

    !Rescale the left- and right hand side by this diagonal matrix. This makes the diagonal of LHS equal to ones and zeros.
    do i=1,cv_OBS%nobs+cv_PAR%p
        d_XQR%LHS(i,:) = d_XQR%C_S(i) * d_XQR%LHS(i,:)
        d_XQR%LHS(:,i) = d_XQR%C_S(i) * d_XQR%LHS(:,i)
        d_XQR%RHS(i,:) = d_XQR%C_S(i) * d_XQR%RHS(i,:)
    enddo

    end subroutine normalize_LHS_RHS

    !*********************
    ! Rescale solution
    !*********************
    subroutine rescale_solution(cv_A, cv_PAR, cv_OBS, d_XQR)

    type (cv_algorithmic)   , intent(in)    :: cv_A
    type (cv_param)         , intent(in)    :: cv_PAR
    type (cv_observ)        , intent(in)    :: cv_OBS
    type (kernel_XQR)       , intent(inout) :: d_XQR
    integer                                 :: i

    !Allocate for Lambda*Ap (LAp) and M*Ap (MAp)
    allocate(d_XQR%LA_p(cv_OBS%nobs,(cv_A%K+cv_PAR%p)))
    allocate(d_XQR%MA_p(cv_PAR%p,(cv_A%K+cv_PAR%p)))
    
    !Using the scale matrix for rescaling the solution
    do i=1,cv_OBS%nobs
        d_XQR%LA_p(i,:) = d_XQR%C_S(i) * d_XQR%RHS(i,:)
    enddo
    do i=cv_OBS%nobs+1,cv_OBS%nobs+cv_PAR%p
        d_XQR%MA_p(i-cv_OBS%nobs,:) = d_XQR%C_S(i) * d_XQR%RHS(i,:)
    enddo

    end subroutine rescale_solution

    end module PCGA_solve_linear_equations