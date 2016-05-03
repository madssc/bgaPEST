    module PCGA_objective_function

    use bayes_pest_control
    use utilities
    include 'mkl_blas.fi' !MC: For ddot: https://software.intel.com/en-us/forums/intel-math-kernel-library/topic/478276
    contains

    !**********************************************************************
    ! This function computes the objective function for the PCGA version
    !**********************************************************************
    subroutine PCGA_cal_ob_funcs(d_XQR, d_S, cv_PAR, cv_OBS, d_OBS,  cv_A, d_PAR, b_ind, p_ind)

    implicit none
    type(kernel_XQR),       intent(inout)   :: d_XQR
    type(d_struct),         intent(inout)   :: d_S
    type(cv_param),         intent(in)      :: cv_PAR
    type(cv_algorithmic),   intent(in)      :: cv_A
    type(d_param),          intent(inout)   :: d_PAR
    type(cv_observ),        intent(in)      :: cv_OBS
    type(d_observ),         intent(in)      :: d_OBS
    integer,                intent(in)      :: b_ind, p_ind
    double precision,       allocatable     :: U_Zs(:), CU_Zs(:), work(:), IPIV(:)
    integer                                 :: i, j, k, info_LU, work_size

    !*********************************************************************
    !Regularization objective function phi_R = 1/2 s^T*U_Z*C^-1*U_Z^T*s
    !*********************************************************************
    
    !Computing inverse of truncated C if C is new
    if (p_ind .eq. 0) then
        if (b_ind > 1) deallocate(d_XQR%invtC) !K can change from one b_ind to another
        allocate(d_XQR%invtC(cv_A%K,cv_A%K))
        allocate(IPIV(cv_A%K))

        d_XQR%invtC = d_XQR%tC

        !LU decomposition:
        !dgetrf (  m=cv_A%K ,  n=cv_A%K ,  a=d_XQR%invtC ,  lda=cv_A%K ,  ipiv=IPIV ,  info=info_LU )
        call dgetrf (  cv_A%K ,  cv_A%K ,  d_XQR%invtC ,  cv_A%K ,  IPIV ,  info_LU )

        if (info_LU /= 0) then
            stop 'Truncated C-matrix is numerically singular!'
        end if

        !Invert matrix
        !dgetri (  n=cv_A%K ,  a=d_XQR%invtC ,  lda=cv_A%K ,  ipiv=IPIV ,  work=bestwork ,  lwork=-1 ,  info=info_LU )
        allocate(work(1))
        call dgetri (  cv_A%K ,  d_XQR%invtC ,  cv_A%K ,  IPIV ,  work ,  -1 ,  info_LU )
        work_size = int(work(1))
        deallocate(work)
        allocate(work(work_size))
        call dgetri (  cv_A%K ,  d_XQR%invtC ,  cv_A%K ,  IPIV ,  work ,  work_size ,  info_LU )
        deallocate(work)
        deallocate(IPIV)
    endif

    !Calculate tU_Z^T * s
    allocate(U_Zs(cv_A%K))
    allocate(CU_Zs(cv_A%K))
    !dgemv ( trans='t' ,  m=cv_A%K ,  n=cv_PAR%npar ,  alpha=1.d0 ,  a=d_XQR%tU_Z ,  lda=cv_A%K ,  x=d_PAR%pars ,  incx=1 ,  beta=0.d0 ,  y=sU_Z ,  incy=1 )
    call dgemv ( 't' ,  cv_PAR%npar, cv_A%K ,  1.d0 ,  d_XQR%tU_Z ,  cv_PAR%npar ,  d_PAR%pars ,  1 ,  0.d0 ,  U_Zs ,  1 )

    !Calculate tC * tU_Z^T * s
    !dgemv ( trans='n' ,  m=cv_A%K ,  n=cv_A%K ,  alpha=1.d0 ,  a=d_XQR%tC ,  lda=cv_A%K ,  x=U_Zs ,  incx=1 ,  beta=0.d0 ,  y=CU_Zs ,  incy=1 )
    call dgemv ( 'n' ,  cv_A%K ,  cv_A%K ,  1.d0 ,  d_XQR%invtC ,  cv_A%K ,  U_Zs ,  1 ,  0.d0 ,  CU_Zs ,  1 )

    !Calculate (tU_Z^T * s)^T * (tC * tU_Z^T * s)
    !res =  ddot ( n=cv_A%K ,  x=sU_ZC ,  incx=1 ,  y=sU_Z ,  incy=1)
    d_PAR%phi_R = 0.5 * ddot( cv_A%K ,  U_Zs ,  1 ,  CU_Zs ,  1)

    deallocate(U_Zs)
    deallocate(CU_Zs)

    !*********************************************************************
    !Misfit objective function phi_M = 1/2* (y-h(s))t * R^-1 * (y-h(s))
    !*********************************************************************
    
    !We use just a loop because R0 is diagonal
    if (associated(d_PAR%phi_M_vec)) deallocate(d_PAR%phi_M_vec)
    allocate(d_PAR%phi_M_vec(cv_OBS%nobsgp))
    
    !Calculate phi_M for each observational group. This is done so the misfit in each data set can be displayed.
    d_PAR%phi_M_vec = 0.d0
    do i = 1, cv_OBS%nobs
        k = 0
        do j = 1,cv_OBS%nobsgp
            if (d_OBS%group(i) == cv_OBS%grp_name(j)) k = j
        enddo
        if (k == 0) write (*,*) "Mismatch between group name in observation_groups and GroupName in observation_data"
        d_PAR%phi_M_vec(k) = d_PAR%phi_M_vec(k) + (1./(d_S%sig * d_XQR%R0(i,i)))*((d_OBS%obs(i)  - d_OBS%h(i))**2)
    enddo
    d_PAR%phi_M_vec = 0.5*d_PAR%phi_M_vec

    !Add all misfits together
    d_PAR%phi_M = 0.d0
    do i = 1,cv_OBS%nobsgp
        d_PAR%phi_M = d_PAR%phi_M + d_PAR%phi_M_vec(i)
    enddo
    
    !Total objective function phi_T= phi_R +phi_M
    d_PAR%phi_T = d_PAR%phi_R + d_PAR%phi_M

    end subroutine PCGA_cal_ob_funcs
    end module PCGA_objective_function