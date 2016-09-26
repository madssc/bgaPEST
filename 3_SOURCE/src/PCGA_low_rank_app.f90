    module PCGA_low_rank_app

    use bayes_pest_control
    use utilities
    use error_message
    use PCGA_non_complete_grid_modul
    use mkl_dfti !For fast Fourier transformation
    use mkl_vsl  !For random number generator
    include 'mkl_blas.fi'

    contains

    !********************************************************************************************************************************************************
    ! This subroutine construc a low rank approximation of the full covariance matrix.
    ! The method uses fast fourier transformations and random samples for this approximation.
    ! The method is inspired by:
    ! Kitanidis, P. K., and J. Lee (2014), Principal Component Geostatistical Approach for large-dimensional inverse problems, Water Resour. Res.
    ! Nowak, W., Tenkleve, S., and Cirpka, O. A., Efficient Computation of Linearized Cross-Covariance and Auto-Covariance Matrices of Interdependent Quantities, Mathematical Geology, Vol. 35, No. 1, January 2003
    !********************************************************************************************************************************************************

    subroutine Low_rank_app(cv_A, cv_PAR, cv_S, d_ANI, d_PAR, d_S, d_XQR, Q0_All)

    implicit none
    type(cv_algorithmic),  intent(inout)    :: cv_A
    type(cv_param),        intent(inout)    :: cv_PAR
    type(cv_struct),       intent(in)       :: cv_S
    type(d_anisotropy),    intent(in)       :: d_ANI
    type (d_param),        intent(in)       :: d_PAR
    type (d_struct),       intent(in)       :: d_S
    type(kernel_XQR),      intent(inout)    :: d_XQR
    type(Q0_compr),        intent(inout), pointer   :: Q0_All(:)
    type(DFTI_DESCRIPTOR), pointer          :: handle_FFT, handle_IFFT
    type(VSL_STREAM_STATE)                  :: stream_Gaussian
    double precision,      allocatable      :: work_Z(:), Sigma_Z(:), work_C(:), Sigma_C(:), UTz(:,:)
    double precision                        :: sqrt_eig, K_old
    complex(kind=selected_real_kind(30))    :: result_C
    real(kind=selected_real_kind(30))       :: result_C_real1, result_C_real2!, result_C_imag1, result_C_imag2
    integer                                 :: i, iii, j, k, p, lc, rc, cc, tmp_npar, tmp_Nrow, tmp_Ncol, tmp_Nlay, K_s_tmp
    integer                                 :: work_Z_size, info_Z, work_C_size, info_C, status, count_K, count_K_s, count_parameters, start_v, end_v, Nrow, Ncol, Nlay, NAR, NAC, NAL
    integer                                 :: extend_factor, extend_factor_max
    logical, save                           :: first_entry = .true. !Is used for deallocation
    double precision, allocatable           :: zr(:,:,:), zr2(:,:,:), zi(:,:,:), zi2(:,:,:), Z_sample_selected(:,:), eigenvalues(:), SCC_copy(:,:,:)
    integer, parameter                      :: out_unit=20

    !************************************************
    ! Low rank approximation of covariance matrix
    !************************************************

    !Nullify
    if (first_entry) then
        nullify(d_XQR%tU_Z)
        nullify(d_XQR%tC)
        nullify(d_XQR%SCC)
        nullify(d_XQR%SCCi)
        nullify(d_XQR%Z_sample_L)
        nullify(d_XQR%Z_sample)
        nullify(d_XQR%Ap)
        do p = 1, cv_PAR%p
            nullify(Q0_All(p)%tU_Z)
            nullify(Q0_All(p)%tC)
        enddo
        first_entry = .false.
    endif

    !Open file for eigenvalue test
    if (cv_A%eigenvalue_test == 1) then
        open (unit=out_unit,file="eigenvalue_test.txt",action="write",status="replace")
    endif

    !Constructing low rank approximation for each beta association
    count_parameters = 1
    do p = 1, cv_PAR%p

        !Non-complete grid:
        !This function is used if the grid spacing is not constant, if there is holes in the grid or if the edge is not rectangular
        if (cv_A%non_complete_grid == 1) then
            !Checks
            if (cv_A%store_Q /= .true.) stop "Store Q should be true if non regular grid should be used"
            if (cv_A%Q_compression_flag /= 1) stop "Q compression flag should be 1 if non regular grid should be used"
            if (cv_PAR%ndim > 3 .or. cv_PAR%ndim < 1) stop "The number of dimensions are incorrect if non regular grid should be used"

            call PCGA_non_complete_grid(cv_A, cv_PAR, d_ANI, d_PAR, Q0_All, p)

            !Saving the original variables
            tmp_npar = Q0_All(p)%npar
            tmp_Nrow = Q0_All(p)%Nrow
            tmp_Ncol = Q0_All(p)%Ncol
            tmp_Nlay = Q0_All(p)%Nlay

            !Incerting the fictive values
            Q0_All(p)%Nrow = Q0_All(p)%Nfrow
            Q0_All(p)%Ncol = Q0_All(p)%Nfcol
            Q0_All(p)%Nlay = Q0_All(p)%Nflay
            Q0_All(p)%npar = Q0_All(p)%fpar
        endif

        !Number of columns, rows and layers in the p'th model
        Nrow = Q0_All(p)%Nrow
        Ncol = Q0_All(p)%Ncol
        Nlay = Q0_All(p)%Nlay

        NAR  = 2*Nrow - 2 !Number of rows after appending
        NAC  = 2*Ncol - 2 !Number of columns after appending
        NAL  = 2*Nlay - 2 !Number of layers after appending

        !There is no appending in a non existing dimension
        if (Nrow < 2) NAR = 1
        if (Ncol < 2) NAC = 1
        if (Nlay < 2) NAL = 1

        !Transfer K and K_s values
        Q0_All(p)%K   = cv_PAR%K(p)
        Q0_All(p)%K_s = cv_PAR%K_s(p)

        !It is possible to automaticly use the maximum k-values
        if (cv_A%seed_Gaussian <= 0) then
            Q0_All(p)%K_s = NAR*NAC*NAL
            Q0_All(p)%K   = Q0_All(p)%npar
        endif

        !Make Q0_All(p)%K_s even
        if (mod(Q0_All(p)%K_s,2) == 1) Q0_All(p)%K_s = Q0_All(p)%K_s + 1
        
        !Construct the symtric circulant matrix (SC) or symmetric block circulant matrices with circulant blocks matrix (SCC). Saved in d_XQR%SCC
        !Note SCC is saved in a three dimensional matrix instead of one dimensional. This is convinient for the FFT.
        call Construct_SCC(d_XQR, cv_S, cv_PAR, cv_A, d_PAR, d_S, Q0_All, p)

        nullify(handle_FFT)
        nullify(handle_IFFT)

        !If the correlation length is to long compared with the model domain the eigenvalues of the covariance matrix can become negative.
        !The solution was to extend the model domain until all eigenvalues are positive.
        !At the moment the model domain is extended maximum ten times.
        !Note: It was experienced that a rectangular (non square) model did not improve with this technique. Maybe it would work better if the smallest dimension is extended or the model domain is made square and then extended.
        extend_factor = 1
        extend_factor_max = 10 !just a choise
        do while (extend_factor < extend_factor_max)

            !The imaginary part for the Fourier transformation
            if (associated(d_XQR%SCCi)) deallocate(d_XQR%SCCi)
            allocate(d_XQR%SCCi(NAR,NAC,NAL))
            d_XQR%SCCi = 0.d0

            !Setup for the Fourier transformations
            call setup_FFT(handle_FFT, handle_IFFT, NAR, NAC, NAL)

            !Performing a Fourier transformation of the SCC matrix. Then SCC becomes the (approximate) eigenvalues of SCC.
            !The introduction of (:,1,1) is made since FFT and IFFT works in one dimensional arrays, however SCC and SCCi are three dimensional. The structure of the 3D matrix is set in DFTI_INPUT_STRIDES whereby the FFT only need the memory adress of the first element. This is given by (:,1,1)
            status = DftiComputeForward(handle_FFT, d_XQR%SCC(:,1,1), d_XQR%SCCi(:,1,1))
            !d_XQR%SCC = sqrt(dble(NAR*NAC*NAL))*d_XQR%SCC !Eigenvalues should be rescaled according to Nowak & Cirpka 2013 equation 8. They assume the FFT/IFFT algorithm to have a rescaling constant in front of the summation (A1-A2) however that is not the case for this FFT/IFFT. Here there are no scaling, wherefore this rescaling should not be used.

            !SCCi is zero after the transformation because SCC is symmetric before the transformation.
            deallocate(d_XQR%SCCi)

            !If the eigenvalues are positive, then they are fine.
            if (minval(d_XQR%SCC) > 0) then
                exit
            else !Else the model domain is increased in all directions.
                if (extend_factor == 1) then
                    allocate(SCC_copy(NAR,NAC,NAL))
                    SCC_copy = d_XQR%SCC
                endif

                extend_factor = extend_factor + 1

                !The following is only in use if extend_factor < extend_factor_max.
                if (extend_factor < extend_factor_max) then
                    write (6,122) minval(d_XQR%SCC), p, extend_factor
122                 format("Note: The eigenvalues of the approximation of the covariance matrix Qss was negative ( ",EN10.1," ) for beta association ",I2,". The model domain is fictively increased by a factor of ",I2," in the sampling process for the approximation. This will not have an influence on the forward modelling, however will increase the workload for the computations related to the approxmation of Qss.")

                    !Reset the values
                    Nrow = Q0_All(p)%Nrow
                    Ncol = Q0_All(p)%Ncol
                    Nlay = Q0_All(p)%Nlay

                    !Increase those dimensions which exist
                    if (Nrow > 1) Nrow = extend_factor * Nrow
                    if (Ncol > 1) Ncol = extend_factor * Ncol
                    if (Nlay > 1) Nlay = extend_factor * Nlay

                    !Set the new number of Rows/Columns/Layers after appending
                    NAR  = 2*Nrow - 2
                    NAC  = 2*Ncol - 2
                    NAL  = 2*Nlay - 2

                    if (Q0_All(p)%Nrow < 2) NAR = 1
                    if (Q0_All(p)%Ncol < 2) NAC = 1
                    if (Q0_All(p)%Nlay < 2) NAL = 1

                    !Construct a SCC matrix which has the new size.
                    call Construct_extended_SCC(d_XQR, cv_S, cv_PAR, cv_A, d_PAR, d_S, Q0_All, p, extend_factor, Nrow, Ncol, Nlay, NAR, NAC, NAL) !remember to remove the extra parts

                    status = DftiFreeDescriptor(handle_FFT)
                    status = DftiFreeDescriptor(handle_IFFT)
                endif
            endif
        end do

        !If positive eigenvalues are not achieved then the algorithm cannot proceed. extend_factor_max could be increased.
        if (extend_factor == extend_factor_max) then
            deallocate(d_XQR%SCC)
            write (*,123) p, (extend_factor_max-1)
123         format('Warning: The eigenvalues of the approximation of the covariance matrix Qss for beta association ',I2,' is still negative after the model domain has been increased by a factor of ',I2,' in all directions. If this happens using the exponential variogram, try to lower the range. The program cannot continue and will therefore exit.')
            cv_A%eigenvalue_test = 1 !Will make the program exit just as if we were testing for the eigenvalues.
            exit
        endif

        !Sample from a unit normal distribution into d_XQR%Z_sample_L
        if (cv_A%seed_Gaussian > 0) then
            call sample_unit_normal_distribution(Q0_All, cv_A, d_XQR, p, NAR, NAC, NAL)
        elseif (cv_A%seed_Gaussian <= 0) then !set d_XQR%Z_sample_L to the identity matrix. Used for testing.
            allocate(d_XQR%Z_sample_L(NAR*NAC*NAL, Q0_All(p)%K_s))
            d_XQR%Z_sample_L = 0.d0
            do iii = 1, NAR*NAC*NAL
                d_XQR%Z_sample_L(iii,iii) = 1.d0
            enddo
        endif

        !Allocate for fouriertransformation of Z_sample and Z_sample_L
        allocate(zr (NAR,NAC,NAL))
        allocate(zi (NAR,NAC,NAL))

        if (associated(d_XQR%Z_sample)) deallocate(d_XQR%Z_sample)
        allocate(d_XQR%Z_sample(Q0_All(p)%npar, Q0_All(p)%K_s))
        d_XQR%Z_sample = 0.d0

        !Z_sample_L is made of samples for a unit (identity covariance matrix) normal distribution. This do-loop will make the random numbers follow the covariance matrix Qss
        do i = 1, Q0_All(p)%K_s, 2 !Two at a time.
            zr = 0.d0
            zi = 0.d0

            !Reshape NAR*NAC*NAL random numbers into a 3D matrix. zr and zi will not be correlated.
            call Reshape_PCGA(d_XQR%Z_sample_L(:,i),   zr, NAR, NAC, NAL)
            call Reshape_PCGA(d_XQR%Z_sample_L(:,i+1), zi, NAR, NAC, NAL)

            !Originally a FFT was performed however since zr and zi are uncorrelated and the Fourier tranformation is a orthogonal transformation the vectors the transformation would just result in other uncorrelated Gaussian white random vector.
            !This means that the process is not required. Information is found at https://en.wikipedia.org/wiki/White_noise -> White noise vector -> "Therefore, any orthogonal (Fourier) transformation of the vector will result in a Gaussian white random vector".
            !status = DftiComputeForward(handle_FFT, zr(:,1,1), zi(:,1,1))

            !Multiply the square-root of the eigenvalues of SCC onto the vectors zr and zi.
            do lc = 1, NAL
                do cc = 1, NAC
                    do rc = 1, NAR
                        sqrt_eig = sqrt(d_XQR%SCC(rc,cc,lc))
                        zr(rc,cc,lc) = sqrt_eig * zr(rc,cc,lc)
                        zi(rc,cc,lc) = sqrt_eig * zi(rc,cc,lc)
                    enddo
                enddo
            enddo

            !IFFT
            !zr and zi are not correlated after the transformation due to above wiki link.
            status = DftiComputeBackward(handle_IFFT, zr(:,1,1), zi(:,1,1))

            !Rescaling of FFT since the scaling-factor is 1.0 for both FFT and IFFT
            zr = zr/sqrt(dble(NAR*NAC*NAL))
            zi = zi/sqrt(dble(NAR*NAC*NAL))

            !Using reverse zero padding (Nrow, Ncol, Nlay instead of NAR, NAC, NAL) and incerting the result in Z_sample.
            call IReshape_PCGA(d_XQR%Z_sample(:,i),   zr, Q0_All(p)%Nrow, Q0_All(p)%Ncol, Q0_All(p)%Nlay)
            call IReshape_PCGA(d_XQR%Z_sample(:,i+1), zi, Q0_All(p)%Nrow, Q0_All(p)%Ncol, Q0_All(p)%Nlay)

            !Rescaling due to the random process.
            if (cv_A%seed_Gaussian > 0) then
                d_XQR%Z_sample(:,i) = d_XQR%Z_sample(:,i) / sqrt(dble(Q0_All(p)%K_s))
            endif
        end do

        !Reincert the non-extended SCC matrix into d_XQR%SCC and the parameters. Setting the Fourier transformation settings.
        if (extend_factor > 1) then
            Nrow = Q0_All(p)%Nrow
            Ncol = Q0_All(p)%Ncol
            Nlay = Q0_All(p)%Nlay

            NAR  = 2*Nrow - 2
            NAC  = 2*Ncol - 2
            NAL  = 2*Nlay - 2

            if (Nrow < 2) NAR = 1
            if (Ncol < 2) NAC = 1
            if (Nlay < 2) NAL = 1

            !NAR, NAC and NAL has changed, wherefore the Fourier settings should be reset.
            status = DftiFreeDescriptor(handle_FFT)
            status = DftiFreeDescriptor(handle_IFFT)
            call setup_FFT(handle_FFT, handle_IFFT, NAR, NAC, NAL)

            !Reincert the original SCC
            if (associated(d_XQR%SCC)) deallocate(d_XQR%SCC)
            allocate(d_XQR%SCC(NAR,NAC,NAL))
            d_XQR%SCC = SCC_copy
            deallocate(SCC_copy)
        endif

        !Picking out those nodes which really exist.
        if (cv_A%non_complete_grid == 1) then
            Q0_All(p)%Nrow = tmp_Nrow
            Q0_All(p)%Ncol = tmp_Ncol
            Q0_All(p)%Nlay = tmp_Nlay
            Q0_All(p)%npar = tmp_npar

            allocate(Z_sample_selected(Q0_All(p)%npar, Q0_All(p)%K_s))

            !Using findex (fictive index) to pick out those nodes which really exist.
            do i = 1, Q0_All(p)%npar
                Z_sample_selected(i,:) = d_XQR%Z_sample(Q0_All(p)%findex(i),:)
            enddo

            deallocate(d_XQR%Z_sample)
            allocate(d_XQR%Z_sample(Q0_All(p)%npar, Q0_All(p)%K_s))
            d_XQR%Z_sample = Z_sample_selected
            deallocate(Z_sample_selected)
        endif

        !Detrend Z_sample with respect to X using P = I - U_X*U_X^T.
        !UTz = U_X^T * Z_sample
        !dgemm ( transa='t' ,  transb='n' ,  m=cv_PAR%p ,  n=Q0_All(p)%K_s ,  k=Q0_All(p)%npar ,  alpha=1.d0 ,  a=d_XQR%U_X ,  lda=Q0_All(p)%npar ,  b=d_XQR%Z_sample ,  ldb=Q0_All(p)%npar ,  beta=0.d0 ,  c=UTz ,  ldc=cv_PAR%p )
        allocate(UTz(cv_PAR%p,Q0_All(p)%K_s))
        UTz = 0.d0
        call dgemm ( 't' ,  'n' ,  cv_PAR%p ,  Q0_All(p)%K_s ,  Q0_All(p)%npar ,  1.d0 ,  d_XQR%U_X(count_parameters:(count_parameters+Q0_All(p)%npar-1),:) ,  Q0_All(p)%npar ,  d_XQR%Z_sample(:,:Q0_All(p)%K_s) ,  Q0_All(p)%npar ,  0.d0 ,  UTz ,  cv_PAR%p )

        !Z_sample = Z_sample - U * UTz
        !This assumes that different beta associations do not share rows. This could be nice to avoid however, there are difficulties when computing C. This could however be overcome by using Toeplitz multiplication instead of FFT.
        !dgemm ( transa='n' ,  transb='n' ,  m=Q0_All(p)%npar ,  n=Q0_All(p)%K_s ,  k=cv_PAR%p ,  alpha=1.d0 ,  a=d_XQR%U_X ,  lda=Q0_All(p)%npar ,  b=UTz ,  ldb=cv_PAR%p ,  beta=0.d0 ,  c=UUTz ,  ldc=Q0_All(p)%npar )
        call dgemm ( 'n' ,  'n' ,  Q0_All(p)%npar ,  Q0_All(p)%K_s ,  cv_PAR%p ,  -1.d0 ,  d_XQR%U_X(count_parameters:(count_parameters+Q0_All(p)%npar-1),:) ,  Q0_All(p)%npar ,  UTz ,  cv_PAR%p ,  1.d0 ,  d_XQR%Z_sample(:,:Q0_All(p)%K_s) ,  Q0_All(p)%npar )
        deallocate(UTz)

        !SVD of Z:
        !d_XQR%Z_sample becomes U_Z, a matrix with the eigenvectors with the Q0_All(p)%K_s highest eigenvalues
        allocate(Sigma_Z(Q0_All(p)%npar))
        allocate(work_Z(1))
        !dgesvd ( jobu='O' ,  jobvt='N' ,  m=Q0_All(p)%npar ,  n=Q0_All(p)%K_s ,  a=d_XQR%Z_sample ,  lda=Q0_All(p)%npar ,  s=Sigma_Z ,  u=unused 0.d0 ,  ldu=Q0_All(p)%npar ,  vt=unused 0.d0 ,  ldvt=1 ,  work=work_Z ,  lwork=-1 ,  info=info_Z )
        call dgesvd ( 'O' ,  'N' ,  Q0_All(p)%npar ,  Q0_All(p)%K_s ,  d_XQR%Z_sample,  Q0_All(p)%npar ,  Sigma_Z , 0 ,  Q0_All(p)%npar ,  0 ,  Q0_All(p)%K_s ,  work_Z ,  -1 ,  info_Z )
        work_Z_size = int(work_Z(1))
        deallocate(work_Z)
        allocate(work_Z(work_Z_size))
        call dgesvd ( 'O' ,  'N' ,  Q0_All(p)%npar ,  Q0_All(p)%K_s ,  d_XQR%Z_sample,  Q0_All(p)%npar ,  Sigma_Z , 0.d0 ,  Q0_All(p)%npar ,  0.d0 ,  Q0_All(p)%K_s ,  work_Z ,  work_Z_size ,  info_Z )

        !If some of the eigenvalues are very small the eigenvectors contains almost no information. Therefore are they excluded.
        if (Q0_All(p)%K_s > Q0_All(p)%npar) Q0_All(p)%K_s = Q0_All(p)%npar
        K_s_tmp = Q0_All(p)%K_s
        do i = Q0_All(p)%K_s,1,-1
            if ((Sigma_Z(i)/Sigma_Z(1))<1.d-08) K_s_tmp = i-1
        enddo
        Q0_All(p)%K_s = K_s_tmp
        if (Q0_All(p)%K > Q0_All(p)%K_s) Q0_All(p)%K = Q0_All(p)%K_s
        
        deallocate(work_Z)
        deallocate(Sigma_Z)

        !Saving a truncated version of the eigenvectors.
        if (associated(Q0_All(p)%tU_Z)) deallocate(Q0_All(p)%tU_Z)
        allocate(Q0_All(p)%tU_Z(Q0_All(p)%npar,Q0_All(p)%K))
        Q0_All(p)%tU_Z = d_XQR%Z_sample(:,:Q0_All(p)%K) !Truncate to K instead of K_s

        !Extending Z_sample to the number of fictive nodes in the non-complete grid
        if (cv_A%non_complete_grid == 1) then
            allocate(Z_sample_selected(Q0_All(p)%fpar,Q0_All(p)%K_s))
            Z_sample_selected = 0.d0
            do i = 1, Q0_All(p)%npar
                Z_sample_selected(Q0_All(p)%findex(i),:) = d_XQR%Z_sample(i,:Q0_All(p)%K_s)
            enddo
            deallocate(d_XQR%Z_sample)
            allocate(d_XQR%Z_sample(Q0_All(p)%fpar,Q0_All(p)%K_s))
            d_XQR%Z_sample = Z_sample_selected
            deallocate(Z_sample_selected)
        endif

        !Computing C = U_Z * Qss * U_Z by use of fast Fourier transformation
        allocate(Q0_All(p)%C(Q0_All(p)%K_s,Q0_All(p)%K_s))
        allocate(zr2(NAR,NAC,NAL))
        allocate(zi2(NAR,NAC,NAL))
        do i = 1, Q0_All(p)%K_s
            zr = 0.d0
            zi = 0.d0

            !Reshape Npar random numbers into a 3D matrix, with zeros around
            !Using the i'th column in Z_sample
            call Reshape_PCGA(d_XQR%Z_sample(:,i), zr, Nrow, Ncol, Nlay)

            !Fast Fourier transformtion of zr and zi and rescaling due to FFT.
            status = DftiComputeForward(handle_FFT, zr(:,1,1), zi(:,1,1))
            zr = zr/sqrt(dble(NAR*NAC*NAL))
            zi = zi/sqrt(dble(NAR*NAC*NAL))

            do j = i, Q0_All(p)%K_s
                if (j /= i) then
                    zr2 = 0.d0
                    zi2 = 0.d0

                    !Reshape Npar random numbers into a 3D matrix, with zeros around
                    !Using the j'th column in Z_sample
                    call Reshape_PCGA(d_XQR%Z_sample(:,j), zr2, Nrow, Ncol, Nlay)

                    !Fast Fourier transformtion of zr2 and zi2 and rescaling due to FFT.
                    status = DftiComputeForward(handle_FFT, zr2(:,1,1), zi2(:,1,1))
                    zr2 = zr2/sqrt(dble(NAR*NAC*NAL))
                    zi2 = zi2/sqrt(dble(NAR*NAC*NAL))
                end if

                !Originally the imaginary part of the product was also computed, however since it always was close to zero is it omitted.
                !One could choose to compute this by cdot, however the two part are almost the same numbers with different sign. Then when the numbers are added round off errors are accumulated. This error can be much larger (e.g. 10^-8) than the correct answer (e.g. 10^-16).
                !Therefore the numbers are computed separatly with high precision, whereafter they are added. In this way the error is as small as the correct answer (10^-16).
                result_C_real1  = 0.0
                result_C_real2  = 0.0
                !result_C_imag1 = 0.0
                !result_C_imag2 = 0.0

                !FFT(u_1)^H * Lambda * FFT(u_2)
                if (j /= i) then
                    do lc = 1, NAL
                        do cc = 1, NAC
                            do rc = 1, NAR
                                result_C_real1  = result_C_real1 + zr(rc,cc,lc) * d_XQR%SCC(rc,cc,lc) * zr2(rc,cc,lc)
                                result_C_real2  = result_C_real2 + zi(rc,cc,lc) * d_XQR%SCC(rc,cc,lc) * zi2(rc,cc,lc)
                                !result_C_imag1 = result_C_imag1 + zr(rc,cc,lc) * d_XQR%SCC(rc,cc,lc) * zi2(rc,cc,lc)
                                !result_C_imag2 = result_C_imag2 + zi(rc,cc,lc) * d_XQR%SCC(rc,cc,lc) * zr2(rc,cc,lc)
                            enddo
                        enddo
                    enddo
                elseif (j == i) then !FFT(u_1)^H * Lambda * FFT(u_1)
                    do lc = 1, NAL
                        do cc = 1, NAC
                            do rc = 1, NAR
                                result_C_real1 = result_C_real1 + d_XQR%SCC(rc,cc,lc) * (zr(rc,cc,lc)**2 + zi(rc,cc,lc)**2)
                            enddo
                        enddo
                    enddo
                endif

                result_C_real1 = result_C_real1 + result_C_real2
                !result_C_imag1 = result_C_imag1 - result_C_imag2
                !result_C = cmplx(result_C_real1, result_C_imag1)

                !The C matrix is symmetric and therefore the value is saved twice. The real function is probably not needed.
                Q0_All(p)%C(i,j) = real(result_C_real1)
                Q0_All(p)%C(j,i) = real(result_C_real1)
            end do
        end do

        !If the eigenvalues of the covariance matrix Qss is set to be printed these values are found by computing the eigenvalues of the C-matrix.
        if (cv_A%eigenvalue_test == 1) then
            !Computing the SVD of C
            allocate(Sigma_C(Q0_All(p)%K_s))
            allocate(work_C(1))
            !dgesvd ( jobu='O' ,  jobvt='N' ,  m=Q0_All(p)%npar ,  n=Q0_All(p)%K_s ,  a=d_XQR%Z_sample ,  lda=Q0_All(p)%npar ,  s=Sigma_Z ,  u=unused 0.d0 ,  ldu=Q0_All(p)%npar ,  vt=unused 0.d0 ,  ldvt=1 ,  work=work_Z ,  lwork=-1 ,  info=info_Z )
            call dgesvd ( 'N' ,  'N' ,  Q0_All(p)%K_s ,  Q0_All(p)%K_s ,  Q0_All(p)%C,  Q0_All(p)%K_s ,  Sigma_C , 0 ,  Q0_All(p)%K_s ,  0 ,  Q0_All(p)%K_s ,  work_C ,  -1 ,  info_C )
            work_C_size = int(work_C(1))
            deallocate(work_C)
            allocate(work_C(work_C_size))
            call dgesvd ( 'N' ,  'N' ,  Q0_All(p)%K_s ,  Q0_All(p)%K_s ,  Q0_All(p)%C,  Q0_All(p)%K_s ,  Sigma_C , 0 ,  Q0_All(p)%K_s ,  0 ,  Q0_All(p)%K_s ,  work_C ,  work_C_size ,  info_C )
            deallocate(work_C)

            !Normalize the singular values with respect to the highest value. Following this the values are printed to out_unit.
            Sigma_C = Sigma_C/maxval(Sigma_C)
            write (out_unit,'(*(D18.10))') (Sigma_C)
            deallocate(Sigma_C)
        endif

        !Truncate C
        if (associated(Q0_All(p)%tC))   deallocate(Q0_All(p)%tC)
        allocate(Q0_All(p)%tC(Q0_All(p)%K,Q0_All(p)%K))
        Q0_All(p)%tC   = Q0_All(p)%C(:Q0_All(p)%K,:Q0_All(p)%K)

        !Setting the counter of number of passed parameters
        count_parameters = count_parameters + Q0_All(p)%npar

        status = DftiFreeDescriptor(handle_FFT)
        status = DftiFreeDescriptor(handle_IFFT)

        deallocate(zr)
        deallocate(zr2)
        deallocate(zi)
        deallocate(zi2)
        deallocate(Q0_All(p)%C)
        deallocate(d_XQR%SCC)
        deallocate(d_XQR%Z_sample)
        deallocate(d_XQR%Z_sample_L)
    end do

    if (cv_A%eigenvalue_test == 1) then
        close (out_unit) !Closing the file out_unit
    else
        !Counting the total number of K
        cv_A%K = 0
        do p = 1, cv_PAR%p
            cv_A%K = cv_A%K + Q0_All(p)%K
        enddo

        !All tU_Z and tC is incerted into d_XQR
        if (associated(d_XQR%tU_Z)) deallocate(d_XQR%tU_Z)
        if (associated(d_XQR%tC))   deallocate(d_XQR%tC)
        allocate(d_XQR%tU_Z(cv_PAR%npar,cv_A%K))
        allocate(d_XQR%tC(cv_A%K, cv_A%K))
        count_K = 1
        count_parameters = 1
        d_XQR%tU_Z = 0.d0
        d_XQR%tC = 0.d0
        do p = 1,cv_PAR%p
            d_XQR%tU_Z(count_parameters:count_parameters+Q0_All(p)%npar-1, count_K:count_K+Q0_All(p)%K-1) = Q0_All(p)%tU_Z
            d_XQR%tC(count_K:count_K+Q0_All(p)%K-1,count_K:count_K+Q0_All(p)%K-1)                 = Q0_All(p)%tC
            count_K      = count_K      + Q0_All(p)%K
            count_parameters = count_parameters + Q0_All(p)%npar
            if (associated(Q0_All(p)%tU_Z)) deallocate(Q0_All(p)%tU_Z)
            if (associated(Q0_All(p)%tC))   deallocate(Q0_All(p)%tC)
        enddo

        !Declairing Ap = [tU_Z, U_X]
        if (associated(d_XQR%Ap)) deallocate(d_XQR%Ap)
        allocate(d_XQR%Ap(cv_PAR%npar,cv_A%K+cv_PAR%p))
        d_XQR%Ap(:,:cv_A%K) = d_XQR%tU_Z(:,:)
        d_XQR%Ap(:,(cv_A%K+1):) = d_XQR%U_X
    endif

    end subroutine Low_rank_app

    !*****************************************************************************
    ! Construct the Symmetric Circulant matrix of Symmetric Circulant matrices
    !*****************************************************************************
    subroutine Construct_SCC(d_XQR, cv_S, cv_PAR, cv_A, d_PAR, d_S, Q0_All, p)

    type(kernel_XQR),       intent(inout)           :: d_XQR
    type(cv_struct),        intent(in)              :: cv_S
    type(cv_param),         intent(in)              :: cv_PAR
    type(cv_algorithmic),   intent(in)              :: cv_A
    type (d_param),         intent(in)              :: d_PAR
    type(d_struct),         intent(in)              :: d_S
    type(Q0_compr),         intent(in), pointer     :: Q0_All(:)
    integer,                intent(in)              :: p
    double precision,       allocatable             :: A(:)
    integer                                         :: Nrow, Ncol, Nlay, NAR, NAC, NAL, lc, cc, place, i
    double precision                                :: theta_1, theta_2, Lmax

    !The variance of the covariance matrix
    theta_1 = d_S%theta(Q0_All(p)%BetaAss,1)

    Nrow = Q0_All(p)%Nrow
    Ncol = Q0_All(p)%Ncol
    Nlay = Q0_All(p)%Nlay

    NAR  = 2*Nrow - 2; !Number of rows after appending
    NAC  = 2*Ncol - 2; !Number of columns after appending
    NAL  = 2*Nlay - 2; !Number of layers after appending

    if (Nrow < 2) NAR = 1
    if (Ncol < 2) NAC = 1
    if (Nlay < 2) NAL = 1

    !Calculate the first column of the covariance matrix
    allocate(A(Q0_All(p)%npar))
    select case (cv_S%var_type(Q0_All(p)%BetaAss)) !Variogram type
    case (0) !Nugget:
        A(:) = 0.d0
        do i = 1,(Nrow*Ncol*Nlay),Nrow !Changed from Ncol after 1. pub
            A(i) = theta_1
        enddo

    case (1) !Linear:
        Lmax = d_XQR%L
        A(:) = theta_1 * Lmax * exp(-Q0_All(p)%Q0_C(:,1)/Lmax)

    case (2) !Exponential:
        theta_2 = d_S%theta(Q0_All(p)%BetaAss,2)
        A(:) = theta_1 * exp(-Q0_All(p)%Q0_C(:,1)/theta_2)

    case (3) !Anisotropic exponential correlation function:
        if (cv_A%par_anisotropy == 1) write(*,*) "anisotropy is not implemented in anisotropic correlation function: var_type 3)"
        Q0_All(p)%Q0_C(j,1) = 0.d0
        do j=2, Q0_All(p)%npar
            do k = 1,cv_PAR%ndim ! calculate the distances
                Q0_All(p)%Q0_C(j,1) = Q0_All(p)%Q0_C(j,1) - abs(d_PAR%lox(Q0_All(p)%Beta_Start,k) - d_PAR%lox(Q0_All(p)%Beta_Start+j-1,k))/d_S%theta(Q0_All(p)%BetaAss,k+1)
            enddo
        enddo
        A(:) = theta_1 * exp(Q0_All(p)%Q0_C(j,1))

        case default
        write (*,*) "ERROR: Variogram type is incorrect"
        stop
    end select

    !Allocating SCC and reshaping the first column of the covariance matrix into it
    if (associated(d_XQR%SCC)) deallocate(d_XQR%SCC)
    allocate(d_XQR%SCC(NAR,NAC,NAL))
    d_XQR%SCC = 0.d0
    call Reshape_PCGA(A, d_XQR%SCC, Nrow, Ncol, Nlay)
    deallocate(A)

    !Make the circulance for rows
    if (Nrow > 2) then
        do lc = 1,Nlay
            do cc = 1,Ncol
                place = Nrow+1
                do i = Nrow-1,2,-1
                    d_XQR%SCC(place,cc,lc) = d_XQR%SCC(i,cc,lc)
                    place = place + 1
                enddo
            enddo
        enddo
    endif

    !Make the circulance for columns
    if (Ncol > 2) then
        do lc = 1,Nlay
            place = Ncol+1
            do i = Ncol-1,2,-1
                d_XQR%SCC(:,place,lc) = d_XQR%SCC(:,i,lc)
                place = place + 1
            enddo
        enddo
    endif

    !Make the circulance for layers
    if (Nlay > 2) then
        place = Nlay+1
        do i = Nlay-1,2,-1
            d_XQR%SCC(:,:,place) = d_XQR%SCC(:,:,i)
            place = place + 1
        enddo
    endif

    end subroutine Construct_SCC

    !******************************************************************************************************************
    !*** An extended version of the SCC matrix when the eigenvalues estimated through FFT computations are negative
    !******************************************************************************************************************
    subroutine Construct_extended_SCC(d_XQR, cv_S, cv_PAR, cv_A, d_PAR, d_S, Q0_All, p, extend_factor, Nrow_mul, Ncol_mul, Nlay_mul, NAR_mul, NAC_mul, NAL_mul)

    type(kernel_XQR),       intent(inout)           :: d_XQR
    type(cv_struct),        intent(in)              :: cv_S
    type(cv_param),         intent(in)              :: cv_PAR
    type(cv_algorithmic),   intent(in)              :: cv_A
    type (d_param),         intent(in)              :: d_PAR
    type(d_struct),         intent(in)              :: d_S
    type(Q0_compr),         intent(in), pointer     :: Q0_All(:)
    integer,                intent(in)              :: p, extend_factor !>1
    double precision,       allocatable             :: A(:)
    integer,                intent(in)              :: Nrow_mul, Ncol_mul, Nlay_mul, NAR_mul, NAC_mul, NAL_mul
    integer                                         :: Nrow, Ncol, Nlay, lc, cc, place, i, ii, j, jj, k, kk
    double precision                                :: theta_1, theta_2, Lmax
    double precision                                :: row_min, col_min, lay_min, row_max, col_max, lay_max, p_row, p_col, p_lay
    integer, parameter                              :: out_unit2=22

    size_mul_row = extend_factor
    size_mul_col = extend_factor
    size_mul_lay = extend_factor

    Nrow = Q0_All(p)%Nrow
    Ncol = Q0_All(p)%Ncol
    Nlay = Q0_All(p)%Nlay

    if (Nrow < 2) size_mul_row = 1
    if (Ncol < 2) size_mul_col = 1
    if (Nlay < 2) size_mul_lay = 1

    if (associated(d_XQR%SCC)) deallocate(d_XQR%SCC)
    allocate(d_XQR%SCC(NAR_mul,NAC_mul,NAL_mul))

    !Finding the coordinates of the first parameter and the distance to the last parameter. It is assumed that they are the two points furthest away from each other.
    if (Nrow > 1) then
        row_min     = d_PAR%lox(Q0_All(p)%Beta_Start,1)
        row_max     = d_PAR%lox(Q0_All(p)%Beta_Start+Q0_All(p)%npar-1,1) + d_PAR%lox(Q0_All(p)%Beta_Start+1,1) - 2*row_min
    endif
    if (Ncol > 1) then
        col_min     = d_PAR%lox(Q0_All(p)%Beta_Start,2)
        col_max     = d_PAR%lox(Q0_All(p)%Beta_Start+Q0_All(p)%npar-1,2) + d_PAR%lox(Q0_All(p)%Beta_Start+Nrow+1,2) - 2*col_min
    endif
    if (Nlay > 1) then
        lay_min     = d_PAR%lox(Q0_All(p)%Beta_Start,3)
        lay_max     = d_PAR%lox(Q0_All(p)%Beta_Start+Q0_All(p)%npar-1,3) + d_PAR%lox(Q0_All(p)%Beta_Start+Nrow*Ncol+1,3) - 2*lay_min
    endif

    !Computing the distances between all points in the extended model domain
    if (cv_S%var_type(Q0_All(p)%BetaAss) /= 0 .or. cv_S%var_type(Q0_All(p)%BetaAss) /= 3) then
        d_XQR%SCC = 0.d0
        do kk = 0,size_mul_lay-1
            do k  = 1,Nlay
                do jj = 0,size_mul_col-1
                    do j  = 1,Ncol
                        do ii = 0,size_mul_row-1
                            do i  = 1,Nrow
                                place = i + (j-1)*Nrow + (k-1)*Ncol*Nrow
                                if (Nrow > 1) then
                                    p_row = d_PAR%lox(Q0_All(p)%Beta_Start+place-1,1)
                                    d_XQR%SCC(i + ii*Nrow, j + jj*Ncol, k + kk*Nlay) = d_XQR%SCC(i + ii*Nrow, j + jj*Ncol, k + kk*Nlay) + (p_row-row_min + dble(ii)*row_max)**2
                                endif
                                if (Ncol > 1) then
                                    p_col = d_PAR%lox(Q0_All(p)%Beta_Start+place-1,2)
                                    d_XQR%SCC(i + ii*Nrow, j + jj*Ncol, k + kk*Nlay) = d_XQR%SCC(i + ii*Nrow, j + jj*Ncol, k + kk*Nlay) + (p_col-col_min + dble(jj)*col_max)**2
                                endif
                                if (Nlay > 1) then
                                    p_lay = d_PAR%lox(Q0_All(p)%Beta_Start+place-1,3)
                                    d_XQR%SCC(i + ii*Nrow, j + jj*Ncol, k + kk*Nlay) = d_XQR%SCC(i + ii*Nrow, j + jj*Ncol, k + kk*Nlay) + (p_lay-lay_min + dble(kk)*lay_max)**2
                                endif
                                d_XQR%SCC(i + ii*Nrow, j + jj*Ncol, k + kk*Nlay) = sqrt(d_XQR%SCC(i + ii*Nrow, j + jj*Ncol, k + kk*Nlay))
                            enddo
                        enddo
                    enddo
                enddo
            enddo
        enddo
    endif

    !Convert the distances in SCC to the covariance matrix (first column).
    theta_1 = d_S%theta(Q0_All(p)%BetaAss,1)
    select case (cv_S%var_type(Q0_All(p)%BetaAss)) !Variogram type
    case (0) !Nugget:
        d_XQR%SCC(:Nrow_mul,:Ncol_mul,:Nlay_mul) = theta_1

    case (1) !Linear:
        Lmax = d_XQR%L
        d_XQR%SCC(:Nrow_mul,:Ncol_mul,:Nlay_mul) = theta_1 * Lmax * exp(-d_XQR%SCC(:Nrow_mul,:Ncol_mul,:Nlay_mul)/Lmax)

    case (2) !Exponential:
        theta_2 = d_S%theta(Q0_All(p)%BetaAss,2)
        d_XQR%SCC(:Nrow_mul,:Ncol_mul,:Nlay_mul) = theta_1 * exp(-d_XQR%SCC(:Nrow_mul,:Ncol_mul,:Nlay_mul)/theta_2)

    case (3) !Anisotropy exponential correlation function
        write(*,*) "Warning: Error, the anisotropy correlation function should never end up here.. It should always have positive eigenvalues according to Dietrick and Newsam 1997 p. 1098."

        case default
        write (*,*) "ERROR: Variogram type is incorrect"
        stop
    end select

    !Make the circulance for rows
    if (Nrow_mul > 2) then
        do lc = 1,Nlay_mul
            do cc = 1,Ncol_mul
                place = Nrow_mul+1
                do i = Nrow_mul-1,2,-1
                    d_XQR%SCC(place,cc,lc) = d_XQR%SCC(i,cc,lc)
                    place = place + 1
                enddo
            enddo
        enddo
    endif

    !Make the circulance for columns
    if (Ncol_mul > 2) then
        do lc = 1,Nlay_mul
            place = Ncol_mul+1
            do i = Ncol_mul-1,2,-1
                d_XQR%SCC(:,place,lc) = d_XQR%SCC(:,i,lc)
                place = place + 1
            enddo
        enddo
    endif

    !Make the circulance for layers
    if (Nlay_mul > 2) then
        place = Nlay_mul+1
        do i = Nlay_mul-1,2,-1
            d_XQR%SCC(:,:,place) = d_XQR%SCC(:,:,i)
            place = place + 1
        enddo
    endif

    end subroutine Construct_extended_SCC

    !**************************************
    ! Reshape: vector to matrix
    !**************************************
    subroutine Reshape_PCGA(vectorin, matrixout, NR, NC, NL)

    double precision,   intent(in)      :: vectorin(:) !in
    double precision,   intent(inout)   :: matrixout(:,:,:)
    integer, 			intent(in) 	    :: NR, NC, NL
    integer :: i, j, place

    !Reshaping the elements of the vector into the matrix
    place = 1
    do j = 1,NL
        do i = 1,NC
            matrixout(:NR,i,j) = vectorin(place:place+NR-1)
            place = place + NR
        enddo
    enddo

    end subroutine Reshape_PCGA

    !**************************************
    ! Inverse reshape: matrix to vector
    !**************************************
    subroutine IReshape_PCGA(vectorout, matrixin, NR, NC, NL)

    double precision,   intent(inout)   :: vectorout(:)
    double precision,   intent(in)      :: matrixin(:,:,:)
    integer, 			intent(in) 	    :: NR, NC, NL
    integer :: i, j, place

    !Reshaping the elements of the matrix into the vector
    place = 1
    do j = 1,NL
        do i = 1,NC
            vectorout(place:place+NR-1) = matrixin(:NR,i,j)
            place = place + NR
        enddo
    enddo

    end subroutine IReshape_PCGA

    !**********************************************************************************
    ! Setup for the Fourier transformations and the inverse Fourier transformations
    !**********************************************************************************
    subroutine setup_FFT(handle_FFT, handle_IFFT, NAR, NAC, NAL)

    type(DFTI_DESCRIPTOR), intent(inout), pointer :: handle_FFT, handle_IFFT
    integer, intent(in)   :: NAR, NAC, NAL
    integer, dimension(3) :: L

    !Declair length of individual dimensions of the matrix which has to be fourier transformed
    L = (/ NAR, NAC, NAL/)

    !FFT setup
    status = DftiCreateDescriptor(handle_FFT, DFTI_DOUBLE, DFTI_COMPLEX, 3, L)          !Doule precission, complex values are used as in and output, 3 dimensions ,<dfti_length> = L describe the size of the 3 dimensions.
    status = DftiSetValue(handle_FFT, DFTI_COMPLEX_STORAGE,  DFTI_REAL_REAL)            !This can become less memory demanding by using DFTI_CONJUGATE_EVEN_STORAGE and DFTI_PACKED_FORMAT
    status = DftiSetValue(handle_FFT, DFTI_PLACEMENT, DFTI_INPLACE)                     !In and output is saved in the same variable
    !These are just the default settings:
    !status = DftisetValue(handle_FFT,  DFTI_INPUT_STRIDES, (/0, 1, NAR, NAR*NAC/))     !Zero offset s0 = 0, data is continues in the first dimension s1 = 1, for the second dimension there are Nrow points between every data point in the matrix s2 = Nrow and Nrow*Ncol in the third dimension s3 = Nrow*Ncol.
    !status = DftisetValue(handle_FFT,  DFTI_OUTPUT_STRIDES, (/0, 1, Nrow, Nrow*Ncol/)  !For in-place transforms ( DFTI_PLACEMENT=DFTI_INPLACE ), the configuration set by  DFTI_OUTPUT_STRIDES is ignored when the element types in the forward and backward domains are the same.
    status = DftiCommitDescriptor(handle_FFT)                                           !Make the descriptor ready for the transformation

    !IFFT setup
    status = DftiCreateDescriptor(handle_IFFT, DFTI_DOUBLE, DFTI_COMPLEX, 3, L)         !<dfti_length> = L
    status = DftiSetValue(handle_IFFT, DFTI_COMPLEX_STORAGE,  DFTI_REAL_REAL)
    status = DftiSetValue(handle_IFFT, DFTI_PLACEMENT, DFTI_INPLACE)
    !These are just the default settings:
    !status = DftisetValue(handle_IFFT,  DFTI_INPUT_STRIDES, (/0, 1, NAR, NAR*NAC/))    !Sero offset s0 = 0, data is continues in the first dimension s1 = 1, for the second dimension there are Nrow points between every data point in the matrix s2 = Nrow and Nrow*Ncol in the third dimension s3 = Ncol*Nrow.
    !status = DftisetValue(handle_IFFT,  DFTI_OUTPUT_STRIDES, (/0, 1, Nrow, Nrow*Ncol/) !For in-place transforms ( DFTI_PLACEMENT=DFTI_INPLACE ), the configuration set by  DFTI_OUTPUT_STRIDES is ignored when the element types in the forward and backward domains are the same.
    status = DftiCommitDescriptor(handle_IFFT)                                          !Make the descriptor ready for the inverse transformation

    !Note: The scaling factor for both FFT and IFFT is 1.0 by default
    end subroutine setup_FFT

    !**********************************************************************************************************************
    ! This function sample random values from a gaussian distribution with mean = 0 and variance = 1 and covariance = 0
    !**********************************************************************************************************************
    subroutine sample_unit_normal_distribution(Q0_All, cv_A, d_XQR, p, NAR, NAC, NAL)

    type(cv_algorithmic),  intent(inout)    :: cv_A
    type(Q0_compr),   intent(in),   pointer :: Q0_All(:)
    type(kernel_XQR), intent(inout)         :: d_XQR
    integer,          intent(in)            :: p, NAR, NAC, NAL
    double precision, allocatable           :: mu(:), diagonal_covariance(:)
    integer                                 :: status

    !Mean value
    allocate(mu(NAR*NAC*NAL))
    mu = 0.d0

    !Variance. Since diagonal_covariance only is a vector it is implicit assumed that the covariance is zero.
    allocate(diagonal_covariance(NAR*NAC*NAL))
    diagonal_covariance = 1.d0

    !Matrix for saving random values
    if (associated(d_XQR%Z_sample_L)) deallocate(d_XQR%Z_sample_L)
    allocate(d_XQR%Z_sample_L(NAR*NAC*NAL, Q0_All(p)%K_s))
    d_XQR%Z_sample_L = 0.d0

    !Sample random values
    !status =  vdrnggaussianmv (  method=method_Gaussian ,  stream=stream_Gaussian ,  n=Q0_All(p)%K_s ,  r=d_XQR%Z_sample_L ,  dimen=NAR*NAC*NAL ,  mstorage=VSL_MATRIX_STORAGE_DIAGONAL ,  a=d_XQR%mu ,  t=Diagonal_covariance )
    status = vdrnggaussianmv (  cv_A%method_Gaussian ,  cv_A%stream_Gaussian ,  Q0_All(p)%K_s ,  d_XQR%Z_sample_L ,  NAR*NAC*NAL ,  VSL_MATRIX_STORAGE_DIAGONAL ,  mu ,  Diagonal_covariance )

    deallocate(mu)
    deallocate(diagonal_covariance)

    end subroutine sample_unit_normal_distribution

    !************************************
    ! Seeding for the random sampling
    !************************************
    subroutine seed_random_samples(cv_A)

    type(cv_algorithmic),  intent(inout)    :: cv_A
    integer(KIND=4)                         :: clock_system

    !Seed value for the random sampling
    !If seed_Gaussian = 123 it is not changed i.e. samples are the same every time the program is executed
    if (abs(cv_A%seed_Gaussian) /= 123) then !not equal
        !call cpu_time(clock) !Alternative method
        !call itime(time_val) !Alternative method
        call system_clock(count=clock_system)
        cv_A%seed_Gaussian = clock_system
    end if

    !Initialize stream
    if     (cv_A%brng_Gaussian == 1) then
        cv_A%brng_Gaussian=VSL_BRNG_MCG31
    elseif (cv_A%brng_Gaussian == 2) then
        cv_A%brng_Gaussian=VSL_BRNG_R250
    elseif (cv_A%brng_Gaussian == 3) then
        cv_A%brng_Gaussian=VSL_BRNG_MRG32K3A
    elseif (cv_A%brng_Gaussian == 4) then
        cv_A%brng_Gaussian=VSL_BRNG_MCG59
    elseif (cv_A%brng_Gaussian == 5) then
        cv_A%brng_Gaussian=VSL_BRNG_WH
    elseif (cv_A%brng_Gaussian == 6) then
        cv_A%brng_Gaussian=VSL_BRNG_MT19937
    elseif (cv_A%brng_Gaussian == 7) then
        cv_A%brng_Gaussian=VSL_BRNG_MT2203
    elseif (cv_A%brng_Gaussian == 8) then
        cv_A%brng_Gaussian=VSL_BRNG_SFMT19937
    elseif (cv_A%brng_Gaussian == 9) then
        cv_A%brng_Gaussian=VSL_BRNG_SOBOL
    elseif (cv_A%brng_Gaussian == 10) then
        cv_A%brng_Gaussian=VSL_BRNG_NIEDERR
        !elseif (cv_A%brng_Gaussian == 11) then !Of unknown reasons does this not work
        !    cv_A%brng_Gaussian=VSL_BRNG_IABSTRACT
        !elseif (cv_A%brng_Gaussian == 12) then !Of unknown reasons does this not work
        !    cv_A%brng_Gaussian=VSL_BRNG_DABSTRACT
        !elseif (cv_A%brng_Gaussian == 13) then !Of unknown reasons does this not work
        !    cv_A%brng_Gaussian=VSL_BRNG_SABSTRACT
    elseif (cv_A%brng_Gaussian == 14) then
        cv_A%brng_Gaussian=VSL_BRNG_NONDETERM
    elseif (cv_A%brng_Gaussian == 15) then
        cv_A%brng_Gaussian=VSL_BRNG_PHILOX4X32X10
    elseif (cv_A%brng_Gaussian == 16) then
        cv_A%brng_Gaussian=VSL_BRNG_ARS5
    else
        cv_A%brng_Gaussian=VSL_BRNG_MCG31
    endif

    !Gaussian sample generation methods
    if     (cv_A%method_Gaussian  <= 1) then
        cv_A%method_Gaussian=VSL_RNG_METHOD_GAUSSIANMV_BOXMULLER
    elseif (cv_A%method_Gaussian  == 2) then
        cv_A%method_Gaussian=VSL_RNG_METHOD_GAUSSIANMV_BOXMULLER2
    elseif (cv_A%method_Gaussian  == 3) then
        vmethod_Gaussian=VSL_RNG_METHOD_GAUSSIANMV_ICDF
    elseif (cv_A%method_Gaussian > 3) then
        cv_A%method_Gaussian=VSL_RNG_METHOD_GAUSSIANMV_BOXMULLER2
    endif

    !Stream initialization
    !vslnewstream(  stream ,  brng ,  seed )
    status = vslnewstream(cv_A%stream_Gaussian, cv_A%brng_Gaussian, cv_A%seed_Gaussian)

    end subroutine seed_random_samples

    end module PCGA_low_rank_app
