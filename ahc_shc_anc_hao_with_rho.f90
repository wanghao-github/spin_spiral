program ahc_shc_anc_hao
    use mpi
    use pauli  ! 使用 pauli 模块，其中定义了 dp
    implicit none

    ! 现在，dp 是从 pauli 模块中获取的，无需在主程序中重新定义
    ! integer, parameter :: dp = kind(1.0d0)  ! 已移除

    ! 声明变量
    complex(kind=dp), allocatable :: hops(:,:,:)
    complex(kind=dp), allocatable :: rsnabla(:,:,:,:)
    complex(kind=dp), allocatable :: rspauli(:,:,:,:)
    complex(kind=dp), allocatable :: pauli_matrix(:,:,:)  ! 避免与模块同名
    complex(kind=dp), allocatable :: paulifft(:,:,:) 
    complex(kind=dp), allocatable :: paulifft2(:,:,:) 
    complex(kind=dp), allocatable :: rho(:,:), rho_mpi(:,:)
    real(kind=dp) :: rdum, idum 
    integer :: ix, iy, iz, band1, band2, h, num_wann, m, n
    integer :: ik1, ik2, ik3 
    real(kind=dp) :: twopi, temperature
    real(kind=dp) :: phas 
    complex(kind=dp) :: fac, fac2 
    complex(kind=dp), allocatable :: ham(:,:) 
    real(kind=dp) :: vl, vu 
    integer :: ne, j 
    real(kind=dp) :: abstol, time_start, time_end, time_start1, time_end1
    real, allocatable :: eigvals(:) 
    complex(kind=dp), allocatable :: eigvecs(:,:), eigvecs_dag(:,:), mat_temp(:,:)
    integer :: info 
    complex(kind=dp), allocatable :: work(:) 
    integer :: lwork 
    integer, allocatable :: iwork(:) 
    real, allocatable :: rwork(:) 
    integer, allocatable :: ifail(:) 
    real :: kpoints(3)
    real :: scale, Beta_fake, mu
    integer :: maxhopx2, maxhopy2, maxhopz2, dire, ik
    real, allocatable :: fermienergy(:) 
    real, allocatable :: deviation(:) 
    integer :: grid, i1, i2, i3, i4, orb, Nk1, Nk2, Nk3, knv3, ikx, iky, ikz
    real, allocatable :: conductivity(:), conductivity2(:), sigma_tensor_ahc_x(:)
    real, allocatable :: conductivity13(:), conductivity23(:) 
    real, allocatable :: conductivity_ahe(:), sigma_tensor_ahc_y(:), sigma_tensor_ahc_z(:)
    real, allocatable :: conductivity13_ahe(:), conductivity23_ahe(:) 
    real, allocatable :: conductivity_fsur(:) 
    real, allocatable :: conductivity13_fsur(:) 
    real, allocatable :: conductivity23_fsur(:), sigma_tensor_ahc_mpi(:,:), sigma_tensor_ahc_mpi2(:,:)
    real, allocatable :: sigma_tensor_shc_mpi(:,:), sigma_tensor_shc_mpi2(:,:), conductivity_she(:), conductivity13_she(:), conductivity23_she(:) 
    integer :: ierr, isize, irank, kp1, kp2, kp3, num_steps_tot
    integer :: ix1, ix2, ix3, num_occ 
    complex(kind=dp), allocatable :: momentum(:,:,:) 
    complex(kind=dp), allocatable :: momentum2(:,:,:) 

    complex(kind=dp), allocatable :: spinmomentum(:,:,:) 
    complex(kind=dp), allocatable :: spinmomentum2(:,:,:) 

    complex(kind=dp) :: berry 

    integer, allocatable :: sortarray(:), mag_wann_orbs_index(:)
    integer :: n1, n2, n3, n4, dir, mag_wann_num
    real, allocatable :: occupation(:), occupation2(:) 
    integer, allocatable :: nrpts(:) 
    real :: occupation_number 
    integer :: step, i, ii, num_steps 
    real :: fermi_min, fermi_max, efermi
    logical :: l_tb, l_nabla, l_bfield, l_mag_vec
    real :: bfield(3) 
    integer :: rvecnum, num_lines, length 
    integer, allocatable :: irvec(:,:)
    real, allocatable :: kpts(:,:), crvec(:,:)
    real :: kder, amat(3,3) 
    real :: volume, cross_vec(3), mag_field_x, mag_field_y, mag_field_z
    real :: bohrincm, condq, mag_strength, mag_theta, mag_phi
    integer :: maxdim, num_kpts 
    real :: minenerg, maxenerg, gamma, kbT
    logical :: l_bandstruc, l_fermisurf 
    real, allocatable :: magnetic(:,:) 
    complex(kind=dp), allocatable :: Omega_x(:), Omega_y(:), Omega_z(:), Omega_x_shc(:), Omega_y_shc(:), Omega_z_shc(:)
    complex(kind=dp), allocatable :: Omega_x_t(:), Omega_y_t(:), Omega_z_t(:)
    complex(kind=dp), allocatable :: Omega_x_t_shc(:), Omega_y_t_shc(:), Omega_z_t_shc(:)
    real(kind=dp) :: fermi

    real(kind=dp) :: K3D_start_cube(3)
    real(kind=dp) :: K3D_vec1_cube(3)
    real(kind=dp) :: K3D_vec2_cube(3)
    real(kind=dp) :: K3D_vec3_cube(3)

    real(kind=dp) :: pauli_result(4)
    real(kind=dp) :: trace_value

    real :: test_matrix(4, 4)

    integer, dimension(MPI_STATUS_SIZE) :: stt                                          
    call MPI_INIT(ierr) 
                                                                            
    call MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierr) 
    call MPI_COMM_SIZE(MPI_COMM_WORLD, isize, ierr) 
                                                                            
    abstol = 2.0_dp * tiny(abstol) 
    twopi = 2.0_dp * 3.141592654_dp 
                                                                            
    bohrincm = 0.5291772_dp * 1.0e-8_dp 
    condq = 38.7405_dp * 1.0e-6_dp 

    if (irank == 0) then 
        open(unit=300, file='ahe_inp') 
        read(300,*) amat(1,:) 
        read(300,*) amat(2,:) 
        read(300,*) amat(3,:) 
        read(300,*) fermi_min, fermi_max, num_steps 
        read(300,*) efermi
        read(300,*) Nk1, Nk2, Nk3 
        read(300,*) maxdim 
        read(300,*) occupation_number 
        read(300,*) temperature
        read(300,*) l_mag_vec         
        read(300,*) mag_strength                
        read(300,*) mag_theta                        
        read(300,*) mag_phi                          
        read(300,*) mag_wann_num
        allocate(mag_wann_orbs_index(mag_wann_num))
        read(300,*) mag_wann_orbs_index
        close(unit=300) 
        !print*, mag_wann_orbs_index
    endif 

    grid = Nk1
    kbT = 0.0861733_dp * 1.0e-3_dp * temperature
    amat = amat * 0.5291772_dp

    call MPI_BCAST(amat, 9, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)                             
    call MPI_BCAST(fermi_min, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)                             
    call MPI_BCAST(occupation_number, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)                             
    call MPI_BCAST(fermi_max, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)                             
    call MPI_BCAST(num_steps, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)                             
    call MPI_BCAST(maxdim, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)                             
    call MPI_BCAST(grid, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)    
    call MPI_BCAST(Nk1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)   
    call MPI_BCAST(Nk2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)   
    call MPI_BCAST(Nk3, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)                            
    call MPI_BCAST(bfield, 3, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)        
    call MPI_BCAST(l_mag_vec, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(mag_strength, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(mag_theta, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(mag_phi, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    mag_field_x = mag_strength * cos(mag_theta)
    mag_field_y = mag_strength * sin(mag_theta) * cos(mag_phi)
    mag_field_z = mag_strength * sin(mag_theta) * sin(mag_phi)

    call MPI_BCAST(mag_field_x, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(mag_field_y, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(mag_field_z, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
          
    if (irank == 0) then 
        open(unit=200, file='hopping.1') 
        num_lines = 0 
        num_wann = 0 
        do 
            read(200, fmt=*, end=311) ix, iy, iz, band1, band2, rdum, idum 
            num_lines = num_lines + 1 
            num_wann = max(num_wann, band1) 
        end do 
    311 continue 
        rvecnum = num_lines / (num_wann * num_wann) 
        allocate(hops(1:num_wann, 1:num_wann, rvecnum)) 
        allocate(irvec(3, rvecnum)) 
        hops = (0.0_dp, 0.0_dp) 
        rewind(unit=200) 
        num_lines = 0 
        do 
            read(200, fmt=*, end=300) ix, iy, iz, band1, band2, rdum, idum 
            num_lines = num_lines + 1 
            rvecnum = (num_lines - 1) / (num_wann * num_wann) + 1 
            irvec(1, rvecnum) = ix 
            irvec(2, rvecnum) = iy 
            irvec(3, rvecnum) = iz 
            hops(band1, band2, rvecnum) = cmplx(rdum, idum) 
        end do 
    300 continue 
        close(unit=200) 

        allocate(crvec(3, rvecnum))
        do ii = 1, rvecnum
            crvec(:, ii) = amat(1, :) * irvec(1, ii) + amat(2, :) * irvec(2, ii) + amat(3, :) * irvec(3, ii)
        end do
            
        allocate(rspauli(1:num_wann, 1:num_wann, 3, rvecnum)) 
        open(unit=400, file='./rspauli.1') 
        num_lines = 0 
        do 
            read(400, fmt=*, end=500) ix, iy, iz, band1, band2, dir, rdum, idum 
            num_lines = num_lines + 1 
            rvecnum = (num_lines - 1) / (num_wann * num_wann * 3) + 1 
            rspauli(band1, band2, dir, rvecnum) = cmplx(rdum, idum) 
        end do 
    500 continue 
        close(unit=400)  

        allocate(nrpts(rvecnum)) 
        open(unit=14, file='nrpts_inp') 
        do j = 1, rvecnum / 15 
            read(14, '(15I5)') (nrpts(15*(j-1)+i), i = 1, 15) 
        end do 
        read(14, '(I5)') (nrpts(15*(rvecnum/15)+i), i = 1, mod(rvecnum, 15))           
        close(unit=14)                                   
    endif 
                                                                         
    call MPI_BCAST(num_wann, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)                             
    call MPI_BCAST(rvecnum, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)                             
                                                                            
    if (.not. allocated(nrpts)) allocate(nrpts(rvecnum)) 
    call MPI_BCAST(nrpts, rvecnum, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)                             
                                                                            
    if (.not. allocated(irvec)) allocate(irvec(3, rvecnum)) 
    length = 3 * rvecnum 
    call MPI_BCAST(irvec, length, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)                             
            
    if (.not. allocated(crvec)) allocate(crvec(3, rvecnum))
    call MPI_BCAST(crvec, length, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    
    if (.not. allocated(hops)) allocate(hops(num_wann, num_wann, rvecnum)) 
    length = num_wann * num_wann * rvecnum 
    call MPI_BCAST(hops, length, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)                             
          
    if (.not. allocated(pauli_matrix)) allocate(pauli_matrix(num_wann, num_wann, 3))
    length = num_wann * num_wann * 3
    call MPI_BCAST(pauli_matrix, length, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
    
    if (.not. allocated(rspauli)) allocate(rspauli(num_wann, num_wann, 3, rvecnum))
    length = num_wann * num_wann * 3 * rvecnum
    call MPI_BCAST(rspauli, length, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)
    
    if (l_bfield) then 
        allocate(magnetic(3, num_steps)) 
        allocate(paulifft(num_wann, num_wann, 3)) 
        allocate(paulifft2(num_wann, num_wann, 3)) 
        if (.not. allocated(rspauli)) allocate(rspauli(num_wann, num_wann, 3, rvecnum))              
        length = num_wann * num_wann * rvecnum * 3 
        call MPI_BCAST(rspauli, length, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)                             
    endif 
                                                                                                

    cross_vec(1) = amat(1, 2) * amat(2, 3) - amat(1, 3) * amat(2, 2) 
    cross_vec(2) = amat(1, 3) * amat(2, 1) - amat(1, 1) * amat(2, 3) 
    cross_vec(3) = amat(1, 1) * amat(2, 2) - amat(1, 2) * amat(2, 1) 
                                                                            
    volume = cross_vec(1) * amat(3, 1) + cross_vec(2) * amat(3, 2) + cross_vec(3) * amat(3, 3) 

    K3D_start_cube = (/ 0.0_dp,  0.0_dp,   0.0_dp /)
    K3D_vec1_cube = (/ 1.0_dp,  0.0_dp,   0.0_dp /)
    K3D_vec2_cube = (/ 0.0_dp,  1.0_dp,   0.0_dp /)
    K3D_vec3_cube = (/ 0.0_dp,  0.0_dp,   1.0_dp /)    

    ! 初始化测试矩阵
    test_matrix = reshape((/1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, &
                           5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp, &
                           9.0_dp, 10.0_dp, 11.0_dp, 12.0_dp, &
                           13.0_dp, 14.0_dp, 15.0_dp, 16.0_dp/), (/4, 4/))
    print *, "The 4x4 test_matrix is:"
    print *, test_matrix
    call MPI_BCAST(test_matrix, size(test_matrix), MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! 分配数组
    allocate(sortarray(num_steps)) 
    allocate(deviation(num_steps)) 
    allocate(occupation(num_steps)) 
    allocate(occupation2(num_steps)) 
    allocate(conductivity(num_steps)) 
    allocate(conductivity2(num_steps)) 
    allocate(conductivity13(num_steps)) 
    allocate(conductivity23(num_steps)) 
    allocate(conductivity_ahe(num_steps)) 
    allocate(conductivity13_ahe(num_steps)) 
    allocate(conductivity23_ahe(num_steps)) 
    allocate(conductivity_she(num_steps)) 
    allocate(conductivity13_she(num_steps)) 
    allocate(conductivity23_she(num_steps)) 
    allocate(fermienergy(num_steps)) 
          
    allocate(sigma_tensor_ahc_x(num_steps))
    allocate(sigma_tensor_ahc_y(num_steps))
    allocate(sigma_tensor_ahc_z(num_steps))
    allocate(sigma_tensor_ahc_mpi(3, num_steps))
    allocate(sigma_tensor_ahc_mpi2(3, num_steps))
    allocate(Omega_x(num_wann), Omega_y(num_wann), Omega_z(num_wann))
    allocate(Omega_x_t(num_wann), Omega_y_t(num_wann), Omega_z_t(num_wann))
    allocate(sigma_tensor_shc_mpi(3, num_steps))
    allocate(sigma_tensor_shc_mpi2(3, num_steps))
    allocate(Omega_x_shc(num_wann), Omega_y_shc(num_wann), Omega_z_shc(num_wann))
    allocate(Omega_x_t_shc(num_wann), Omega_y_t_shc(num_wann), Omega_z_t_shc(num_wann))

    allocate(rho(num_wann, num_wann))
    allocate(rho_mpi(num_wann, num_wann)) 

    rho = (0.0_dp, 0.0_dp)
    rho_mpi = (0.0_dp, 0.0_dp)

    magnetic = 0.0_dp 
    occupation = 0.0_dp 
    conductivity = 0.0_dp 
    conductivity13 = 0.0_dp 
    conductivity23 = 0.0_dp 
    conductivity_ahe = 0.0_dp 
    conductivity13_ahe = 0.0_dp 
    conductivity23_ahe = 0.0_dp 

    conductivity_she = 0.0_dp 
    conductivity13_she = 0.0_dp 
    conductivity23_she = 0.0_dp 

    allocate(ham(num_wann, num_wann)) 
    ! allocate(pauli_matrix(num_wann, num_wann, 3)) 
    allocate(momentum(num_wann, num_wann, 3)) 
    allocate(momentum2(num_wann, num_wann, 3)) 
    allocate(spinmomentum(num_wann, num_wann, 3)) 
    allocate(spinmomentum2(num_wann, num_wann, 3))    

    allocate(eigvals(num_wann)) 
    allocate(eigvecs(num_wann, num_wann)) 
    allocate(eigvecs_dag(num_wann, num_wann))
    allocate(mat_temp(num_wann, num_wann))
    lwork = 12 * num_wann 
    allocate(work(lwork)) 
    allocate(rwork(17 * num_wann)) 
    allocate(iwork(15 * num_wann)) 
    allocate(ifail(15 * num_wann)) 

    if (irank == 0) then
        if (l_mag_vec) then
            do ii = 1, rvecnum
                do i = 1, num_wann
                    do j = 1, num_wann
                        hops(j, i, ii) = hops(j, i, ii) + &
                            mag_field_x * 0.5_dp * rspauli(j, i, 1, ii) + &
                            mag_field_y * 0.5_dp * rspauli(j, i, 2, ii) + &
                            mag_field_z * 0.5_dp * rspauli(j, i, 3, ii)
                    end do ! j
                end do ! i
            end do ! ii
        endif
    endif

    ! time_start = 0.0
    knv3 = Nk1 * Nk2 * Nk3
    ! call now(time_start)
    do ik = 1 + irank, knv3, isize
        ! 时间估算的注释代码
        ! if (irank == 0 .and. mod(ik/isize, 1) == 0) then
        !     call now(time_end) 
        !     write(*, '(a, i18, "/", i18, a, f10.2, "s")') 'ik/knv3', &
        !          ik, knv3, '  time left', (knv3-ik)*(time_end-time_start)/isize
        !     time_start = time_end
        ! endif

        ikx = (ik - 1) / (Nk2 * Nk3) + 1
        iky = ((ik - 1 - (ikx - 1) * Nk2 * Nk3) / Nk3) + 1
        ikz = ik - (iky - 1) * Nk3 - (ikx - 1) * Nk2 * Nk3
        kpoints = K3D_start_cube + K3D_vec1_cube * (ikx - 1) / dble(Nk1)  &
                  + K3D_vec2_cube * (iky - 1) / dble(Nk2)  &
                  + K3D_vec3_cube * (ikz - 1) / dble(Nk3)

        ham = (0.0_dp, 0.0_dp) 
        pauli_matrix = (0.0_dp, 0.0_dp) 
        Beta_fake = 1.0_dp / kbT                                                        
        do ii = 1, rvecnum 
            ix = irvec(1, ii) 
            iy = irvec(2, ii) 
            iz = irvec(3, ii) 
            phas = iy * kpoints(2) 
            phas = phas + iz * kpoints(3) 
            phas = phas + ix * kpoints(1) 
            phas = phas * twopi 
            fac = cmplx(cos(phas), sin(phas)) 
            ham(:,:) = ham(:,:) + fac * hops(:,:,ii) / nrpts(ii)               
            pauli_matrix(:,:,1) = pauli_matrix(:,:,1) + fac * rspauli(:,:,1,ii)
            pauli_matrix(:,:,2) = pauli_matrix(:,:,2) + fac * rspauli(:,:,2,ii)
            pauli_matrix(:,:,3) = pauli_matrix(:,:,3) + fac * rspauli(:,:,3,ii)
        end do 

        call zheevx('V', 'A', 'U', num_wann, ham, num_wann, vl, vu, 1, num_wann, abstol, ne, &
                   eigvals, eigvecs, num_wann, work, lwork, rwork, iwork, ifail, info)                        
        if (info /= 0) stop 'zheevx'                              

        eigvecs_dag = conjg(transpose(eigvecs))
        rho = rho + matmul(eigvecs * fermi(eigvals - efermi, Beta_fake), eigvecs_dag) * (1.0_dp / knv3)
        write(*,*) "knv3", knv3
        write(*,*) "eigvecs= ", eigvecs 
        write(*,*) "eigvals= ", eigvals
        write(*,*) "efermi= ", efermi
        write(*,*) "Beta_fake= ", Beta_fake
        write(*,*) "fermi(eigvals-efermi, Beta_fake) is ", fermi(eigvals - efermi, Beta_fake)
        write(*,*) "irank = ", irank, "rho is ", rho 

        momentum = (0.0_dp, 0.0_dp) 
        do ii = 1, rvecnum 
            ix = irvec(1, ii) 
            iy = irvec(2, ii) 
            iz = irvec(3, ii) 
            phas = iy * kpoints(2) 
            phas = phas + iz * kpoints(3) 
            phas = phas + ix * kpoints(1) 
            phas = phas * twopi                       
            fac = cmplx(-sin(phas), cos(phas))                    
            momentum(:,:,1) = momentum(:,:,1) + fac * crvec(1, ii) * hops(:,:,ii) / nrpts(ii) 
            momentum(:,:,2) = momentum(:,:,2) + fac * crvec(2, ii) * hops(:,:,ii) / nrpts(ii)
            momentum(:,:,3) = momentum(:,:,3) + fac * crvec(3, ii) * hops(:,:,ii) / nrpts(ii)
        end do 

        momentum2 = (0.0_dp, 0.0_dp)
        call mat_mul(num_wann, momentum(:,:,1), eigvecs, mat_temp)
        call mat_mul(num_wann, eigvecs_dag, mat_temp, momentum2(:,:,1))
        call mat_mul(num_wann, momentum(:,:,2), eigvecs, mat_temp)
        call mat_mul(num_wann, eigvecs_dag, mat_temp, momentum2(:,:,2))
        call mat_mul(num_wann, momentum(:,:,3), eigvecs, mat_temp)
        call mat_mul(num_wann, eigvecs_dag, mat_temp, momentum2(:,:,3))

        spinmomentum = (0.0_dp, 0.0_dp)
        do dir = 1, 3
            spinmomentum(:,:,dir) = (matmul(pauli_matrix(:,:,3), momentum(:,:,dir)) + matmul(momentum(:,:,dir), pauli_matrix(:,:,3))) / 2.0_dp
        end do

        spinmomentum2 = (0.0_dp, 0.0_dp)
        call mat_mul(num_wann, spinmomentum(:,:,1), eigvecs, mat_temp)
        call mat_mul(num_wann, eigvecs_dag, mat_temp, spinmomentum2(:,:,1))
        call mat_mul(num_wann, spinmomentum(:,:,2), eigvecs, mat_temp)
        call mat_mul(num_wann, eigvecs_dag, mat_temp, spinmomentum2(:,:,2))
        call mat_mul(num_wann, spinmomentum(:,:,3), eigvecs, mat_temp)
        call mat_mul(num_wann, eigvecs_dag, mat_temp, spinmomentum2(:,:,3))

        Omega_x = (0.0_dp, 0.0_dp)
        Omega_y = (0.0_dp, 0.0_dp)
        Omega_z = (0.0_dp, 0.0_dp)
        do m = 1, num_wann
            do n = 1, num_wann
                if (abs(eigvals(m) - eigvals(n)) > 1.0e-8_dp) then
                    Omega_z(m) = Omega_z(m) - aimag(momentum2(m,n,1) * conjg(momentum2(m,n,2))) / (eigvals(m) - eigvals(n))**2
                    Omega_y(m) = Omega_y(m) - aimag(momentum2(m,n,1) * conjg(momentum2(m,n,3))) / (eigvals(m) - eigvals(n))**2
                    Omega_x(m) = Omega_x(m) - aimag(momentum2(m,n,2) * conjg(momentum2(m,n,3))) / (eigvals(m) - eigvals(n))**2
                endif
            end do ! n
        end do ! m

        Omega_x = Omega_x * 2.0_dp
        Omega_y = Omega_y * 2.0_dp
        Omega_z = Omega_z * 2.0_dp   
             
        do step = 1, num_steps                                                             
            fermienergy(step) = fermi_min + (fermi_max - fermi_min) * real(step) / real(num_steps)        
            mu = fermienergy(step)

            do m = 1, num_wann
                Omega_x_t(m) = Omega_x(m) * fermi(eigvals(m) - mu, Beta_fake)
                Omega_y_t(m) = Omega_y(m) * fermi(eigvals(m) - mu, Beta_fake)
                Omega_z_t(m) = Omega_z(m) * fermi(eigvals(m) - mu, Beta_fake)
            end do
            sigma_tensor_ahc_mpi(1, step) = sigma_tensor_ahc_mpi(1, step) + real(sum(Omega_x_t))
            sigma_tensor_ahc_mpi(2, step) = sigma_tensor_ahc_mpi(2, step) + real(sum(Omega_y_t))
            sigma_tensor_ahc_mpi(3, step) = sigma_tensor_ahc_mpi(3, step) + real(sum(Omega_z_t))
        end do

        Omega_x_shc = (0.0_dp, 0.0_dp)
        Omega_y_shc = (0.0_dp, 0.0_dp)
        Omega_z_shc = (0.0_dp, 0.0_dp)
        do n = 23, num_wann
            do m = 1, 22
                if (abs(eigvals(m) - eigvals(n)) > 1.0e-8_dp) then
                    Omega_z_shc(m) = Omega_x_shc(m) - aimag(momentum2(m,n,1) * conjg(spinmomentum2(m,n,2))) / (eigvals(m) - eigvals(n))**2
                    Omega_y_shc(m) = Omega_y_shc(m) - aimag(momentum2(m,n,1) * conjg(spinmomentum2(m,n,3))) / (eigvals(m) - eigvals(n))**2
                    Omega_x_shc(m) = Omega_z_shc(m) - aimag(momentum2(m,n,2) * conjg(spinmomentum2(m,n,3))) / (eigvals(m) - eigvals(n))**2
                endif
            end do ! n
        end do ! m

        Omega_x_shc = Omega_x_shc * 2.0_dp
        Omega_y_shc = Omega_y_shc * 2.0_dp
        Omega_z_shc = Omega_z_shc * 2.0_dp   
             
        do step = 1, num_steps                                                             
            fermienergy(step) = fermi_min + (fermi_max - fermi_min) * real(step) / real(num_steps)        
            mu = fermienergy(step)

            do m = 1, num_wann         
                Omega_x_t_shc(m) = Omega_x_shc(m)  !* fermi(eigvals(m) - mu, Beta_fake)
                Omega_y_t_shc(m) = Omega_y_shc(m)  !* fermi(eigvals(m) - mu, Beta_fake)
                Omega_z_t_shc(m) = Omega_z_shc(m)  !* fermi(eigvals(m) - mu, Beta_fake)
            end do
            sigma_tensor_shc_mpi(1, step) = sigma_tensor_shc_mpi(1, step) + real(sum(Omega_x_t_shc))
            sigma_tensor_shc_mpi(2, step) = sigma_tensor_shc_mpi(2, step) + real(sum(Omega_y_t_shc))
            sigma_tensor_shc_mpi(3, step) = sigma_tensor_shc_mpi(3, step) + real(sum(Omega_z_t_shc))
        end do

    end do

    num_steps_tot = num_steps * 3       
       
    call MPI_REDUCE(rho, rho_mpi, size(rho), MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, ierr) 

    call MPI_REDUCE(sigma_tensor_ahc_mpi, sigma_tensor_ahc_mpi2, size(sigma_tensor_ahc_mpi), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call MPI_REDUCE(sigma_tensor_shc_mpi, sigma_tensor_shc_mpi2, size(sigma_tensor_shc_mpi), MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    if (irank == 0) then
        sigma_tensor_ahc_mpi = sigma_tensor_ahc_mpi2
        sigma_tensor_ahc_mpi = sigma_tensor_ahc_mpi / knv3 / volume * condq * 1.0e8_dp * twopi
    
        sigma_tensor_shc_mpi = sigma_tensor_shc_mpi2
        sigma_tensor_shc_mpi = sigma_tensor_shc_mpi / knv3 / volume * condq * 1.0e8_dp * twopi

        conductivity_ahe = sigma_tensor_ahc_mpi(3, :) 
        conductivity13_ahe = sigma_tensor_ahc_mpi(2, :)   
        conductivity23_ahe = sigma_tensor_ahc_mpi(1, :)
        
        conductivity_she = sigma_tensor_shc_mpi(3, :) 
        conductivity13_she = sigma_tensor_shc_mpi(2, :)   
        conductivity23_she = sigma_tensor_shc_mpi(1, :)
    
        write(*,*) "rho_mpi = ", rho_mpi 
        call pauli_block_all(real(rho_mpi), num_wann, pauli_result)

        open(unit=123, file='output_ahe_condquant', recl=10000) 
        do step = 1, num_steps 
            write(123,*) "fermienergy=", fermienergy(step)                                                          
            write(123,*) "occupation=", occupation(step)                                
            write(123,*) "conductivity=", conductivity_ahe(step) / condq * amat(3,3) * 1.0e-8  ! 2D direct convert to integer
            write(123,*) "conductivity13=", conductivity13_ahe(step) / condq * amat(3,3) * 1.0e-8                 
            write(123,*) "conductivity23=", conductivity23_ahe(step) / condq * amat(3,3) * 1.0e-8
            ! write(123,*) "conductivity=", conductivity_ahe(step)   ! unit is same with wannier_tool  S/cm
            ! write(123,*) "conductivity13=", conductivity13_ahe(step)
            ! write(123,*) "conductivity23=", conductivity23_ahe(step)                      
            write(123,*) "************************************"                                                  
        end do 
        close(unit=123)        
             
        open(unit=321, file='output_she_condquant', recl=10000) 
        do step = 1, num_steps 
            write(321,*) "fermienergy=", fermienergy(step)                                                          
            write(321,*) "occupation=", occupation(step)                                
            write(321,*) "conductivity=", conductivity_she(step) / condq * amat(3,3) * 1.0e-8  ! 2D direct convert to integer
            write(321,*) "conductivity13=", conductivity13_she(step) / condq * amat(3,3) * 1.0e-8                 
            write(321,*) "conductivity23=", conductivity23_she(step) / condq * amat(3,3) * 1.0e-8
            ! write(321,*) "conductivity_shc=", conductivity_she(step)   ! unit is same with wannier_tool  S/cm
            ! write(321,*) "conductivity13_shc=", conductivity13_she(step)
            ! write(321,*) "conductivity23_shc=", conductivity23_she(step)                      
            write(321,*) "************************************"                                                  
        end do 
        close(unit=321)  
    
        open(unit=458, file='charge_on_each_atom', recl=10000) 
        do m = 1, 4
            write(458,*) "m= ", m
            write(458,*) "charge=", pauli_result(m)
            write(458,*) "************************************" 
        end do
        close(unit=458) 
    endif 

    call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
    call MPI_FINALIZE(ierr) 

end program ahc_shc_anc_hao

! fermi函数定义
function fermi(omega, Beta_fake) result(value)
    use pauli  ! 使用 pauli 模块以获取 dp
    implicit none
    real(kind=dp), intent(in) :: omega
    real(kind=dp), intent(in) :: Beta_fake
    real(kind=dp) :: value 
    if (Beta_fake * omega >= 20.0_dp) then
        value = 0.0_dp
    elseif (Beta_fake * omega <= -20.0_dp) then
        value = 1.0_dp
    else
        value = 1.0_dp / (1.0_dp + exp(Beta_fake * omega))
    endif
    return
end function fermi    

! 时间获取子程序
subroutine now(time_now)
    implicit none
    integer :: time_new(8)
    real(kind=dp), intent(out) :: time_now
    call Date_and_time(values=time_new)
    time_now = time_new(3) * 24.0_dp * 3600.0_dp + time_new(5) * 3600.0_dp + &
               time_new(6) * 60.0_dp + time_new(7) + time_new(8) / 1000.0_dp  
    return
end subroutine now

! 矩阵乘法子程序
subroutine mat_mul(nmatdim, A, B, C)
    use pauli  ! 使用 pauli 模块以获取 dp
    implicit none
    integer, intent(in) :: nmatdim
    complex(kind=dp), intent(in) :: A(nmatdim, nmatdim)
    complex(kind=dp), intent(in) :: B(nmatdim, nmatdim)
    complex(kind=dp), intent(out) :: C(nmatdim, nmatdim)
    complex(kind=dp) :: ALPHA, BETA

    ALPHA = (1.0_dp, 0.0_dp)
    BETA = (0.0_dp, 0.0_dp)

    C = (0.0_dp, 0.0_dp)

    call ZGEMM('N', 'N', nmatdim, nmatdim, nmatdim, ALPHA, A, nmatdim, B, nmatdim, BETA, C, nmatdim)

    return
end subroutine mat_mul
