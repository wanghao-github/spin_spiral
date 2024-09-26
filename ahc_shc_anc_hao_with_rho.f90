      program ahc_shc_anc_hao
!*************************************************                      
!rewrite by HaoWang March 7 2023 replace momentum2 by wt
!*************************************************     
      use mpi     
      use pauli_comp
      use fermi_module            
      implicit none 
                                                                        
      complex,allocatable:: hops   (:,:,:) 
      complex,allocatable:: rsnabla(:,:,:,:) 
      complex,allocatable:: rspauli(:,:,:,:) 
      complex,allocatable:: pauli(:,:,:) 
      complex,allocatable:: paulifft(:,:,:) 
      complex,allocatable:: paulifft2(:,:,:) 
      complex,allocatable:: rho(:,:),rho_mpi(:,:)
      real               :: rdum,idum 
      integer            :: ix,iy,iz,band1,band2,h,num_wann,m,n
      integer            :: ik1,ik2,ik3 
      real               :: twopi,temperature
      real               :: phas 
      complex            :: fac,fac2 
      complex,allocatable:: ham(:,:) 
      real               :: vl,vu 
      integer            :: ne,j 
      real               :: abstol,time_start,time_end ,time_start1,time_end1
      real,allocatable   :: eigvals(:) 
      complex,allocatable:: eigvecs(:,:) ,eigvecs_dag(:,:),mat_temp(:,:),eigvecs_f(:,:),temp(:,:)
      integer            :: info 
      complex,allocatable:: work(:) 
      integer            :: lwork 
      integer,allocatable:: iwork(:) 
      real,allocatable   :: rwork(:) 
      integer,allocatable:: ifail(:) 
      real               :: kpoints(3)
      real               :: scale ,Beta_fake,mu
      integer            :: maxhopx2,maxhopy2,maxhopz2,dire,ik
      real,allocatable   :: fermienergy(:) 
      real,allocatable   :: deviation(:) 
      integer            :: grid,i1,i2,i3,i4,orb,Nk1,Nk2,Nk3,knv3,ikx,iky,ikz
      real,allocatable   :: conductivity(:),conductivity2(:),sigma_tensor_ahc_x(:)
      real,allocatable   :: conductivity13(:),conductivity23(:) 
      real,allocatable   :: conductivity_ahe(:) ,sigma_tensor_ahc_y(:),sigma_tensor_ahc_z(:)
      real,allocatable   :: conductivity13_ahe(:),conductivity23_ahe(:) 
      real,allocatable   :: conductivity_fsur(:) 
      real,allocatable   :: conductivity13_fsur(:) 
      real,allocatable   :: conductivity23_fsur(:),sigma_tensor_ahc_mpi(:,:),sigma_tensor_ahc_mpi2(:,:)
      real,allocatable   :: sigma_tensor_shc_mpi(:,:),sigma_tensor_shc_mpi2(:,:),conductivity_she(:),conductivity13_she(:),conductivity23_she(:) 
      integer            :: ierr,isize,irank,kp1,kp2,kp3,num_steps_tot
      integer            :: ix1,ix2,ix3,num_occ 
      complex,allocatable:: momentum(:,:,:) 
      complex,allocatable:: momentum2(:,:,:) 
                                                                        
      complex,allocatable:: spinmomentum(:,:,:) 
      complex,allocatable:: spinmomentum2(:,:,:) 
      

      ! integer,parameter :: Dp=kind(1.0d0)      
      complex            :: berry 
                                                                        
      integer,allocatable:: sortarray(:),mag_wann_orbs_index(:)
      integer            :: n1,n2,n3,n4,dir,mag_wann_num
      real,allocatable   :: occupation(:),occupation2(:) 
      integer,allocatable   :: nrpts(:) 
      real               :: occupation_number 
      integer            :: step,i,ii,num_steps 
      real               :: fermi_min,fermi_max,efermi
      logical            :: l_tb,l_nabla,l_bfield,l_mag_vec
      real               :: bfield(3) 
      integer            :: rvecnum,num_lines,length 
      integer,allocatable:: irvec(:,:)
      real,allocatable   :: kpts(:,:) ,crvec(:,:)
      real               :: kder,amat(3,3) 
      real               :: volume,cross(3),mag_field_x,mag_field_y,mag_field_z
      real               :: bohrincm,condq,mag_strength,mag_theta,mag_phi
      integer            :: maxdim,num_kpts 
      real               :: minenerg,maxenerg,gamma,kbT
      logical            :: l_bandstruc,l_fermisurf 
      real,allocatable   :: magnetic(:,:) 
      complex(kind(1.0d0)), allocatable :: Omega_x(:), Omega_y(:), Omega_z(:), Omega_x_shc(:), Omega_y_shc(:), Omega_z_shc(:)
      complex(kind(1.0d0)), allocatable :: Omega_x_t(:), Omega_y_t(:), Omega_z_t(:)
      complex(kind(1.0d0)), allocatable :: Omega_x_t_shc(:), Omega_y_t_shc(:), Omega_z_t_shc(:)
      ! real(kind(1.0d0)) :: fermi
      ! real(kind(1.0d0)) :: fermi_array
      real(kind(1.0d0)) :: K3D_start_cube(3)
      real(kind(1.0d0)) :: K3D_vec1_cube(3)
      real(kind(1.0d0)) :: K3D_vec2_cube(3)
      real(kind(1.0d0)) :: K3D_vec3_cube(3)

      complex(kind(1.0d0)) :: pauli_result(4)
      complex(kind(1.0d0)) :: trace_value
      real(kind=8), dimension(:), allocatable :: fermi_values

      real :: test_matrix(4, 4)

      ! INCLUDE 'mpif.h' 
      ! integer stt(MPI_STATUS_SIZE) 
      integer, dimension(MPI_STATUS_SIZE) :: stt                                          
      CALL MPI_INIT(ierr) 
                                                                        
      CALL MPI_COMM_RANK (MPI_COMM_WORLD,irank,ierr) 
      CALL MPI_COMM_SIZE (MPI_COMM_WORLD,isize,ierr) 
                                                                        
                                                                        
      abstol=2.0*tiny(abstol) 
      twopi=2*3.141592654 
                                                                        
      bohrincm=0.5291772*1.e-8 
      condq=38.7405*1.e-6 
      if(irank.eq.0)then 
         open(300,file='ahe_inp') 
         read(300,*)amat(1,:) 
         read(300,*)amat(2,:) 
         read(300,*)amat(3,:) 
         read(300,*)fermi_min,fermi_max,num_steps 
         read(300,*)efermi
         read(300,*)Nk1,Nk2,Nk3 
         read(300,*)maxdim 
         read(300,*)occupation_number 
         read(300,*)temperature
         read(300,*)l_mag_vec         
         read(300,*)mag_strength                
         read(300,*)mag_theta                        
         read(300,*)mag_phi                          
         read(300,*)mag_wann_num
         allocate(mag_wann_orbs_index(mag_wann_num))
         read(300,*)mag_wann_orbs_index(mag_wann_num)
         close(300) 
         !print*, mag_wann_orbs_index
      endif 
      grid = Nk1
      kbT = 0.0861733*1e-3*temperature
      amat= amat*0.5291772
      call mpi_bcast(amat,9,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr)                             
      call mpi_bcast(fermi_min,1,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr)                             
      call mpi_bcast(occupation_number,1,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr)                             
      call mpi_bcast(fermi_max,1,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr)    
      call mpi_bcast(efermi,1,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr)                     
      call mpi_bcast(num_steps,1,MPI_INTEGER, 0,mpi_comm_world,ierr)                             
      call mpi_bcast(maxdim,1,MPI_INTEGER, 0,mpi_comm_world,ierr)                             
      call mpi_bcast(grid,1,MPI_INTEGER,0,mpi_comm_world,ierr)    
      call mpi_bcast(Nk1,1,MPI_INTEGER, 0,mpi_comm_world,ierr)   
      call mpi_bcast(Nk2,1,MPI_INTEGER, 0,mpi_comm_world,ierr)   
      call mpi_bcast(Nk3,1,MPI_INTEGER, 0,mpi_comm_world,ierr)                            
      !call mpi_bcast(l_tb,1,MPI_LOGICAL,0,mpi_comm_world,ierr)                             
      !call mpi_bcast(l_nabla,1,MPI_LOGICAL,0,mpi_comm_world,ierr)                             
      !call mpi_bcast(l_bfield,1,MPI_LOGICAL, 0,mpi_comm_world,ierr)                             
      !call mpi_bcast(l_fermisurf,1,MPI_LOGICAL,0,mpi_comm_world,ierr)                             
      !call mpi_bcast(bfield,3,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)                             
      !call mpi_bcast(gamma,1,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr)                             
      call mpi_bcast(bfield,3,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)        
      call mpi_bcast(l_mag_vec,1,MPI_LOGICAL,0,mpi_comm_world,ierr)
      call mpi_bcast(mag_strength,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)
      call mpi_bcast(mag_theta,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)
      call mpi_bcast(mag_phi,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)
      !call mpi_bcast(mag_wann_num,1,MPI_INTEGER, 0,mpi_comm_world,ierr)
      !call mpi_bcast(mag_wann_orbs_index,mag_wann_num,MPI_INTEGER,0,mpi_comm_world,ierr)
      call mpi_bcast(temperature,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)
      call mpi_bcast(kbT,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)

      mag_field_x = mag_strength * cos(mag_theta)
      mag_field_y = mag_strength * sin(mag_theta) * cos(mag_phi)
      mag_field_z = mag_strength * sin(mag_theta) * sin(mag_phi)

      call mpi_bcast(mag_field_x,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)
      call mpi_bcast(mag_field_y,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)
      call mpi_bcast(mag_field_z,1,MPI_DOUBLE_PRECISION, 0,mpi_comm_world,ierr)
      
      if(irank.eq.0)then 
         open(200,file='hopping.1') 
         num_lines=0 
         num_wann=0 
         do 
            read(200,fmt=*,end=311)ix,iy,iz,band1,band2,rdum,idum 
            num_lines=num_lines+1 
            num_wann=max(num_wann,band1) 
         enddo 
  311    continue 
         rvecnum=num_lines/(num_wann*num_wann) 
         allocate( hops(1:num_wann,1:num_wann,rvecnum) ) 
         allocate( irvec(3,rvecnum) ) 
         hops=0.0 
         rewind(200) 
         num_lines=0 
         do 
            read(200,fmt=*,end=300)ix,iy,iz,band1,band2,rdum,idum 
            num_lines=num_lines+1 
            rvecnum=(num_lines-1)/(num_wann*num_wann)+1 
            irvec(1,rvecnum)=ix 
            irvec(2,rvecnum)=iy 
            irvec(3,rvecnum)=iz 
            hops( band1,band2,rvecnum )=cmplx(rdum,idum) 
         enddo 
  300    continue 
         close(200) 

        allocate(crvec(3,rvecnum))
        do ii=1,rvecnum
           crvec(:,ii)= amat(1,:)*irvec(1,ii)+amat(2,:)* irvec(2,ii)+amat(3,:)* irvec(3,ii)
        enddo
        
        allocate(rspauli(1:num_wann, 1:num_wann, 3, rvecnum)) 
        open(400,file='./rspauli.1') 
        num_lines=0 
        Do 
           read(400, fmt=*,end=500) ix,iy,iz,band1,band2,dir,rdum,idum 
           num_lines=num_lines+1 
           rvecnum=(num_lines-1)/(num_wann*num_wann*3)+1 
           rspauli(band1, band2, dir, rvecnum)=cmplx(rdum,idum) 
        End Do 
  500   continue 
        close(400)  
        allocate(nrpts(rvecnum)) 
        open(14,file='nrpts_inp') 
        do j=1,rvecnum/15 
                read(14,'(15I5)') (nrpts(15*(j-1)+i) ,i=1,15) 
        enddo 
        read(14,'(<mod(rvecnum,15)>I5)') (nrpts(15*(rvecnum/15)+i),i=1,mod(rvecnum,15))           
        close(14)                                   
      endif 
                                                                     
                                                              
      call mpi_bcast(num_wann,1,MPI_INTEGER,0,mpi_comm_world,ierr)                             
      call mpi_bcast(rvecnum,1,MPI_INTEGER,0,mpi_comm_world,ierr)                             
                                                                        
                                                                        
      if(.not.allocated(nrpts))allocate(nrpts(rvecnum)) 
      call mpi_bcast(nrpts,rvecnum,MPI_INTEGER,0,mpi_comm_world,ierr)                             
                                                                        
      if(.not.allocated(irvec))allocate(irvec(3,rvecnum)) 
      length=3*rvecnum 
      call mpi_bcast(irvec,length,MPI_INTEGER,0,mpi_comm_world,ierr)                             
        
      if(.not.allocated(crvec))allocate(crvec(3,rvecnum))
      call mpi_bcast(crvec,length,MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr)

      if(.not.allocated(hops))allocate(hops(num_wann,num_wann,rvecnum)) 
      length=num_wann*num_wann*rvecnum 
      call mpi_bcast(hops,length,MPI_DOUBLE_COMPLEX,0,mpi_comm_world,ierr)                             
      
      if(.not.allocated(pauli))allocate(pauli(num_wann,num_wann,3))
      length=num_wann*num_wann*3
      call mpi_bcast(pauli,length,MPI_DOUBLE_COMPLEX,0,mpi_comm_world,ierr)


      if(.not.allocated(rspauli)) allocate(rspauli(num_wann,num_wann,3,rvecnum))
      length=num_wann*num_wann*rvecnum*3
      call mpi_bcast(rspauli,length,MPI_DOUBLE_COMPLEX,0,mpi_comm_world,ierr)

      if(l_bfield)then 
         allocate(magnetic(3,num_steps)) 
         allocate(paulifft(num_wann,num_wann,3)) 
         allocate(paulifft2(num_wann,num_wann,3)) 
         if(.not.allocated(rspauli)) allocate(rspauli(num_wann,num_wann,3,rvecnum))              
         length=num_wann*num_wann*rvecnum*3 
         call mpi_bcast(rspauli,length,MPI_DOUBLE_COMPLEX,0,mpi_comm_world,ierr)                             
      endif 
                                                                                                                        
      cross(1)=amat(1,2)*amat(2,3)-amat(1,3)*amat(2,2) 
      cross(2)=amat(1,3)*amat(2,1)-amat(1,1)*amat(2,3) 
      cross(3)=amat(1,1)*amat(2,2)-amat(1,2)*amat(2,1) 
                                                                        
      volume=cross(1)*amat(3,1)+cross(2)*amat(3,2)+cross(3)*amat(3,3) 

      K3D_start_cube= (/ 0.0,  0.0,   0.0/)
      K3D_vec1_cube = (/ 1.0,  0.0,   0.0/)
      K3D_vec2_cube = (/ 0.0,  1.0,   0.0/)
      K3D_vec3_cube = (/ 0.0,  0.0,   1.0/)    


      test_matrix = reshape((/1.0, 2.0, 3.0, 4.0, &
      5.0, 6.0, 7.0, 8.0, &
      9.0, 10.0, 11.0, 12.0, &
      13.0, 14.0, 15.0, 16.0/), (/4, 4/))
      ! print *, "The 4x4 test_matrix is:"
      ! print *, test_matrix
      call mpi_bcast(test_matrix,size(test_matrix),MPI_DOUBLE_PRECISION,0,mpi_comm_world,ierr)

      allocate( sortarray(num_steps) ) 
      allocate( deviation(num_steps) ) 
      allocate( occupation(num_steps) ) 
      allocate( occupation2(num_steps) ) 
      allocate( conductivity(num_steps) ) 
      allocate( conductivity2(num_steps) ) 
      allocate( conductivity13(num_steps) ) 
      allocate( conductivity23(num_steps) ) 
      allocate( conductivity_ahe(num_steps) ) 
      allocate( conductivity13_ahe(num_steps) ) 
      allocate( conductivity23_ahe(num_steps) ) 
      allocate( conductivity_she(num_steps) ) 
      allocate( conductivity13_she(num_steps) ) 
      allocate( conductivity23_she(num_steps) ) 
      allocate( fermienergy(num_steps) ) 
      
      allocate(sigma_tensor_ahc_x(num_steps))
      allocate(sigma_tensor_ahc_y(num_steps))
      allocate(sigma_tensor_ahc_z(num_steps))
      allocate(sigma_tensor_ahc_mpi(3,num_steps))
      allocate(sigma_tensor_ahc_mpi2(3,num_steps))
      allocate(Omega_x(num_wann), Omega_y(num_wann), Omega_z(num_wann))
      allocate(Omega_x_t(num_wann), Omega_y_t(num_wann), Omega_z_t(num_wann))
      allocate(sigma_tensor_shc_mpi(3,num_steps))
      allocate(sigma_tensor_shc_mpi2(3,num_steps))
      allocate(Omega_x_shc(num_wann), Omega_y_shc(num_wann), Omega_z_shc(num_wann))
      allocate(Omega_x_t_shc(num_wann), Omega_y_t_shc(num_wann), Omega_z_t_shc(num_wann))

      allocate(rho(num_wann,num_wann))
      allocate(rho_mpi(num_wann,num_wann)) 


      rho = 0.0
      rho_mpi  = 0.0


      magnetic=0.0 
      occupation=0.0 
      conductivity=0.0 
      conductivity13=0.0 
      conductivity23=0.0 
      conductivity_ahe=0.0 
      conductivity13_ahe=0.0 
      conductivity23_ahe=0.0 

      conductivity_she=0.0 
      conductivity13_she=0.0 
      conductivity23_she=0.0 

                                                                        
      allocate(        ham(num_wann,num_wann)) 
     ! allocate(    pauli(num_wann,num_wann,3)) 
      allocate(momentum (num_wann,num_wann,3)) 
      allocate(momentum2(num_wann,num_wann,3)) 
      allocate(spinmomentum (num_wann,num_wann,3)) 
      allocate(spinmomentum2(num_wann,num_wann,3))    
      

      allocate(eigvals(num_wann)) 
      allocate(eigvecs(num_wann,num_wann)) 
      allocate(eigvecs_f(num_wann,num_wann)) 
      allocate(temp(num_wann,num_wann)) 
      allocate(eigvecs_dag(num_wann,num_wann))
      allocate(mat_temp(num_wann,num_wann))
      lwork=12.0*num_wann 
      allocate( work(lwork) ) 
      allocate( rwork(17*num_wann) ) 
      allocate( iwork(15*num_wann) ) 
      allocate( ifail(15*num_wann) ) 

      if(irank == 0)then
       if(l_mag_vec)then
        do ii=1,rvecnum
          do i=1,num_wann
            do j=1,num_wann
                hops(j,i,ii) = hops(j,i,ii) + mag_field_x * 0.5* rspauli(j,i,1,ii) + &
                mag_field_y * 0.5* rspauli(j,i,2,ii) + mag_field_z * 0.5* rspauli(j,i,3,ii)
            enddo !j
          enddo !i
        enddo
     endif
     endif

      time_start = 0.0
      knv3= Nk1*Nk2*Nk3
      call now(time_start)
       do ik= 1+ irank, knv3, isize
         if (irank .eq. 0 .and. mod(ik/isize, 1) .eq. 0) then
          call now(time_end) 
         !   write(*, '(a, i18, "/", i18, a, f10.2, "s")') 'ik/knv3', &
            ! ik, knv3, '  time left', (knv3-ik)*(time_end-time_start)/isize
            time_start= time_end
         endif

         ikx= (ik-1)/(nk2*nk3)+1
         iky= ((ik-1-(ikx-1)*Nk2*Nk3)/nk3)+1
         ikz= (ik-(iky-1)*Nk3- (ikx-1)*Nk2*Nk3)
         kpoints= K3D_start_cube+ K3D_vec1_cube*(ikx-1)/dble(nk1)  &
          + K3D_vec2_cube*(iky-1)/dble(nk2)  &
          + K3D_vec3_cube*(ikz-1)/dble(nk3)

         ham=0.0 
         pauli=cmplx(0.d0,0.d0) 
         Beta_fake = 1.0/kbT                                                        
         do ii=1,rvecnum 
            ix=irvec(1,ii) 
            iy=irvec(2,ii) 
            iz=irvec(3,ii) 
            phas=     iy*kpoints(2) 
            phas=phas+iz*kpoints(3) 
            phas=phas+ix*kpoints(1) 
            phas=phas*twopi 
            fac=cmplx(cos(phas),sin(phas)) 
            ham(:,:)=ham(:,:)+fac*hops(:,:,ii)/nrpts(ii)               
            pauli(:,:,1)=pauli(:,:,1)+fac*rspauli(:,:,1,ii)
            pauli(:,:,2)=pauli(:,:,2)+fac*rspauli(:,:,2,ii)
            pauli(:,:,3)=pauli(:,:,3)+fac*rspauli(:,:,3,ii)
         enddo 
 
      call zheevx('V','A','U',num_wann,ham,num_wann,vl,vu,1,num_wann,abstol,ne,&
                  eigvals,eigvecs,num_wann,work,lwork,rwork,iwork,ifail,info)                        
         if(info.ne.0)stop 'zheevx'                              
         

         eigvecs_dag=conjg(transpose(eigvecs))
         ! rho(:,:) = rho(:,:) +MATMUL(eigvecs * fermi(eigvals-efermi, Beta_fake), eigvecs_dag)*(1/knv3)
         ! rho(:,:) = rho(:,:) + MATMUL(eigvecs * fermi(eigvals-efermi, Beta_fake), eigvecs_dag)/knv3
         ! 确保 eigvecs 和 fermi_array 的维度匹配
         ! rho(:,:) = rho(:,:) + matmul(eigvecs * diag(fermi_array(eigvals-efermi, Beta_fake)), eigvecs_dag) / knv3
         
         fermi_values = fermi_array(eigvals-efermi, Beta_fake)
         print *, "fermi_values:", fermi_values
         print *, "here no problem0"
         do j = 1, num_wann
            do i = 1, num_wann
               eigvecs_f(i, j) = eigvecs(i, j) * fermi_values(j)
            end do
         end do

         print *, "here no problem1"
         temp = matmul(eigvecs_f, eigvecs_dag) / knv3
         print *, "here no problem2"
         write(*,*) "eigvecs_f:", eigvecs_f
         write(*,*) "temp:", temp
         rho = rho + temp
         ! rho(:,:) = rho(:,:) + MATMUL(eigvecs, eigvecs_dag)/knv3
         ! write(*,*) "no idea"
         ! write(*,*) "knv3", knv3
         ! write(*,*) "eigvecs= ", eigvecs 
         ! write(*,*) "eigvals= ", eigvals
         ! write(*,*) "eigvals-efermi", eigvals-efermi
         ! write(*,*) "efermi= ", efermi
         ! write(*,*) "Beta_fake= ", Beta_fake
         ! write(*,*) "fermi_array(eigvals-efermi, Beta_fake) is ", fermi_array(eigvals-efermi, Beta_fake)
         ! write(*,*) "irank = ", irank
         ! write(*,*) "1/knv3 is " , 1/knv3 
         ! write(*,*) "rho is " , rho 

         momentum=0.0 
         do ii=1,rvecnum 
            ix=irvec(1,ii) 
            iy=irvec(2,ii) 
            iz=irvec(3,ii) 
            phas=     iy*kpoints(2) 
            phas=phas+iz*kpoints(3) 
            phas=phas+ix*kpoints(1) 
            phas=phas*twopi                       
            fac=cmplx(-sin(phas),cos(phas))                    
            momentum(:,:,1)= momentum(:,:,1)+fac*crvec(1,ii)*hops(:,:,ii)/nrpts(ii) 
            momentum(:,:,2)= momentum(:,:,2)+fac*crvec(2,ii)*hops(:,:,ii)/nrpts(ii)
            momentum(:,:,3)= momentum(:,:,3)+fac*crvec(3,ii)*hops(:,:,ii)/nrpts(ii)
          enddo 

          momentum2= 0.0
          call mat_mul(num_wann,momentum(:,:,1),eigvecs,mat_temp)
          call mat_mul(num_wann,eigvecs_dag,mat_temp,momentum2(:,:,1))
          call mat_mul(num_wann,momentum(:,:,2),eigvecs,mat_temp)
          call mat_mul(num_wann,eigvecs_dag,mat_temp,momentum2(:,:,2))
          call mat_mul(num_wann,momentum(:,:,3),eigvecs,mat_temp)
          call mat_mul(num_wann,eigvecs_dag,mat_temp,momentum2(:,:,3))


         spinmomentum=0.0
         Do dir=1,3
            spinmomentum(:,:,dir)= (MATMUL(pauli(:,:,3),momentum(:,:,dir)) +MATMUL(momentum(:,:,dir),pauli(:,:,3)))/2.d0
         End Do
  
          spinmomentum2=0.0
          call mat_mul(num_wann,spinmomentum(:,:,1),eigvecs,mat_temp)
          call mat_mul(num_wann,eigvecs_dag,mat_temp,spinmomentum2(:,:,1))
          call mat_mul(num_wann,spinmomentum(:,:,2),eigvecs,mat_temp)
          call mat_mul(num_wann,eigvecs_dag,mat_temp,spinmomentum2(:,:,2))
          call mat_mul(num_wann,spinmomentum(:,:,3),eigvecs,mat_temp)
          call mat_mul(num_wann,eigvecs_dag,mat_temp,spinmomentum2(:,:,3))

          

         Omega_x=0d0;Omega_y=0d0; Omega_z=0d0
         do m=1,num_wann
            do n=1,num_wann
               if (abs(eigvals(m)-eigvals(n)) .gt. 0.00000001d0) then
                  Omega_z(m)= Omega_z(m)-aimag(momentum2(m,n,1)*conjg(momentum2(m,n,2)))/(eigvals(m)-eigvals(n))**2
                  Omega_y(m)= Omega_y(m)-aimag(momentum2(m,n,1)*conjg(momentum2(m,n,3)))/(eigvals(m)-eigvals(n))**2
                  Omega_x(m)= Omega_x(m)-aimag(momentum2(m,n,2)*conjg(momentum2(m,n,3)))/(eigvals(m)-eigvals(n))**2
               endif
            enddo !n
         enddo  !m

         Omega_x = Omega_x*2.0
         Omega_y = Omega_y*2.0
         Omega_z = Omega_z*2.0   
         
         do step=1,num_steps                                                             
            fermienergy(step) = fermi_min + (fermi_max-fermi_min)*real(step)/real(num_steps)        
            mu = fermienergy(step)

            do m= 1, num_wann
 !              if (-1.0*(eigvals(m)-mu)/kbT > 20) then
 !                  Omega_x_t(m)= Omega_x(m)*(1.0/(temperature))*((eigvals(m)-mu)*fermi(eigvals(m)-mu, Beta_fake)+(kbT*(-(eigvals(m)-mu)/kbT)))
 !                  Omega_y_t(m)= Omega_y(m)*(1.0/(temperature))*((eigvals(m)-mu)*fermi(eigvals(m)-mu, Beta_fake)+(kbT*(-(eigvals(m)-mu)/kbT)))
 !                  Omega_z_t(m)= Omega_z(m)*(1.0/(temperature))*((eigvals(m)-mu)*fermi(eigvals(m)-mu, Beta_fake)+(kbT*(-(eigvals(m)-mu)/kbT)))
 !              else
 !                  Omega_x_t(m)= Omega_x(m)*(1.0/(temperature))*((eigvals(m)-mu)*fermi(eigvals(m)-mu, Beta_fake)+(kbT*log(1+exp(-(eigvals(m)-mu)/kbT))))
 !                  Omega_y_t(m)= Omega_y(m)*(1.0/(temperature))*((eigvals(m)-mu)*fermi(eigvals(m)-mu, Beta_fake)+(kbT*log(1+exp(-(eigvals(m)-mu)/kbT))))
 !                  Omega_z_t(m)= Omega_z(m)*(1.0/(temperature))*((eigvals(m)-mu)*fermi(eigvals(m)-mu, Beta_fake)+(kbT*log(1+exp(-(eigvals(m)-mu)/kbT))))
 !               endif
            
                Omega_x_t(m)= Omega_x(m)*fermi(eigvals(m)-mu, Beta_fake)
                Omega_y_t(m)= Omega_y(m)*fermi(eigvals(m)-mu, Beta_fake)
                Omega_z_t(m)= Omega_z(m)*fermi(eigvals(m)-mu, Beta_fake)
            enddo
               sigma_tensor_ahc_mpi(1, step)= sigma_tensor_ahc_mpi(1, step)+ real(sum(Omega_x_t))
               sigma_tensor_ahc_mpi(2, step)= sigma_tensor_ahc_mpi(2, step)+ real(sum(Omega_y_t))
               sigma_tensor_ahc_mpi(3, step)= sigma_tensor_ahc_mpi(3, step)+ real(sum(Omega_z_t))
           enddo



      Omega_x_shc=0d0;Omega_y_shc=0d0; Omega_z_shc=0d0
         do n=23,num_wann
            do m=1,22
               if (abs(eigvals(m)-eigvals(n)) .gt. 0.00000001d0) then
                  Omega_z_shc(m)= Omega_x_shc(m)-aimag(momentum2(m,n,1)*conjg(spinmomentum2(m,n,2)))/(eigvals(m)-eigvals(n))**2
                  Omega_y_shc(m)= Omega_y_shc(m)-aimag(momentum2(m,n,1)*conjg(spinmomentum2(m,n,3)))/(eigvals(m)-eigvals(n))**2
                  Omega_x_shc(m)= Omega_z_shc(m)-aimag(momentum2(m,n,2)*conjg(spinmomentum2(m,n,3)))/(eigvals(m)-eigvals(n))**2
               endif
            enddo !n
         enddo  !m

         Omega_x_shc = Omega_x_shc*2.0
         Omega_y_shc = Omega_y_shc*2.0
         Omega_z_shc = Omega_z_shc*2.0   
         
         do step=1,num_steps                                                             
            fermienergy(step) = fermi_min + (fermi_max-fermi_min)*real(step)/real(num_steps)        
            mu = fermienergy(step)

            do m= 1, num_wann         
                Omega_x_t_shc(m)= Omega_x_shc(m)!*fermi(eigvals(m)-mu, Beta_fake)
                Omega_y_t_shc(m)= Omega_y_shc(m)!*fermi(eigvals(m)-mu, Beta_fake)
                Omega_z_t_shc(m)= Omega_z_shc(m)!*fermi(eigvals(m)-mu, Beta_fake)
            enddo
               sigma_tensor_shc_mpi(1, step)= sigma_tensor_shc_mpi(1, step)+ real(sum(Omega_x_t_shc))
               sigma_tensor_shc_mpi(2, step)= sigma_tensor_shc_mpi(2, step)+ real(sum(Omega_y_t_shc))
               sigma_tensor_shc_mpi(3, step)= sigma_tensor_shc_mpi(3, step)+ real(sum(Omega_z_t_shc))
           enddo



   enddo

   num_steps_tot =num_steps*3       
   
   ! call MPI_REDUCE(rho,rho_mpi,size(rho),MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr) 
   call MPI_REDUCE(rho, rho_mpi, size(rho), MPI_DOUBLE_COMPLEX, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
   call MPI_REDUCE(sigma_tensor_ahc_mpi,sigma_tensor_ahc_mpi2,num_steps_tot,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)
   call MPI_REDUCE(sigma_tensor_shc_mpi,sigma_tensor_shc_mpi2,num_steps_tot,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm_world,ierr)

   sigma_tensor_ahc_mpi=sigma_tensor_ahc_mpi2
   sigma_tensor_ahc_mpi=sigma_tensor_ahc_mpi/knv3/volume*condq*1.0e8*twopi

   sigma_tensor_shc_mpi=sigma_tensor_shc_mpi2
   sigma_tensor_shc_mpi=sigma_tensor_shc_mpi/knv3/volume*condq*1.0e8*twopi


   conductivity_ahe=sigma_tensor_ahc_mpi(3,:) 
   conductivity13_ahe=sigma_tensor_ahc_mpi(2,:)   
   conductivity23_ahe=sigma_tensor_ahc_mpi(1,:)
    
   conductivity_she=sigma_tensor_shc_mpi(3,:) 
   conductivity13_she=sigma_tensor_shc_mpi(2,:)   
   conductivity23_she=sigma_tensor_shc_mpi(1,:)

   if(irank.eq.0)then 
      ! write(*,*) "rho_mpi = ",rho_mpi

      write(*,*) '矩阵的实部:'
      do i = 1, num_wann
         do j = 1, num_wann
            ! 使用格式化输出，每个元素占10个字符宽，小数点后3位
            write(*, '(F10.3)', advance='no') real(rho_mpi(i, j))
         end do
      write(*,*)  ! 换行
      end do


      call pauli_block_all2(rho_mpi,num_wann,pauli_result)
   endif

   if(irank.eq.0)then 
      open(123,file='output_ahe_condquant',recl=10000) 
         do step=1,num_steps 
            write(123,*)"fermienergy=",fermienergy(step)                                                          
            write(123,*)"occupation=",occupation(step)                                
             write(123,*)"conductivity=",conductivity_ahe(step)/condq*amat(3,3)*1.0e-8  !2D direct convert to integer
             write(123,*)"conductivity13=",conductivity13_ahe(step)/condq*amat(3,3)*1.0e-8                 
             write(123,*)"conductivity23=",conductivity23_ahe(step)/condq*amat(3,3)*1.0e-8
  !           write(123,*)"conductivity=",conductivity_ahe(step)   !unit is same with wannier_tool  S/cm
  !           write(123,*)"conductivity13=",conductivity13_ahe(step)
  !           write(123,*)"conductivity23=",conductivity23_ahe(step)                      
            write(123,*)"************************************"                                                  
            enddo 
         close(123)        
         
         open(321,file='output_she_condquant',recl=10000) 
         do step=1,num_steps 
            write(321,*)"fermienergy=",fermienergy(step)                                                          
            write(321,*)"occupation=",occupation(step)                                
             write(321,*)"conductivity=",conductivity_she(step)/condq*amat(3,3)*1.0e-8  !2D direct convert to integer
             write(321,*)"conductivity13=",conductivity13_she(step)/condq*amat(3,3)*1.0e-8                 
             write(321,*)"conductivity23=",conductivity23_she(step)/condq*amat(3,3)*1.0e-8
  !           write(321,*)"conductivity_shc=",conductivity_she(step)   !unit is same with wannier_tool  S/cm
  !           write(321,*)"conductivity13_shc=",conductivity13_she(step)
  !           write(321,*)"conductivity23_shc=",conductivity23_she(step)                      
            write(321,*)"************************************"                                                  
            enddo 
         close(321)  

         open(458,file='charge_on_each_atom',recl=10000) 
         do m=1,4
             write(458,*)"m= ", m
             write(458,*)"charge=",pauli_result(m)
             write(458,*)"************************************" 
         enddo
         close(458) 


      endif 
      call mpi_barrier(mpi_comm_world,ierr) 
      call MPI_Finalize(ierr) 
                                                                                                                             
end program ahc_shc_anc_hao
     
! function fermi(omega, Beta_fake) result(value)
!    implicit none
!    real(kind(1.0d0)), intent(in) :: omega
!    real(kind(1.0d0)), intent(in) :: Beta_fake
!    real(kind(1.0d0)) :: value 
!    if (Beta_fake*omega .ge. 20d0) then
!       value = 0.0
!    elseif (Beta_fake*omega .le. -20d0)then
!       value = 1.0
!    else
!       value= 1.0/(1.0+exp(Beta_fake*omega))
!    endif
!    return
! end function fermi    


! function fermi_array(omega, Beta_fake) result(value)
!    implicit none
!    real(kind(1.0d0)), intent(in) :: omega(:)  ! omega 改为数组
!    real(kind(1.0d0)), intent(in) :: Beta_fake
!    real(kind(1.0d0)) :: value(size(omega))   ! 返回值也改为数组
!    integer :: i

!    do i = 1, size(omega)
!       if (Beta_fake*omega(i) .ge. 20d0) then
!          value(i) = 0.0
!       elseif (Beta_fake*omega(i) .le. -20d0) then
!          value(i) = 1.0
!       else
!          value(i) = 1.0/(1.0 + exp(Beta_fake*omega(i)))
!       endif
!    end do
!    return
! end function fermi_array


subroutine now(time_now)

   implicit none
   integer   :: time_new(8)
   real      :: time_now
   call Date_and_time(values=time_new)
   time_now= time_new(3)*24*3600+time_new(5)*3600+&
             time_new(6)*60+time_new(7)+time_new(8)/1000d0  
   return
end subroutine now

subroutine mat_mul(nmatdim,A,B,C)

   implicit none
   integer,intent(in) :: nmatdim

   complex(kind(1.0d0)) :: ALPHA
   complex(kind(1.0d0)) :: BETA

   complex(kind(1.0d0)), intent(in)  :: A(nmatdim ,nmatdim)
   complex(kind(1.0d0)), intent(in)  :: B(nmatdim ,nmatdim)
   complex(kind(1.0d0)), intent(out) :: C(nmatdim,nmatdim)

   ALPHA=1.0d0
   BETA=0.0D0

   C(:,:)=(0.0d0,0.0d0)

   call ZGEMM('N','N',nmatdim,nmatdim,nmatdim,ALPHA,A,nmatdim,B,nmatdim,BETA,C,nmatdim)

   return
end subroutine mat_mul

! subroutine pauli_block_all(M,size,pauli_result)
    
!    implicit none
!    integer,  intent(in)  :: size
!    integer :: nwann_2
!    complex(kind(1.0d0)), intent(in)  :: M(size,size)
!    complex(kind(1.0d0)), intent(out) :: pauli_result(4)
!    complex(kind(1.0d0)) :: MI(size/2,size/2)
!    complex(kind(1.0d0)) :: Mx(size/2,size/2)
!    complex(kind(1.0d0)) :: My(size/2,size/2)
!    complex(kind(1.0d0)) :: Mz(size/2,size/2)
!    complex(kind(1.0d0)) :: trace_value
!    nwann_2 = size/2
!    MI = 0.0
!    Mx = 0.0
!    My = 0.0
!    Mz = 0.0

!    MI = (M(1:nwann_2:1, 1:nwann_2:1) + M(nwann_2+1:2*nwann_2:1, nwann_2+1:2*nwann_2:1)) / 2.0
!    Mx = (M(nwann_2+1:2*nwann_2:1, 1:nwann_2:1) + M(1:nwann_2:1, nwann_2+1:2*nwann_2:1)) / 2.0
!    My = (M(nwann_2+1:2*nwann_2:1, 1:nwann_2:1) - M(1:nwann_2:1, nwann_2+1:2*nwann_2:1)) * 0.5 * (0.0, 1.0)
!    Mz = (M(1:nwann_2:1, 1:nwann_2:1) - M(nwann_2+1:2*nwann_2:1, nwann_2+1:2*nwann_2:1)) / 2.0

!   call trace(MI, trace_value)
!   pauli_result(1) = trace_value
!   call trace(Mx, trace_value)
!   pauli_result(2) = trace_value
!   call trace(My, trace_value)
!   pauli_result(3) = trace_value
!   call trace(Mz, trace_value)
!   pauli_result(4) = trace_value

! end subroutine pauli_block_all

! subroutine trace(M, trace_value)
!      implicit none
!      complex(kind(1.0d0)), intent(in) :: M(:,:)
!      complex(kind(1.0d0)), intent(out) :: trace_value
!      integer :: i, n

!      n = size(M, 1)
!      trace_value = 0.0

!     do i = 1, n
!         trace_value = trace_value + M(i, i)
!     end do
! end subroutine trace
