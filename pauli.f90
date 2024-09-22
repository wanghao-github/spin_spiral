module pauli

!*************************************************                      
!write by HaoWang September 22 2024
!*************************************************                      
      implicit none 



      subroutine pauli_block_all(M,size,pauli_result)
    
        implicit none
        integer,  intent(in)  :: size
        integer :: nwann_2
        real(kind(1.0d0)), intent(in)  :: M(size,size)
        real(kind(1.0d0)), intent(out) :: pauli_result(4)
        real(kind(1.0d0)) :: MI(size/2,size/2)
        real(kind(1.0d0)) :: Mx(size/2,size/2)
        real(kind(1.0d0)) :: My(size/2,size/2)
        real(kind(1.0d0)) :: Mz(size/2,size/2)
        real(kind(1.0d0)) :: trace_value
        nwann_2 = size/2
        MI = 0.0
        Mx = 0.0
        My = 0.0
        Mz = 0.0
     
        MI = (M(1:nwann_2:1, 1:nwann_2:1) + M(nwann_2+1:2*nwann_2:1, nwann_2+1:2*nwann_2:1)) / 2.0
        Mx = (M(nwann_2+1:2*nwann_2:1, 1:nwann_2:1) + M(1:nwann_2:1, nwann_2+1:2*nwann_2:1)) / 2.0
        My = (M(nwann_2+1:2*nwann_2:1, 1:nwann_2:1) - M(1:nwann_2:1, nwann_2+1:2*nwann_2:1)) * 0.5 * (0.0, 1.0)
        Mz = (M(1:nwann_2:1, 1:nwann_2:1) - M(nwann_2+1:2*nwann_2:1, nwann_2+1:2*nwann_2:1)) / 2.0
     
       call trace(MI, trace_value)
       pauli_result(1) = trace_value
       call trace(Mx, trace_value)
       pauli_result(2) = trace_value
       call trace(My, trace_value)
       pauli_result(3) = trace_value
       call trace(Mz, trace_value)
       pauli_result(4) = trace_value
     
     end subroutine pauli_block_all
     
     subroutine trace(M, trace_value)
          implicit none
          real(kind(1.0d0)), intent(in) :: M(:,:)
          real(kind(1.0d0)), intent(out) :: trace_value
          integer :: i, n
     
          n = size(M, 1)
          trace_value = 0.0
     
         do i = 1, n
             trace_value = trace_value + M(i, i)
         end do
     end subroutine trace
     