module pauli
    implicit none
    ! 定义双精度浮点数类型
    integer, parameter :: dp = selected_real_kind(15, 307)

contains

    subroutine pauli_block_all(M, size, pauli_result)
        implicit none
        integer, intent(in) :: size
        complex(kind(1.0d0)), intent(in) :: M(size, size)
        complex(kind(1.0d0)), intent(out) :: pauli_result(4)
        integer :: nwann_2
        complex(kind(1.0d0)) :: MI(size/2, size/2)
        complex(kind(1.0d0)) :: Mx(size/2, size/2)
        complex(kind(1.0d0)) :: My(size/2, size/2)
        complex(kind(1.0d0)) :: Mz(size/2, size/2)
        complex(kind(1.0d0)) :: trace_value

        nwann_2 = size / 2

        MI = (M(1:nwann_2, 1:nwann_2) + M(nwann_2+1:2*nwann_2, nwann_2+1:2*nwann_2)) / 2.0_dp
        Mx = (M(nwann_2+1:2*nwann_2, 1:nwann_2) + M(1:nwann_2, nwann_2+1:2*nwann_2)) / 2.0_dp
        My = (M(nwann_2+1:2*nwann_2, 1:nwann_2) - M(1:nwann_2, nwann_2+1:2*nwann_2)) * 0.5_dp
        Mz = (M(1:nwann_2, 1:nwann_2) - M(nwann_2+1:2*nwann_2, nwann_2+1:2*nwann_2)) / 2.0_dp

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
        complex(kind(1.0d0)), intent(in) :: M(:,:)
        complex(kind(1.0d0)), intent(out) :: trace_value
        integer :: i, n

        n = size(M, 1)
        trace_value = 0.0_dp

        do i = 1, n
            trace_value = trace_value + M(i, i)
        end do
    end subroutine trace

end module pauli
