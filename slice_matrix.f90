module slice_matrix
    implicit none

contains

    subroutine get_atom_matrix_rho(ori_rho, num_orb_on_atom,orbital_index, rho_on_atom)
        implicit none
        integer, intent(in) :: num_orb_on_atom
        integer, intent(in) :: orbital_index(:)
        complex(kind(1.0d0)), intent(in) :: ori_rho(:,:)
        ! output matrix
        complex(kind(1.0d0)), intent(out) :: rho_on_atom(:,:)
        integer :: i, j, m, n


        ! check the dimension of rho_on_atom
        if (size(rho_on_atom, 1) /= num_orb_on_atom .or. &
            size(rho_on_atom, 2) /= num_orb_on_atom) then
            print *, "error: the size of rho_on_atom is not match with num_orb_on_atom"
            stop
        end if


        rho_on_atom = (0.0d0, 0.0d0)
        do i=1, num_orb_on_atom
            m = orbital_index(i)
            do j=1, num_orb_on_atom    
                n = orbital_index(j)
                rho_on_atom(i,j) = ori_rho(m,n)
            enddo
        enddo

    end subroutine get_atom_matrix_rho
end module slice_matrix