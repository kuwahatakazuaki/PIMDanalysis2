
subroutine pbhpo4
  use input_parameter, &
        only: jobtype, Natom, Nbeads, TNstep, Nunit, &
              save_beads, FNameBinary1, lattice, &
              r
!  use calc_parameter
  use calc_histogram1D
  use utility
  implicit none
  integer :: i, j, k, atom1, atom2
  real(8) :: lat_inv(3,3)
  real(8), allocatable :: s(:,:,:,:), unit_bead(:,:,:), unit_step(:,:)
  real(8) :: data_max, data_min, data_ave, data_dev, data_err

  call get_inv_mat(lattice,lat_inv,3)

  allocate(s(3,Natom,Nbeads,TNstep))
  allocate(unit_bead(Nunit*2,Nbeads,TNstep))
  allocate(unit_step(Nunit*2,TNstep))

  deallocate(data_beads,data_step)
  allocate(data_beads(Nbeads,TNstep*Nunit*2))
  allocate(data_step(TNstep*Nunit*2))

  do k = 1, TNstep
    do j = 1, Nbeads
      do i = 1, Natom
        s(:,i,j,k) = matmul(r(:,i,j,k), lat_inv(:,:))
      end do
    end do
  end do

  select case(jobtype)
    case(191)
      call oo_bond
    case(192)
      call deltaOH
  end select

contains

  subroutine deltaOH
    do i = 0,2
      call calc_unit_diff(1+i,19+i,1+i,22+i,i+1)
    end do
    do i = 0,1
      call calc_unit_diff(5+i,28+i,5+i,26+i,i+4)
    end do
    call calc_unit_diff(4,30,4,25,6)

    do i = 1, Nunit*2
      data_beads(:,TNstep*(i-1)+1:TNstep*i) = unit_bead(i,:,:)
    end do

    call print_unit_step

!    out_hist="dletaOH.out"
    call calc_1Dhist(out_hist="dletaOH.out")
  end subroutine deltaOH

  subroutine oo_bond
    do i = 0,2
      call calc_unit_bond(19+i,22+i,i+1)
    end do
    do i = 0,1
      call calc_unit_bond(26+i,28+i,i+4)
    end do
    call calc_unit_bond(25,30,6)

    do i = 1, Nunit*2
      data_beads(:,TNstep*(i-1)+1:TNstep*i) = unit_bead(i,:,:)
    end do

    call print_unit_step

    call calc_1Dhist(out_hist="OObond.out")
  end subroutine oo_bond

! +++++++++++++++++++++++++++++++++++++++++++++++++
! +++ bond difference atom1_atom2 - atom3_atom4 +++
! +++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine calc_unit_diff(atom1,atom2,atom3,atom4,Iunit)
    integer, intent(in) :: atom1, atom2, atom3, atom4, Iunit
    real(8) :: r12(3), s12(3), d12, s34(3), r34(3), d34
    integer :: i,j,k
    do k = 1, TNstep
      do j = 1, Nbeads
        s12(:) = s(:,atom1,j,k) - s(:,atom2,j,k)
        s12(:) = s12(:) - nint(s12(:))
        r12(:) = matmul(s12(:),lattice(:,:))
        d12 = dsqrt(sum(r12(:)*r12(:)))

        s34(:) = s(:,atom3,j,k) - s(:,atom4,j,k)
        s34(:) = s34(:) - nint(s34(:))
        r34(:) = matmul(s34(:),lattice(:,:))
        d34 = dsqrt(sum(r34(:)*r34(:)))
        unit_bead(Iunit,j,k) = d12 - d34
      end do
    end do
  end subroutine calc_unit_diff

! +++++++++++++++++++++++++++++++++++++++++
! +++ bond length between atom1 - atom2 +++
! +++++++++++++++++++++++++++++++++++++++++
  subroutine calc_unit_bond(atom1,atom2,Iunit)
    integer, intent(in) :: atom1, atom2, Iunit
    real(8) :: r12(3), s12(3), d12
    integer :: i,j,k
    do k = 1, TNstep
      do j = 1, Nbeads
        s12(:) = s(:,atom1,j,k) - s(:,atom2,j,k)
        s12(:) = s12(:) - nint(s12(:))
        r12(:) = matmul(s12(:),lattice(:,:))
        d12 = dsqrt(sum(r12(:)*r12(:)))
        unit_bead(Iunit,j,k) = d12
      end do
    end do
  end subroutine calc_unit_bond

  subroutine print_unit_step
    integer :: Uout
    do k = 1, TNstep
      do i = 1,Nunit*2
        unit_step(i,k) = sum(unit_bead(i,:,k))/dble(Nbeads)
      end do
    end do

    open(newunit=Uout,file='step_change.out')
      do k = 1, TNstep
        if ( mod(k,graph_step)==0 ) then
          write(Uout,9998) k, unit_step(:,k)
        end if
      end do
    close(Uout)

    9999  format(6F12.7)
    9998  format(i0,6F11.6)
  end subroutine print_unit_step

end subroutine pbhpo4



