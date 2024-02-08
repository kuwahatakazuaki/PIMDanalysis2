
subroutine calc_bond
  use input_parameter,  &
      only: jobtype, Natom, Nbeads, TNstep, &
            atom1, atom2, atom3, atom4, label,  graph_step, &
            r, data_beads, data_step, save_beads, FNameBinary1
  use calc_histogram1D, only: calc_1Dhist
  use utility,          only: calc_deviation, calc_cumulative, reblock_step
  implicit none
  integer :: i, k, Ifile, Uout
  character(len=128) :: name_hist, name_temp
  real(8) :: data_max, data_min, data_ave, data_dev, data_err

  write(name_temp,'(a,I0,"-",a,I0)')    trim(label(atom1)),atom1,trim(label(atom2)),atom2
  write(name_hist,'("hist_",a,".out")') trim(name_temp)

  select case(jobtype)
    case(1)
      call calc_bond_sub(atom1,atom2)
    case(4)
      call calc_bond_diff(atom1,atom2,atom3,atom4)
    case(5)
      call calc_bond_add(atom1,atom2,atom3,atom4)
  end select


  data_max = maxval(data_beads)
  data_min = minval(data_beads)
  data_ave = sum(data_beads)/size(data_beads)
  call calc_deviation(data_dev, data_err)

  if ( save_beads .eqv. .True. ) then
    open(newunit=Uout,file=FNameBinary1, form='unformatted', access='stream', status='replace')
      do k = 1, TNstep
        do i = 1, Nbeads
          write(Uout) data_beads(i,k)
        end do
      end do
    close(Uout)
  end if

  open(newunit=Uout, file='step_bond.out', status='replace')
    write(Uout,'(a,F13.6)') " # Maximum bond  = ", data_max
    write(Uout,'(a,F13.6)') " # Minimum bond  = ", data_min
    write(Uout,'(a,F13.6)') " # Average bond  = ", data_ave
    write(Uout,'(a,F13.6)') " # St. deviation = ", data_dev
    write(Uout,'(a,F13.6)') " # St. error     = ", data_err
    do k = 1, TNstep
      if (mod(k,graph_step) == 0) then
        write(Uout,'(I7,F10.5)') k, data_step(k)
      end if
    end do
  close(Uout)

  print *, "***** START calculating bond length *****"
  print '("    Maximum bond =", F13.6)', data_max
  print '("    Minimum bond =", F13.6)', data_min
  print '("    Average bond =", F13.6)', data_ave
  print '("    St. deviation=", F13.6)', data_dev
  print '("    St. error    =", F13.6)', data_err
  print *, ""
  call calc_1Dhist(out_hist=trim(name_hist))
!  call calc_cumulative(out_cumulative)
  call reblock_step()
end subroutine calc_bond

subroutine calc_bond_sub(atom1,atom2)
  use input_parameter, only: r, data_beads, data_step, TNstep, Nbeads
  implicit none
  integer :: j, k
  integer, intent(in) :: atom1, atom2 ! data_beads(j=beads,k=step)
  real(8) :: r12(3)

  do k = 1, TNstep
    do j = 1, Nbeads
      r12(:) = r(:,atom1,j,k)-r(:,atom2,j,k)
      data_beads(j,k) = dsqrt( sum(r12(:)*r12(:)) )
    end do
    data_step(k) = sum(data_beads(:,k))/dble(Nbeads)
  end do
end subroutine calc_bond_sub

subroutine calc_bond_diff(atom1,atom2,atom3,atom4)
  use input_parameter, only: r, data_beads, data_step, TNstep, Nbeads
  integer, intent(in) :: atom1, atom2, atom3, atom4 ! data_beads(j=beads,k=step)
  integer :: j, k
  real(8) :: r12(3), r34(3), d12, d34
  do k = 1, TNstep
    do j = 1, Nbeads
      r12(:) = r(:,atom1,j,k)-r(:,atom2,j,k)
      r34(:) = r(:,atom3,j,k)-r(:,atom4,j,k)
      d12 = norm2(r12(:))
      d34 = norm2(r34(:))
      data_beads(j,k) = d12 - d34
    end do
    data_step(k) = sum(data_beads(:,k))/dble(Nbeads)
  end do
end subroutine calc_bond_diff

subroutine calc_bond_add(atom1,atom2,atom3,atom4)
  use input_parameter, only: r, data_beads, data_step, TNstep, Nbeads
  integer, intent(in) :: atom1, atom2, atom3, atom4 ! data_beads(j=beads,k=step)
  stop 'Not Updated'
end subroutine calc_bond_add


