subroutine calc_angle
  use input_parameter,  &
      only: atom1, atom2, atom3, TNstep, graph_step, save_beads, &
            FNameBinary1, Nbeads, data_step, data_beads, label
  use calc_histogram1D, only: calc_1Dhist
  use utility,          only: calc_deviation, calc_cumulative, reblock_step
  implicit none
  integer :: i, k, Uout ! i=atom, j=beads, k=step
  character(len=:), allocatable :: out_name
  character(len=128) :: name_hist, name_temp
  real(8) :: data_max, data_min, data_ave, data_dev, data_err
  out_name = 'step_angle.out'

  write(name_temp,'(a,I0,"-",a,I0,"-",a,I0)') trim(label(atom1)),atom1,trim(label(atom2)),atom2,trim(label(atom3)),atom3
  write(name_hist,'("hist_",a,".out")') trim(name_temp)

  call calc_angle_sub(atom1,atom2,atom3)
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

  open(Uout,file=out_name,status='replace')
    write(Uout, '("# Maximum angle = ", F13.6)') data_max
    write(Uout, '("# Minimum angle = ", F13.6)') data_min
    write(Uout, '("# Average angle = ", F13.6)') data_ave
    write(Uout, '("# Deviation     = ", F13.6)') data_dev
    write(Uout, '("# St. erro      = ", F13.6)') data_err
    do k = 1, TNstep
      if (mod(k,graph_step) == 0) then
        write(Uout,'(I7,F10.5)') k, data_step(k)
      end if
    end do
  close(Uout)

  print '("    Maximum angle = ", F13.6)', data_max
  print '("    Minimum angle = ", F13.6)', data_min
  print '("    Average angle = ", F13.6)', data_ave
  print '("    Deviation     = ", F13.6)', data_dev
  print '("    St. error     = ", F13.6)', data_err
  print '("    Data is saved in ",a,/)', '"'//trim(out_name)//'"'
  call calc_1Dhist(out_hist=trim(name_hist))
  call reblock_step()
!  call calc_cumulative(out_cumulative)
end subroutine calc_angle

subroutine calc_angle_sub(atom_temp1,atom_temp2,atom_temp3)
  use input_parameter,only: Nstep, Nbeads, r, data_beads, data_step, TNstep
  use utility,        only: norm, pi
  implicit none
  integer, intent(in) :: atom_temp1, atom_temp2, atom_temp3
  integer :: j, k
  real(8) :: vec(3,2)

  do k = 1, TNstep
    do j = 1, Nbeads
      vec(:,1) = r(:,atom_temp1,j,k) - r(:,atom_temp2,j,k)
      vec(:,2) = r(:,atom_temp3,j,k) - r(:,atom_temp2,j,k)

      vec(:,1) = vec(:,1)/ norm(vec(:,1))
      vec(:,2) = vec(:,2)/ norm(vec(:,2))

      data_beads(j,k) = dacos(dot_product(vec(:,1),vec(:,2)))
      data_beads(j,k) = data_beads(j,k) *180.0/pi

    end do
    data_step(k) = sum(data_beads(:,k))/dble(Nbeads)
  end do
end subroutine calc_angle_sub

