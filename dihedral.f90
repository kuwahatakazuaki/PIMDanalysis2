subroutine calc_dihedral
  use input_parameter
  use calc_histogram1D, only: calc_1Dhist
  use utility, only: calc_deviation, calc_cumulative, reblock_step
  implicit none
  integer :: i, k, Uout ! r(xyz,i=atom,j=beads,k=step)
  character(len=:), allocatable :: out_name
  character(len=128) :: name_hist, name_temp
  real(8) :: data_max, data_min, data_ave, data_dev


  write(name_temp,'(a,I0,"-",a,I0,"_"a,I0,"-",a,I0)') &
                  trim(label(atom1)),atom1,trim(label(atom2)),atom2,trim(label(atom3)),atom3,trim(label(atom4)),atom4
  write(name_hist,'("hist_",a,".out")') trim(name_temp)
  out_name = 'step_dihed.out'

  call calc_dihedral_sub(atom1,atom2,atom3,atom4)
  data_max = maxval(data_beads)
  data_min = minval(data_beads)
  data_ave = sum(data_beads)/size(data_beads)
  call calc_deviation(data_dev)

  if ( save_beads .eqv. .True. ) then
    open(newunit=Uout,file=FNameBinary1, form='unformatted', access='stream', status='replace')
      do k = 1, TNstep
        do i = 1, Nbeads
          write(Uout) data_beads(i,k)
        end do
      end do
    close(Uout)
  end if


  open(newunit=Uout,file=out_name,status='replace')
!    write(Uout, '("# Angel Histgram of ",a)') trim(angle_name)
    write(Uout, '("# Maximum angle = ", F13.6)') data_max
    write(Uout, '("# Minimum angle = ", F13.6)') data_min
    write(Uout, '("# Average angle = ", F13.6)') data_ave
    write(Uout, '("# Deviation     = ", F13.6)') data_dev
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
  print '("    Data is saved in ",a,/)', '"'//trim(out_name)//'"'
  call calc_1Dhist(out_hist=trim(name_hist))
  call reblock_step()
!  call calc_cumulative(out_cumulative)
end subroutine calc_dihedral


subroutine calc_dihedral_sub(Tatom1, Tatom2, Tatom3, Tatom4)
  use input_parameter, only: r, data_beads, data_step, TNstep, Nbeads
  use utility, only: norm, cross_product, pi
  implicit none
  integer, intent(in) :: Tatom1, Tatom2, Tatom3, Tatom4
  integer :: j, k ! i=atom, j=beads, k=step
  real(8) :: r12(3), r23(3), r34(3)
  real(8) :: vec123(3), vec234(3)

  do k = 1, TNstep
    do j = 1, Nbeads
      r12(:) = r(:,Tatom1,j,k) - r(:,Tatom2,j,k)
      r23(:) = r(:,Tatom2,j,k) - r(:,Tatom3,j,k)
      r34(:) = r(:,Tatom3,j,k) - r(:,Tatom4,j,k)

      vec123(:) = cross_product(r12(:),r23(:))
      vec234(:) = cross_product(r23(:),r34(:))

      vec123(:) = vec123(:) / norm(vec123(:))
      vec234(:) = vec234(:) / norm(vec234(:))

      data_beads(j,k) = dacos( dot_product(vec123(:), vec234(:)) )
      data_beads(j,k) = data_beads(j,k) *180.0/pi
    end do
    data_step(k) = sum(data_beads(:,k))/dble(Nbeads)
  end do
end subroutine calc_dihedral_sub

