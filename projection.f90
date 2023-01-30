subroutine projection
  use input_parameter,  only: atom, atom_num, out_hist, Usave, Nfile, TNstep, Uprint, save_beads, Nbeads, &
      FNameBinary1, graph_step, Nstart, Nstep
  use calc_parameter,   only: r, data_beads, data_step
  use calc_histogram1D, only: calc_1Dhist
  use utility,          only: calc_deviation, calc_cumulative
  implicit none
  integer :: i, j, k, Ifile ! i=atom, j=beads, k=step, l=hist, Ifile=file
  integer :: step
  character(len=128) :: bond_name, out_bond, out_cumulative
  real(8) :: data_max, data_min, data_ave, data_dev, data_err
  real(8) :: e(3), r2(3), rt(3)
  integer :: atom1, atom2, atom3, atom4



  write(bond_name, '(a,I0,"-",a,I0)')   trim(atom(atom_num(1,1))), atom_num(1,1), trim(atom(atom_num(2,1))), atom_num(2,1)
  write(out_bond, '("bond_",a,".out")') trim(bond_name)
  write(out_cumulative, '("cumu_",a,".out")') trim(bond_name)

  if ( trim(out_hist) == "0") write(out_hist, '("hist_",a,".out")') trim(bond_name)

  step = 0
  do Ifile = 1, Nfile
    atom1 = atom_num(1,Ifile)
    atom2 = atom_num(2,Ifile)
    atom3 = atom_num(3,Ifile)
    atom4 = atom_num(4,Ifile)

    do k = Nstart(Ifile), Nstep(Ifile)
    step = step + 1
    do j = 1, Nbeads
      e(:) = r(:,atom4,j,step) - r(:,atom3,j,step)
      e(:) = e(:) / dsqrt( dot_product(e(:),e(:)))

      r2(:) = r(:,atom2,j,step) - r(:,atom1,j,step)
      rt(:) = r2(:) - dot_product(r2(:),e(:)) * e(:)

      data_beads(j,step) = dsqrt( dot_product(rt(:),rt(:)) )
    end do
    data_step(step) = sum(data_beads(:,step)) / dble(Nbeads)
    end do
  end do

!  step = 0
!  do Ifile = 1, Nfile
!    call calc_bond_sub(Ifile,atom_num(1,Ifile),atom_num(2,Ifile),step)
!  end do
  data_max = maxval(data_beads)
  data_min = minval(data_beads)
  data_ave = sum(data_beads)/size(data_beads)
  call calc_deviation(data_dev, data_err)

block
  integer :: Ounit
  if ( save_beads .eqv. .True. ) then
    open(newunit=Ounit,file=FNameBinary1, form='unformatted', access='stream', status='replace')
      do step = 1, TNstep
        do i = 1, Nbeads
          write(Ounit) data_beads(i,step)
        end do
      end do
    close(Ounit)
  end if
end block

  open(Usave, file=out_bond, status='replace')
    do Ifile = 1, Nfile
      write (Usave,'(" # From file",I0, " : "a)') Ifile,trim(bond_name)
    end do
    write(Usave,'(a,F13.6)') " # Maximum bond  = ", data_max
    write(Usave,'(a,F13.6)') " # Minimum bond  = ", data_min
    write(Usave,'(a,F13.6)') " # Average bond  = ", data_ave
    write(Usave,'(a,F13.6)') " # St. deviation = ", data_dev
    write(Usave,'(a,F13.6)') " # St. error     = ", data_err
    do k = 1, TNstep
      if (mod(k,graph_step) == 0) then
        write(Usave,'(I7,F10.5)') k, data_step(k)
      end if
    end do
  close(Usave)

  write(Uprint,*) "*****START calculating bond length*****"
  do Ifile = 1, Nfile
    write (Uprint,'("    From file",I0, " : "a)') Ifile, trim(bond_name)
  end do
  write(Uprint, '("    Maximum bond =", F13.6)') data_max
  write(Uprint, '("    Minimum bond =", F13.6)') data_min
  write(Uprint, '("    Average bond =", F13.6)') data_ave
  write(Uprint, '("    St. deviation=", F13.6)') data_dev
  write(Uprint, '("    St. error    =", F13.6)') data_err
  write(Uprint,*) ""
  call calc_1Dhist(out_hist_ex=out_hist)
!  call calc_cumulative(out_cumulative)
end subroutine projection

