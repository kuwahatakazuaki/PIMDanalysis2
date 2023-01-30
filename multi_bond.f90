subroutine multi_bond
  use input_parameter
  use calc_histogram1D, only: calc_1Dhist
  implicit none
  integer :: i, j, k, Uout
  character(len=128) out_bond, out_hist
  real(8) :: hist_min_ex, hist_max_ex

!  hist_min_ex = hist_min(1)
!  hist_max_ex = hist_max(1)

  select case(jobtype)
    case(11)
      call multi_bond_all
    case(12)
      call multi_bond_sort
    case(13:14)
      call multi_bond_diff
    case default
      stop 'ERROR!!! wrong "Job type" option'
  end select

  if ( save_beads .eqv. .True. ) then
    open(newunit=Uout,file=FNameBinary1, form='unformatted', access='stream', status='replace')
      do k = 1, TNstep
        do i = 1, Nbeads
          write(Uout) data_beads(i,k)
        end do
      end do
    close(Uout)
  end if


  return
contains


!! +++ jobtype == 13 +++
!subroutine multi_bond_diff
!  integer :: Ifile
!  real(8) :: data_multi(Nbeads,TNstep,2), data_step_multi(TNstep,2)
!
!  write(*,*) "*****START calculating the data*****"
!
!  step = 0
!  do Ifile = 1, Nfile
!    call calc_bond_sub(Ifile,atom_num(1,Ifile),atom_num(2,Ifile),step)
!  end do
!  data_multi(:,:,1) = data_beads(:,:)
!  data_step_multi(:,1) = data_step(:)
!
!  step = 0
!  do Ifile = 1, Nfile
!    call calc_bond_sub(Ifile,atom_num(3,Ifile),atom_num(4,Ifile),step)
!  end do
!  data_multi(:,:,2) = data_beads(:,:)
!  data_step_multi(:,2) = data_step(:)
!
!  select case(jobtype)
!    case(13)
!      data_beads(:,:) = data_multi(:,:,1) - data_multi(:,:,2)
!      data_step(:)    = data_step_multi(:,1) - data_step_multi(:,2)
!      write(out_hist,'("hist_diff.out")')
!    case(14)
!      data_beads(:,:) = data_multi(:,:,1) + data_multi(:,:,2)
!      data_step(:)    = data_step_multi(:,1) + data_step_multi(:,2)
!      write(out_hist,'("hist_sum.out")')
!  end select
!
!  if ( Lfolding .eqv. .True. ) then
!    data_beads(:,:) = abs( data_beads(:,:) )
!  end if
!  call calc_1Dhist(out_hist_ex=out_hist)
!
!  write(out_bond, '("bond_multi_sort.out")')
!  open(22,file=trim(out_bond),status='replace')
!    do i = 1, TNstep
!      write(22,'(I7,10F10.5)') i, data_step(i)
!    end do
!  close(22)
!
!end subroutine multi_bond_diff


! +++ jobtype == 11 +++
subroutine multi_bond_all
  use utility
  integer :: Uout
  real(8) :: data_max, data_min, data_ave, data_dev
  do i = 1, Nbond
    data_beads(:,:) = 0.0d0
    data_step(:) = 0.0d0

    write(out_hist, '("hist_",a,I0,"-",a,I0".out")') &
      & trim(atom(atom_multi(1,i))), atom_multi(1,i), trim(atom(atom_multi(2,i))), atom_multi(2,i)
    write(out_bond, '("bond_",a,I0,"-",a,I0".out")') &
      & trim(atom(atom_multi(1,i))), atom_multi(1,i), trim(atom(atom_multi(2,i))), atom_multi(2,i)
!    step=0
    out_bond='multi.out'
    call calc_bond_sub(atom_multi(1,i),atom_multi(2,i))

    data_max = maxval(data_beads)
    data_min = minval(data_beads)
    data_ave = sum(data_beads)/size(data_beads)
    call calc_deviation(data_dev)

    open(22, file=out_bond, status='replace')
      write(22,*) "# ", trim(out_bond)
      write(22,*) "# Maximum bond  = ", data_max
      write(22,*) "# Minimum bond  = ", data_min
      write(22,*) "# Average bond  = ", data_ave
      write(22,*) "# Standard deviation = ", data_dev
      do k = 1, TNstep
        if (mod(k,10) == 0) then
          write(22,'(I7,F10.5)') k, data_step(k)
        end if
      end do
    close(22)
    call calc_1Dhist(0.0d0,0.0d0) ! you need "data_beads"
  end do
end subroutine multi_bond_all
! +++ End jobtype == 11 +++


! +++ jobtype == 12 +++
subroutine multi_bond_sort
  real(8) :: data_multi(Nbeads,TNstep,Nbond), data_step_multi(TNstep,Nbond)
  real(8) :: temp, temp_beads(Nbeads)

  write(*,*) "*****START calculating the data*****"
  do i = 1, Nbond
    step = 0
    call calc_bond_sub(1,atom_multi(1,i),atom_multi(2,i),step)
    data_multi(:,:,i) = data_beads(:,:)
    data_step_multi(:,i) = data_step(:)
  end do

  write(*,*) "*****START sorting*****"
  do step = 1, TNstep
    do i = 1, Nbond-1
      do  j = i+1, Nbond
        if ( data_step_multi(step,i) > data_step_multi(step,j)) then
          temp = data_step_multi(step,i)
          data_step_multi(step,i) = data_step_multi(step,j)
          data_step_multi(step,j) = temp

          temp_beads(:) = data_multi(:,step,i)
          data_multi(:,step,i) = data_multi(:,step,j)
          data_multi(:,step,j) = temp_beads(:)
        end if
      end do
    end do
  end do

  write(out_bond, '("bond_multi_sort.out")')
  open(22,file=trim(out_bond),status='replace')
    do i = 1, TNstep
      write(22,'(I7,10F10.5)') i, data_step_multi(i,:)
    end do
  close(22)

  do i = 1, Nbond
    write(out_hist, '("hist_multisort_",I0,".out")') i
    data_beads(:,:) = data_multi(:,:,i)
    call calc_1Dhist(hist_min_ex, hist_max_ex, out_hist_ex=out_hist)
  end do

end subroutine multi_bond_sort
! +++ End jobtype == 12 +++




end subroutine multi_bond



