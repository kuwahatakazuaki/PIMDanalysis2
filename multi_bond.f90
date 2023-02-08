subroutine multi_bond
  use input_parameter,  &
      only: jobtype, Natom, Nbeads, TNstep, label, FNameBinary1, save_beads, &
            Nbond, Imulti, graph_step, &
            data_beads, data_step
  use calc_histogram1D, only: calc_1Dhist
  implicit none
  integer :: i, j, k, Uout
  character(len=128) out_bond, name_hist

  print '(a,/)', " ***** START Multibond calculation ******"
  select case(jobtype)
    case(11)
      call multi_bond_all
    case(12)
      call multi_bond_sort
    case(13)
      call multi_bond_sum
!    case(14:15)
!      call multi_bond_diff
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

  print '(a)', " ***** End Multibond calculation ******"

  return
contains

! +++ jobtype == 13 +++
subroutine multi_bond_sum
  use utility, only: calc_deviation
  integer :: Uout, Iatom1, Iatom2, OldStep, NewStep
  real(8) :: data_max, data_min, data_ave, data_dev
  real(8) :: data_multi(Nbeads,TNstep,Nbond), data_step_multi(TNstep,Nbond)

  out_bond='step_multi.out'
  name_hist='hist_multi.out'

  do i = 1, Nbond
    Iatom1 = Imulti(1,i)
    Iatom2 = Imulti(2,i)
    call calc_bond_sub(Iatom1,Iatom2)
    data_multi(:,:,i) = data_beads(:,:)
    data_step_multi(:,i) = data_step(:)
  end do

  !data_beads(:,:) = 0.0d0
  !do i = 1, Nbond
  !  data_beads(:,:) = data_beads(:,:) + data_multi(:,:,i)
  !end do
  !data_beads(:,:) = data_beads(:,:) / dble(Nbond)
  data_beads(:,:) = sum(data_multi(:,:,:),dim=3) / dble(Nbond)

  OldStep = TNstep
  NewStep = TNstep * Nbond
  deallocate(data_beads)
  allocate(data_beads(Nbeads,NewStep))
  do i = 1, Nbond
    data_beads(:,1+TNstep*(i-1):1+TNstep*i) = data_multi(:,:,i)
  end do

  data_max = maxval(data_beads)
  data_min = minval(data_beads)
  data_ave = sum(data_beads)/size(data_beads)
  call calc_deviation(data_dev)
  call calc_1Dhist(out_hist=trim(name_hist)) ! you need "data_beads"

  open(newunit=Uout, file=out_bond, status='replace')
    write(Uout,*) "# ", trim(out_bond)
    write(Uout,*) "# Maximum bond  = ", data_max
    write(Uout,*) "# Minimum bond  = ", data_min
    write(Uout,*) "# Average bond  = ", data_ave
    write(Uout,*) "# Standard deviation = ", data_dev
    do k = 1, TNstep
      if (mod(k,graph_step) == 0) then
        write(Uout,'(I7,F10.5)') k, data_step(k)
      end if
    end do
  close(Uout)

  !do i = 1, Nbond
  !  data_beads(:,:) = 0.0d0
  !  data_step(:) = 0.0d0
  !  Iatom1 = Imulti(1,i)
  !  Iatom2 = Imulti(2,i)

  !  write(name_hist, '("hist_",a,I0,"-",a,I0".out")') &
  !    trim(label(Iatom1)),Iatom1,trim(label(Iatom2)),Iatom2
  !  write(out_bond, '("bond_",a,I0,"-",a,I0".out")') &
  !    trim(label(Iatom1)),Iatom1,trim(label(Iatom2)),Iatom2
  !  call calc_bond_sub(Iatom1,Iatom2)
  !end do
end subroutine multi_bond_sum
! +++ End jobtype == 13 +++

! +++ jobtype == 11 +++
subroutine multi_bond_all
  use utility, only: calc_deviation
  integer :: Uout, Iatom1, Iatom2
  real(8) :: data_max, data_min, data_ave, data_dev
  do i = 1, Nbond
    data_beads(:,:) = 0.0d0
    data_step(:) = 0.0d0
    Iatom1 = Imulti(1,i)
    Iatom2 = Imulti(2,i)

    write(name_hist, '("hist_",a,I0,"-",a,I0".out")') &
      trim(label(Iatom1)),Iatom1,trim(label(Iatom2)),Iatom2
    write(out_bond, '("bond_",a,I0,"-",a,I0".out")') &
      trim(label(Iatom1)),Iatom1,trim(label(Iatom2)),Iatom2

    !out_bond='multi.out'
    call calc_bond_sub(Iatom1,Iatom2)

    data_max = maxval(data_beads)
    data_min = minval(data_beads)
    data_ave = sum(data_beads)/size(data_beads)
    call calc_deviation(data_dev)

    open(newunit=Uout, file=out_bond, status='replace')
      write(Uout,*) "# ", trim(out_bond)
      write(Uout,*) "# Maximum bond  = ", data_max
      write(Uout,*) "# Minimum bond  = ", data_min
      write(Uout,*) "# Average bond  = ", data_ave
      write(Uout,*) "# Standard deviation = ", data_dev
      do k = 1, TNstep
        if (mod(k,graph_step) == 0) then
          write(Uout,'(I7,F10.5)') k, data_step(k)
        end if
      end do
    close(Uout)
    call calc_1Dhist(hist_min=0.0d0,hist_max=0.0d0,out_hist=trim(name_hist)) ! you need "data_beads"
  end do
end subroutine multi_bond_all
! +++ End jobtype == 11 +++


! +++ jobtype == 12 +++
subroutine multi_bond_sort
  integer :: Iatom1, Iatom2, step
  real(8) :: data_multi(Nbeads,TNstep,Nbond), data_step_multi(TNstep,Nbond)
  real(8) :: temp, temp_beads(Nbeads)
  real(8) :: hist_min_ex, hist_max_ex

  write(*,*) "*****START calculating the data*****"
  do i = 1, Nbond
    Iatom1 = Imulti(1,i)
    Iatom2 = Imulti(2,i)
    call calc_bond_sub(Iatom1,Iatom2)
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
    write(name_hist, '("hist_multisort_",I0,".out")') i
    data_beads(:,:) = data_multi(:,:,i)
    call calc_1Dhist(hist_min_ex, hist_max_ex, out_hist=name_hist)
  end do

end subroutine multi_bond_sort
! +++ End jobtype == 12 +++


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
!      write(name_hist,'("hist_diff.out")')
!    case(14)
!      data_beads(:,:) = data_multi(:,:,1) + data_multi(:,:,2)
!      data_step(:)    = data_step_multi(:,1) + data_step_multi(:,2)
!      write(name_hist,'("hist_sum.out")')
!  end select
!
!  if ( Lfolding .eqv. .True. ) then
!    data_beads(:,:) = abs( data_beads(:,:) )
!  end if
!  call calc_1Dhist(name_hist_ex=name_hist)
!
!  write(out_bond, '("bond_multi_sort.out")')
!  open(22,file=trim(out_bond),status='replace')
!    do i = 1, TNstep
!      write(22,'(I7,10F10.5)') i, data_step(i)
!    end do
!  close(22)
!
!end subroutine multi_bond_diff



end subroutine multi_bond



