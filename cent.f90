module calc_centoroid
  use input_parameter
  implicit none
  real(8), public, allocatable :: r_cent(:,:,:) ! r_cent(3,Natom,TNstep)
  real(8), public, allocatable :: data_cent(:)  ! data_cent(TNstep)
  real(8), public :: data_max_cent, data_min_cent, data_ave_cent
  integer, private :: i, j, k, l, Ifile, step ! i=atom, j=beads, k=step, l=hist, Ifile=file
contains

subroutine calc_cent
  use calc_parameter
!  integer :: xyz, i, j, k
  character(len=128) out_cent

  allocate(r_cent(3,Natom,TNstep))

  write(out_cent,'("centroid",".xyz")')
  r_cent(:,:,:) = 0.0d0
  do j = 1, Nbeads
    r_cent(:,:,:) = r_cent(:,:,:) + r(:,:,j,:)
  end do
  r_cent(:,:,:) = r_cent(:,:,:) / Nbeads
  open(20,file=out_cent,status='replace')
  do k = 1, TNstep
    write(20,*) Natom
    write(20,*) k
    do i = 1, Natom
      write(20,*)  atom(i), r_cent(:,i,k)
    end do
  end do
  close(20)

  if (jobtype == 31) then
    call calc_cent_1Dhist
  end if
end subroutine calc_cent

subroutine calc_cent_1Dhist
  use calc_histogram1D
  character(len=128) out_bond
  integer :: atom1(Nfile),atom2(Nfile),atom3(Nfile),atom4(Nfile),atom5(Nfile)
  atom1(:) = atom_num(1,:)
  atom2(:) = atom_num(2,:)
  atom3(:) = atom_num(3,:)
  atom4(:) = atom_num(4,:)
  atom5(:) = atom_num(5,:)

  write(out_hist, '("hist_cent_",a,I0,"-",a,I0".out")') trim(atom(atom1(1))), atom1(1), trim(atom(atom2(1))), atom2(1)
  write(out_bond, '("bond_cent_",a,I0,"-",a,I0".out")') trim(atom(atom1(1))), atom1(1), trim(atom(atom2(1))), atom2(1)
  step = 0
  do Ifile = 1, Nfile
!    call calc_cent_1Dhist_sub_old(Ifile,atom1(Ifile),atom2(Ifile),step)
    call calc_cent_1Dhist_sub(Ifile,atom1(Ifile),atom2(Ifile),step)
  end do
  data_max_cent = maxval(data_cent)
  data_min_cent = minval(data_cent)
  data_ave_cent = sum(data_cent)/size(data_cent)
!  call calc_deviation
!  call calc_1Dhist(out_hist)
  call calc_1Dhist(0.0d0,0.0d0)
end subroutine calc_cent_1Dhist

subroutine calc_cent_1Dhist_sub(Ifile,atomA,atomB,step)
  use utility
  integer, intent(in) :: Ifile, atomA, atomB
  integer, intent(inout) :: step
  integer :: k
!  real(8) :: norm
  do k = Nstart(Ifile), Nstart(Ifile)
  step = step + 1
    data_cent = norm( r_cent(:,atomA,step) - r_cent(:,atomB,step))
  end do
end subroutine calc_cent_1Dhist_sub


!subroutine  calc_cent_1Dhist_sub_old(Ifile,atomA,atomB,step)
!  use input_parameter
!  use calc_parameter
!  integer :: i, j, k ! r(xyz,i=atom,j=beads,k=step)
!  integer, intent(in) :: Ifile, atomA, atomB ! data_beads(j=beads,k=step)
!  integer, intent(inout) :: step
!  real(8) :: norm
!  !do k = 1, TNstep
!  do k = Nstart(Ifile), Nstep(Ifile)
!  step = step+1
!    data_beads(1,step) = norm( r_cent(:,atomA,step) - r_cent(:,atomB,step) )
!  end do
!end subroutine calc_cent_1Dhist_sub_old


end module calc_centoroid


