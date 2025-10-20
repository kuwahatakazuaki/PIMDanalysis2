module calc_histogram2D
  use input_parameter, &
      only: Natom, Nbeads, TNstep, Nhist, Nbond, jobtype,  &
            hist_min_inp1, hist_max_inp1, hist_min_inp2, hist_max_inp2, &
            FNameBinary1, FNameBinary2, &
            hist_margin, r, data_beads, delta, &
            atom1, atom2, atom3, atom4, atom5
  implicit none
  private
  real(8) :: hist2D_min(2), hist2D_max(2), Dhist(2)
  real(8) :: hist2D_ave(2), bond_dev2D(2)
  real(8), allocatable :: hist_data(:,:), hist_axis(:,:)
  real(8), allocatable :: hist2D_bead(:,:,:)
  real(8) :: hist_min(2), hist_max(2)
  integer :: i, j, k, l, Ifile, step! i=atom, j=beads, k=step, l=hist, Ifile=file
  integer :: Usave
  character(len=:), allocatable :: out_name

  public :: calc_2Dhist, calc_2Dhist_sub, external_2Dhits_arbitrary, &
            hist2D_bead
contains

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++ Start calc_2Dhist ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine calc_2Dhist()
    character :: axis(2) = ["X","Y"]
    !hist_min(:) = [hist_min1,hist_min2]
    !hist_max(:) = [hist_max1,hist_max2]
    out_name = 'hist2D.out'

    print '(" ***START calculation 2D histgram***")'

!    select case(jobtype)
!      case(21)
!        write(name_axis(1),'(a,I0,"-",a,I0)') trim(atom(atom_num(1,1))),atom_num(1,1), trim(atom(atom_num(2,1))),atom_num(2,1)
!        write(name_axis(2),'(a,I0,"-",a,I0)') trim(atom(atom_num(3,1))),atom_num(3,1), trim(atom(atom_num(4,1))),atom_num(4,1)
!        write(out_hist,'("hist_",a,"and",a,".out")') trim(name_axis(1)), trim(name_axis(2))
!      case(22)
!        write(name_axis(1),'(a,I0,"-",a,I0)') trim(atom(atom_num(1,1))),atom_num(1,1), trim(atom(atom_num(2,1))),atom_num(2,1)
!        write(name_axis(2),'(a,I0,"-",a,I0,"-",a,I0)') &
!          trim(atom(atom_num(3,1))),atom_num(3,1), trim(atom(atom_num(4,1))),atom_num(4,1),trim(atom(atom_num(5,1))),atom_num(5,1)
!        write(out_hist,'("hist_",a,"and",a,".out")') trim(name_axis(1)), trim(name_axis(2))
!      case default
!        stop 'ERROR!!! wrong "Job type" option'
!    end select

    allocate(hist_data(Nhist,Nhist))
    allocate(hist_axis(Nhist,2))

    select case(jobtype)
      case(21:22)
        allocate(hist2D_bead(Nbeads,TNstep,2))
    end select

    !if (jobtype == 21) then ! 2D hist for bond and bond
    select case(jobtype)
      case(21)
        call calc_bond_sub(atom1,atom2)
        hist2D_bead(:,:,1) = data_beads(:,:) ! data_beads(Nbeads,TNstep)
        hist2D_max(1) = maxval(data_beads)
        hist2D_min(1) = minval(data_beads)
        hist2D_ave(1) = sum(data_beads)/size(data_beads)

        call calc_bond_sub(atom3,atom4)
        hist2D_bead(:,:,2) = data_beads(:,:) ! data_beads(Nbeads,TNstep)
        hist2D_max(2) = maxval(data_beads)
        hist2D_min(2) = minval(data_beads)
        hist2D_ave(2) = sum(data_beads)/size(data_beads)

      case(22)
        call calc_bond_sub(atom1,atom2)
        hist2D_bead(:,:,1) = data_beads(:,:) ! data_beads(Nbeads,TNstep)
        hist2D_max(1) = maxval(data_beads)
        hist2D_min(1) = minval(data_beads)
        hist2D_ave(1) = sum(data_beads)/size(data_beads)

        call calc_angle_sub(atom1,atom2,atom3)
        hist2D_bead(:,:,2) = data_beads(:,:) ! data_beads(Nbeads,TNstep)
        hist2D_max(2) = maxval(data_beads)
        hist2D_min(2) = minval(data_beads)
        hist2D_ave(2) = sum(data_beads)/size(data_beads)

      case default
        hist2D_max(1) = maxval(hist2D_bead(:,:,1))
        hist2D_max(2) = maxval(hist2D_bead(:,:,2))
        hist2D_min(1) = minval(hist2D_bead(:,:,1))
        hist2D_min(2) = minval(hist2D_bead(:,:,2))
    end select

    call calc_2Dhist_sub()
    !call calc_2Dhist_sub(hist2D_bead)
! --- START Calculating 2D histgram --- !

    open(Usave, file=out_name, status='replace')
      do i = 1, 2
        write(Usave,'(" # ",a," axis")') axis(i) !,trim(name_axis(i))
        write(Usave,'(" # Maximum hist  ", F13.6)')  hist2D_max(i)
        write(Usave,'(" # Minimum hist  ", F13.6)')  hist2D_min(i)
        write(Usave,'(" # Average hist  ", F13.6)')  hist2D_ave(i)
      end do

      write(Usave,*) "# Hist parameter is as follows"
      do i = 1, 2
        write(Usave, '(" # Range max  ",a,F13.6)') axis(i), hist_max(i)
        write(Usave, '(" # Range min  ",a,F13.6)') axis(i), hist_min(i)
        write(Usave, '(" # Delta hist ",a,F13.6)') axis(i), Dhist(i)
      end do
      write(Usave, '(" # Number hist   ", I6)')     Nhist

      do i = 1, Nhist
        do j = 1, Nhist
          write(Usave,'(2F11.5, E14.5)') hist_axis(i,1), hist_axis(j,2), hist_data(i,j)
        end do
        write(Usave,*) ""
      end do
    close(Usave)

    do i = 1, 2
      print '("    Range max  ",a,F13.6)', axis(i), hist_max(i)
      print '("    Range min  ",a,F13.6)', axis(i), hist_min(i)
      print '("    Delta hist ",a,F13.6)', axis(i),    Dhist(i)
    end do
    print '("    Data is saved in ",a)', '"'//out_name//'"'
    print '(" *****END 2D histgram*****",/)'

    deallocate(hist_data)
    deallocate(hist_axis)
    !deallocate(hist2D_bead)

  end subroutine calc_2Dhist
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++ End calc_2Dhist ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++ Start calc_2Dhist_sub ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine calc_2Dhist_sub()
  !subroutine calc_2Dhist_sub(hist2D_bead)
    integer :: coun(2), UTNstep, UNbead
    !real(8), intent(in) :: hist2D_bead(:,:,:)
    character :: axis(2) = ["X","Y"]
    hist_min(:) = [hist_min_inp1,hist_min_inp2]
    hist_max(:) = [hist_max_inp1,hist_max_inp2]
    UNbead  = ubound(hist2D_bead,dim=1)
    UTNstep = ubound(hist2D_bead,dim=2)

    do i = 1, 2
      if ( hist_min(i) == 0.0d0 .and. hist_max(i) == 0.0d0 ) then
        print '(a,a)', "    Using the margin parameter to ", axis(i)
        hist_min(i) = hist2D_min(i) - hist_margin
        hist_max(i) = hist2D_max(i) + hist_margin
      else
        print '("    Using the hist_min and hist_max")'
      end if

      Dhist(i) = (hist_max(i) - hist_min(i)) / dble(Nhist)
    end do

    do i = 1, 2
      do l = 1, Nhist
        hist_axis(l,i) = hist_min(i) + Dhist(i) * dble(l)
      end do
    end do

! +++++++ hist_data1D ++++++
    hist_data(:,:) = 0.0d0
    !do j = 1, TNstep
    !  do k = 1, Nbeads
    do j = 1, UTNstep
      do k = 1, UNbead
        do i = 1, 2
          coun(i) = int( (hist2D_bead(k,j,i)-hist_min(i))/Dhist(i) ) + 1
          coun(i) = int( (hist2D_bead(k,j,i)-hist_min(i))/Dhist(i) ) + 1
        end do
        hist_data(coun(1),coun(2)) = hist_data(coun(1),coun(2)) + 1.0d0
      end do
    end do

    do i = 1, 2
      hist_axis(:,i) = hist_axis(:,i) - 0.5d0 * Dhist(i)
    end do
    hist_data(:,:) = hist_data(:,:) / (dble(TNstep*Nbeads)*Dhist(1)*Dhist(2))
  end subroutine calc_2Dhist_sub
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++ End calc_2Dhist_sub ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++ Start external_2Dhits_arbitrary  +++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine external_2Dhits_arbitrary
    implicit none
    integer :: Nunit, Nline, coun(2), Ihist1, Ihist2
    real(8), allocatable :: binary_data(:,:)
    character :: axis(2) = ["X","Y"]
    hist_min(:) = [hist_min_inp1,hist_min_inp2]
    hist_max(:) = [hist_max_inp1,hist_max_inp2]

    out_name='hist_2D_external.out'
    block ! Check the number of line between Finput1 and Finput2
      integer :: Nline1, Nline2
      real(8) :: dummyD
      open(newunit=Nunit, file=FNameBinary1, form='unformatted', access='stream', status='old', err=901)
        Nline = 0
        do
          read(Nunit, end=800) dummyD
          Nline = Nline + 1
        end do
        800 continue
      close(Nunit)
      Nline1 = Nline

      open(newunit=Nunit, file=FNameBinary2, form='unformatted', access='stream', status='old', err=902)
        Nline = 0
        do
          read(Nunit, end=801) dummyD
          Nline = Nline + 1
        end do
        801 continue
      close(Nunit)
      Nline2 = Nline

      if ( Nline1 /= Nline2 ) then
        print *, "The number of lines are different between ", FNameBinary1, " and ", FNameBinary2
        stop "ERROR!!"
      end if
    end block
    allocate(binary_data(Nline,2))
    allocate(hist_axis(Nhist,2))
    allocate(hist_data(Nhist,Nhist))


    open(newunit=Nunit, file=FNameBinary1, form='unformatted', access='stream', status='old', err=901)
      do i = 1, Nline
        read(Nunit) binary_data(i,1)
      end do
    close(Nunit)

    open(newunit=Nunit, file=FNameBinary2, form='unformatted', access='stream', status='old', err=901)
      do i = 1, Nline
        read(Nunit) binary_data(i,2)
      end do
    close(Nunit)

    do i = 1, 2
      hist2D_max(i) = maxval(binary_data(:,i))
      hist2D_min(i) = minval(binary_data(:,i))
      hist2D_ave(i) = sum(binary_data(:,i))/dble(Nline)
    end do

! ----- Start call calc_2Dhist_sub -----
    do i = 1, 2
      if ( abs(hist_min(i)) < delta .and. abs(hist_max(i)) < delta ) then
      !if ( hist_min(i) == 0.0d0 .and. hist_max(i) == 0.0d0 ) then
        print '(a,a)', "    Using the margin parameter to ", axis(i)
        hist_min(i) = hist2D_min(i) - hist_margin
        hist_max(i) = hist2D_max(i) + hist_margin
      else
        print '("    Using the hist_min and hist_max")'
      end if

      Dhist(i) = (hist_max(i) - hist_min(i)) / dble(Nhist)
    end do
    do i = 1, 2
      do l = 1, Nhist
        hist_axis(l,i) = hist_min(i) + Dhist(i) * dble(l)
      end do
    end do

    hist_data(:,:) = 0.0d0
    do k = 1, Nline
      do i = 1, 2
        coun(i) = int( (binary_data(k,i)-hist_min(i))/Dhist(i) ) + 1
        coun(i) = int( (binary_data(k,i)-hist_min(i))/Dhist(i) ) + 1
      end do
      hist_data(coun(1),coun(2)) = hist_data(coun(1),coun(2)) + 1.0d0
    end do

    do i = 1, 2
      hist_axis(:,i) = hist_axis(:,i) - 0.5d0 * Dhist(i)
    end do
    hist_data(:,:) = hist_data(:,:) / (dble(Nline)*Dhist(1)*Dhist(2))
! ----- End call calc_2Dhist_sub -----

    open(Usave, file='hist_2D_external.out', status='replace')
      do i = 1, 2
        write(Usave,'(" # Maximum hist  ", F13.6)')  hist2D_max(i)
        write(Usave,'(" # Minimum hist  ", F13.6)')  hist2D_min(i)
        write(Usave,'(" # Average hist  ", F13.6)')  hist2D_ave(i)
      end do

      write(Usave,*) "# Hist parameter is as follows"
      do i = 1, 2
        write(Usave, '(" # Range max  ",a,F13.6)') axis(i), hist_max(i)
        write(Usave, '(" # Range min  ",a,F13.6)') axis(i), hist_min(i)
        write(Usave, '(" # Delta hist ",a,F13.6)') axis(i), Dhist(i)
      end do
      write(Usave, '(" # Number hist   ", I6)')     Nhist

      do i = 1, Nhist
        do j = 1, Nhist
          write(Usave,'(2F11.5, E14.5)') hist_axis(i,1), hist_axis(j,2), hist_data(i,j)
        end do
        write(Usave,*) ""
      end do
    close(Usave)

    do i = 1, 2
      print '("    Range max  ",a,F13.6)', axis(i), hist_max(i)
      print '("    Range min  ",a,F13.6)', axis(i), hist_min(i)
      print '("    Delta hist ",a,F13.6)', axis(i),    Dhist(i)
    end do
    print '("    Data is saved in ",a)', '"'//out_name//'"'
    print '(" *****END 2D histgram*****",/)'

  return
  901 stop "ERROR!! There is no binary input1 file"
  902 stop "ERROR!! There is no binary input2 file"
  end subroutine external_2Dhits_arbitrary
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++ End external_2Dhits_arbitrary  +++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
end module calc_histogram2D

