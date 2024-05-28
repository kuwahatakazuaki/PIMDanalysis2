module calc_histogram1D
  use input_parameter
!  use input_parameter, &
!      only: Natom, Nbeads, TNstep, Nhist, Nfile, Nbond, &
!            hist_min, hist_max, hist_margin
  implicit none
  private
  integer :: j, k, l, Uout
  real(8) :: Dhist, beta
  real(8), allocatable :: histogram(:,:) ! 1:Xaxis, 2:Yaxis
  public :: calc_1Dhist

contains


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++ Start calc_1Dhist                                                            +++
! +++ Need "data_beads(Nbond,TNstep)" and "data_step"                              +++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine calc_1Dhist(hist_min, hist_max,out_hist)
    use utility
    integer :: i
    real(8), intent(in), optional :: hist_min, hist_max
    character(len=*), intent(in), optional :: out_hist
    character(len=:), allocatable :: out_name
    real(8) :: data_max, data_min, data_ave, data_dev, data_err
    real(8) :: hist_umbre(Nhist)


    if (present(hist_min)) hist_min1 = hist_min
    if (present(hist_max)) hist_max1 = hist_max
    if (present(out_hist)) then
      out_name = out_hist
    else
      out_name = 'hist1D.out'
    end if
    allocate(histogram(Nhist,2))

    print '("  *** START calculation 1D histgram ****")'

    data_max = maxval(data_beads)
    data_min = minval(data_beads)
    data_ave = sum(data_beads)/size(data_beads)

    if ( hist_min1 == 0.0d0 .and. hist_max1 == 0.0d0 ) then
      print *, "   Using the margin parameter"
      hist_min1 = data_min - hist_margin
      hist_max1 = data_max + hist_margin
    else
      print *, "   Using the hist_X_min and hist_X_max"
    end if

    Dhist = (hist_max1 - hist_min1) / dble(Nhist)

    print '("    X range max =", F13.6)', hist_max1
    print '("    X range min =", F13.6)', hist_min1
    print '("    Number hist =", I8)',    Nhist
    print '("    Delta hist  =", F13.6)', Dhist

    call calc_1Dhist_sub

    if ( Lfolding .eqv. .True. ) then
      block
      real(8) :: temp(Nhist,2)
      do i = 1, Nhist
        temp(i,2) = histogram(i,2) + histogram(Nhist-i+1,2)
      end do
      histogram(:,2) = temp(:,2) * 0.5d0

      data_ave = 0.0d0
      do i = Nhist/2, Nhist
        data_ave = data_ave + histogram(i,1)*histogram(i,2)
      end do
      end block
    end if

    call calc_deviation(data_dev, data_err)

! ************* Umbrella **************
    if ( umbrella_type > 0 ) then
    block
      real(8) :: poten, normalization
      beta = 1 / ( temperature * KtoAU )
      do l = 1, Nhist
        poten = umbrella_force * histogram(l,1)**2
!        poten = umbrella_force * ( histogram(l,1) * AngtoAU )**2
        hist_umbre(l) = histogram(l,2) * dexp(beta*poten)
      end do
      normalization = sum(hist_umbre(:)) * (histogram(Nhist,1) - histogram(1,1)) / dble(Nhist-1)
      hist_umbre(:) = hist_umbre(:) / normalization
    end block
    end if
! ************* Umbrella **************

    open(newunit=Uout, file=out_name, status='replace')
      write(Uout,'(" # Maximum hist =", F13.6)')  data_max
      write(Uout,'(" # Minimum hist =", F13.6)')  data_min
      write(Uout,'(" # Average hist =", F13.6)')  data_ave
      write(Uout,'(" # St. deviation=", F13.6)')  data_dev
      write(Uout,'(" # St. error    =", F13.6)')  data_err
      write(Uout,'(" # X range max  =", F13.6)')  hist_max1
      write(Uout,'(" # X range min  =", F13.6)')  hist_min1
      write(Uout,'(" # Max of hist  =", F13.6, I3)') maxval(histogram(:,2)), int(maxval(histogram(:,2)))+1 !maxloc(histogram(:,2))
      write(Uout,'(" # Number hist  =", I8)'   )  Nhist
      write(Uout,'(" # Delta hist   =", F13.6)')  Dhist

      if ( umbrella_type == 0 ) then
        do l = 1, Nhist
          write(Uout,'(F13.6,E13.4)') histogram(l,:)
        end do
      else if ( umbrella_type == 1 ) then
        do l = 1, Nhist
          write(Uout,'(F13.6,2E13.4)') histogram(l,:), hist_umbre(l)
        end do
      end if

    close(Uout)
    deallocate(histogram)
    print '(a,a)', "    Hist data is saved in ", '"'//trim(out_name)//'"'
    print '(a,/)', "  *** END calculation 1D histgram ***"
  end subroutine calc_1Dhist
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++ Start calc_1Dhist_sub ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine calc_1Dhist_sub
    integer :: Ihist, UTNstep
    UTNstep = ubound(data_beads,dim=2)
    histogram(:,:) = 0.0d0
    do l = 1, Nhist
      histogram(l,1) = hist_min1 + Dhist*dble(l)
    end do

    do k = 1, UTNstep
      do j = 1, ubound(data_beads, dim=1)
        Ihist = int( (data_beads(j,k)-hist_min1) / Dhist )+1
        histogram(Ihist,2) = histogram(Ihist,2) + 1.0d0
      end do
    end do
    histogram(:,1) = histogram(:,1) - 0.5d0 * Dhist
    histogram(:,2) = histogram(:,2)/(dble(UTNstep*Nbeads)*Dhist)
  end subroutine calc_1Dhist_sub
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine external_1Dhist
    integer :: Nunit

  !block ! Check the number of line between Finput1 and Finput2
    integer :: Nline1, Nline2, Nline
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
  !end block

    901 stop "ERROR!! Tere is no binary input1 file"
    902 stop "ERROR!! Tere is no binary input2 file"
  end subroutine external_1Dhist
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
end module calc_histogram1D

