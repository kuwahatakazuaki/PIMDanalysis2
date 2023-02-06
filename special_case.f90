module mod_special_case
  use input_parameter, only: Natom, Nbeads, TNstep, jobtype, r, &
      save_beads, FNameBinary1, atom1, atom2, atom3, atom4
  use calc_histogram1D
  use utility
  implicit none
  private
  integer :: Uinp, Uout, ierr
  public :: special_case, out_plane

contains

  subroutine special_case

    select case(jobtype)
      case(91)
        call out_plane  ! atom2-atom1-atom3 -> atom1-atom4
    end select
  end subroutine special_case

! +++++++++++++++++++++++
! +++ Start out_plane +++
! +++++++++++++++++++++++
  subroutine out_plane
    integer :: Istep, Ifile, i, j, k
    character(len=:), allocatable :: out_hist
    real(8) :: r21(3), r31(3), r41(3), rt(3)
    real(8) :: data_max, data_min, data_ave, data_dev, data_err
    out_hist='hist_out_plane.out'
    !write(out_hist,'("outplane_",a,I0,"-",a,I0,"-"a,I0,"to",a,I0)') &
    !        & trim(atom(atom2)), atom2, trim(atom(atom1)), atom1, &
    !        & trim(atom(atom3)), atom3, trim(atom(atom4)), atom4
    !if ( trim(out_hist) == "0") write(out_hist, '("hist_",a,".out")') trim(out_hist)

    !Istep = 0
    !do Ifile = 1, Nfile
      !do k = Nstart(Ifile), Nstep(Ifile)
    do k = 1, TNstep
      Istep = Istep + 1
      do j = 1, Nbeads
        r21(:) = r(:,atom2,j,Istep) - r(:,atom1,j,Istep)
        r31(:) = r(:,atom3,j,Istep) - r(:,atom1,j,Istep)
        r41(:) = r(:,atom4,j,Istep) - r(:,atom1,j,Istep)

        r21(:) = r21(:) / norm( r21(:) )
        r31(:) = r31(:) / norm( r31(:) )
        r41(:) = r41(:) / norm( r41(:) )

        rt(:) = outer_product(r21(:),r31(:))
        data_beads(j,Istep) = (0.5 *  pi - acos( dot_product(r41(:),rt(:)) ) ) * 180.0 / pi
      end do
      data_step(Istep) = sum(data_beads(:,Istep))/dble(Nbeads)
    end do
    !end do
    data_max = maxval(data_beads)
    data_min = minval(data_beads)
    data_ave = sum(data_beads)/size(data_beads)
    call calc_deviation(data_dev, data_err)

    if ( save_beads .eqv. .True. ) then
      open(newunit=Uout,file=FNameBinary1, form='unformatted', access='stream', status='replace')
        do Istep = 1, TNstep
          do i = 1, Nbeads
            write(Uout) data_beads(i,Istep)
          end do
        end do
      close(Uout)
    end if

    open(newunit=Uout, file='step_out_plane.out', status='replace')
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

    print *, "*****START calculating bond length*****"
    !do Ifile = 1, Nfile
    !  write (Uprint,'("    From file",I0, " : "a)') Ifile, trim(out_hist)
    !end do
    print '("    Maximum bond =", F13.6)', data_max
    print '("    Minimum bond =", F13.6)', data_min
    print '("    Average bond =", F13.6)', data_ave
    print '("    St. deviation=", F13.6)', data_dev
    print '("    St. error    =", F13.6)', data_err
    print *, ""
    call calc_1Dhist(out_hist=out_hist) ! you need "data_beads"

  end subroutine out_plane
! +++++++++++++++++++++
! +++ End out_plane +++
! +++++++++++++++++++++

end module mod_special_case




