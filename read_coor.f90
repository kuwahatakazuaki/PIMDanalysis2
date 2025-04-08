subroutine read_coor
  use input_parameter, &
      only : Natom, Nbeads, Ncut, Nstep, Nfile, DirResult, FileName, label, Lfirst, &
             r, TNstep, data_step, data_beads
  implicit none
  integer :: i, j, k, check

  print '(a)', " ***** START reading coordinate *****"

  select case(Lfirst)
    case(.True.)
      call read_coor_format()
    case(.False.)
      call read_coor_binary()
  end select

  print '(a,/)', " ***** END reading coordinate *****"

  return
contains

! ++++++++++++++++++++++++++++++++++
! +++++ Start read_coor_binary +++++
! ++++++++++++++++++++++++++++++++++
  subroutine read_coor_binary()
    integer :: Uin, ios
    print '(a)', "  *** Binary reading ***"

    if ( trim(FileName(1)) == "0" ) then
      FileName(1)=trim(DirResult(1))//'/coor.bin'
    end if

    ! --- Reading binary file ---
    open(newunit=Uin, file=trim(FileName(1)), form='unformatted', access='stream', status='old',iostat=ios)
      if ( ios /= 0) then
        print '(a,a)', 'ERROR!!: There is no coordinate of ', FileName; stop
      end if

      read(Uin) Natom, Nbeads, Ncut, Nstep, TNstep
!print *, Natom, Nbeads, Ncut, Nstep, TNstep
      allocate(r(3,Natom,Nbeads,TNstep))
      if ( allocated(label) .eqv. .False.) then
        allocate(label(Natom))
      end if

      read(Uin) label(:)

      do k = 1, TNstep
        do j = 1, Nbeads
          do i = 1, Natom
            read(Uin,end=911) r(:,i,j,k)
          end do
        end do
      end do
    close(Uin)
    ! --- Reading binary file ---

    print '("   Binary file reading ")'
    print '("     FileName = ",a)',  trim(FileName(1))
    print '("     TNstep  = ",I0)', TNstep
    print '("     Nbeads  = ",I0)', Nbeads

    return
    911 print *, "ERROR!!: Reading line exceed the coor lines"; stop
  end subroutine read_coor_binary
! ++++++++++++++++++++++++++++++++++
! +++++ End!! read_coor_binary +++++
! ++++++++++++++++++++++++++++++++++


! ++++++++++++++++++++++++++++++++++
! +++++ Start read_coor_format +++++
! ++++++++++++++++++++++++++++++++++
  subroutine read_coor_format()
    integer :: Uin, ios, Istep, Ifile
    logical :: L0step = .False.

    print '(a)', "  *** Format reading ***"

    ! +++ First reading for Automatically counting Nstep +++
    do Ifile = 1, Nfile
      if ( FileName(Ifile) == "0" ) then
        FileName(Ifile)=trim(DirResult(Ifile))//'/coor.xyz'
      end if

      open(newunit=Uin, file=FileName(Ifile), status='old', iostat=ios)
        if ( ios /= 0 ) then
          print '(a,a)', 'ERROR!!: There is no coordinate of ', FileName(Ifile); stop
        end if

        ! --- Automatically counting Nstep ---
        if ( Nstep(Ifile) == 0 ) then
          L0step = .True.
          print '(a)', '  Nstep is set to be 0'
          do
            read(Uin,'()',end=100)
            read(Uin,'()',end=100)
            do j = 1, Nbeads
              do i = 1, Natom
                read(Uin,'()',end=100)
              end do
            end do
            Nstep(Ifile) = Nstep(Ifile) + 1
          end do
          100 continue

          print '(a,I0)', '  Reading Nstep automatically as ',Nstep(Ifile)
          rewind(Uin)
        end if

        ! --- End Automatically counting Nstep ---

      close(Uin)
    end do

    if ( L0step .eqv. .True.) then
      TNstep = sum(Nstep(:)) - sum(Ncut)
      print '(a, i0)', "   The total number of step = ", TNstep
      deallocate(r)
      allocate(r(3,Natom,Nbeads,TNstep))
    end if
    ! +++ End First reading for Automatically counting Nstep +++

    ! +++ Start Second reading for r(3,Natom,Nbeads,Nstep) +++
    Istep = 0
    do Ifile = 1, Nfile

      ! --- Reading formated file ---
      open(newunit=Uin, file=FileName(Ifile), status='old', iostat=ios)
        ! --- Check the Natom and Nbeads ---
        read(Uin,*) check
          if (check /= Natom * Nbeads) then
            print *, "ERROR!!: The number of atom or beads is wrong"
            stop
          end if
        rewind(Uin)
        ! --- End Check the Natom and Nbeads ---

      ! --- skip reading coor file ---
        do k = 1, Ncut(Ifile)
          read(Uin,'()',end=911)
          read(Uin,'()',end=911)
          do j = 1, Nbeads
            do i = 1, Natom
              read(Uin,'()',end=911)
            end do
          end do
        end do
      ! --- End skip reading coor file ---

        do k = Ncut(Ifile)+1, Nstep(Ifile)
          Istep = Istep + 1
          read(Uin, '()',end=911)
          read(Uin, '()',end=911)
          do j = 1, Nbeads
            do i = 1, Natom
              read(Uin,*,end=911) label(i), r(:,i,j,Istep)
            end do
          end do
        end do
      close(Uin)
      ! --- Reading formated file ---
    end do
    ! +++ End!! Second reading for r(3,Natom,Nbeads,Nstep) +++

    print '("   Formatted file reading ")'
    call compression

    return
    911 print *, "ERROR!!: Reading line exceed the coor lines"; stop
  end subroutine read_coor_format
! ++++++++++++++++++++++++++++++++++
! +++++ End!! read_coor_format +++++
! ++++++++++++++++++++++++++++++++++

  subroutine compression()
    character(len=:),allocatable :: ftemp
    integer :: Uout

    print '("  *** START compression ***")'
    ftemp=trim(DirResult(1))//'/coor.bin'
    open(newunit=Uout,file=ftemp, form='unformatted', access='stream', status='replace')
      write(Uout) Natom, Nbeads, Ncut, Nstep, TNstep
      write(Uout) label(:)
!print *, Natom, Nbeads, Ncut, Nstep, TNstep
      do k = 1, TNstep
        do j = 1, Nbeads
          do i = 1, Natom
            write(Uout) r(:,i,j,k)
          end do
        end do
      end do
    close(Uout)

    print '(a)', '    Binary data is saved as "coor.bin"'
    print '(a)', '      in "'//trim(DirResult(1))//'" directory '
    print '("  *** END compression ***")'
    print *, ""
  end subroutine compression

end subroutine read_coor



