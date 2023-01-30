subroutine read_coor !(Ifile,step)
  use input_parameter, &
      only : Natom, Nbeads, Ncut, Nstep, Nfile, DirResult, FileName, label, Lfirst, &
             r, TNstep, data_step, data_beads
  implicit none
  integer :: i, j, k, check

  select case(Lfirst)
    case(.True.)
      call read_coor_format()
    case(.False.)
      call read_coor_binary()
  end select

  print '(a,/)', " ***** END reading coordinate *****"

  return
contains

  subroutine read_coor_binary()
    integer :: Uin, ios

    if ( trim(FileName(1)) == "0" ) then
      FileName(1)=trim(DirResult(1))//'/coor.bin'
    end if

    ! --- Reading binary file ---
    open(newunit=Uin, file=trim(FileName(1)), form='unformatted', access='stream', status='old',iostat=ios)
      if ( ios /= 0) then
        print '(a,a)', 'ERROR!!: There is no coordinate of ', FileName; stop
      end if

      read(Uin) Natom, Nbeads, Ncut, Nstep, TNstep

      allocate(r(3,Natom,Nbeads,TNstep))
      if ( allocated(label) ) then
      else
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

    return
    911 print *, "ERROR!!: Reading line exceed the coor lines"; stop
  end subroutine read_coor_binary

  subroutine read_coor_format()
    integer :: Uin, ios, Istep, Ifile

    Istep = 0
    do Ifile = 1, Nfile

      if ( FileName(Ifile) == "0" ) then
        FileName(Ifile)=trim(DirResult(Ifile))//'/coor.xyz'
      end if

      ! --- Reading formated file ---
      open(newunit=Uin, file=FileName(Ifile), status='old', iostat=ios)
        if ( ios /= 0 ) then
          print '(a,a)', 'ERROR!!: There is no coordinate of ', FileName(Ifile); stop
        end if

        ! --- Check the Natom and Nbeads ---
        read(Uin,*) check
          if (check /= Natom * Nbeads) then
            print *, "ERROR!!: The number of atom or beads is wrong"
            stop
          end if
        rewind(Uin)
        ! --- End Check the Natom and Nbeads ---

        !! --- Automatically counting Nstep ---
        !if ( Nstep(Ifile) == 0 ) then
        !  print '(a)', '  Nstep is set to be 0'
        !  do
        !    read(Uin,'()',end=100)
        !    read(Uin,'()',end=100)
        !    do j = 1, Nbeads
        !      do i = 1, Natom
        !        read(Uin,'()',end=100)
        !      end do
        !    end do
        !    Nstep(Ifile) = Nstep(Ifile) + 1
        !  end do
        !  100 continue
        !  rewind(Uin)
        !  print '(a,I0)', '  Reading Nstep automatically as ',Nstep(Ifile)
        !end if
        !! --- End Automatically counting Nstep ---

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

    call compression

    return
    911 print *, "ERROR!!: Reading line exceed the coor lines"; stop
  end subroutine read_coor_format

  subroutine compression()
    character(len=:),allocatable :: ftemp
    integer :: Uout

    print '(" ***** START compression *****")'
    ftemp=trim(DirResult(1))//'/coor.bin'
    open(newunit=Uout,file=ftemp, form='unformatted', access='stream', status='replace')
      write(Uout) Natom, Nbeads, Ncut, Nstep, TNstep
      write(Uout) label(:)
      do k = 1, TNstep
        do j = 1, Nbeads
          do i = 1, Natom
            write(Uout) r(:,i,j,k)
          end do
        end do
      end do
    close(Uout)

    print '(a)', '    Binary data is saved in "coor.bin"'
    print '(a)', '    in "'//trim(DirResult(1))//'" directory '
    print '(" ***** END compression *****")'
    print *, ""
  end subroutine compression

end subroutine read_coor



