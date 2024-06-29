module mod_other_quantities
  use input_parameter, &
      only: jobtype, Natom, Nbeads, TNstep, Ncut, Nfile, Nstep, &
            DirResult, save_beads, FNameBinary1, label, &
            atom1, atom2, graph_step, &
            data_beads, data_step
  use calc_histogram1D
  use utility, only: reblock_step
  implicit none
  private
  real(8), allocatable :: charge(:,:,:), dipole(:,:,:), hfcc(:,:,:)
  real(8), allocatable :: force(:,:,:,:)
  character(:), allocatable :: PathFile
  integer :: Uinp, Uout, ierr
  public :: other_quantities

contains

  subroutine  other_quantities

    select case(jobtype)
      case(61:62)
        call charge_analysis
      case(63)
        call dipole_analysis
      case(64)
        call hfcc_analysis
      case(65)
        call force_analysis
    end select
  end subroutine other_quantities

! ++++++++++++++++++++++
! +++ force_analysis +++
! ++++++++++++++++++++++
  subroutine force_analysis
    integer :: Istep, i, j, k, Ifile
    integer :: Iloc(2)
    character(len=:), allocatable :: Finp, Fout
    real(8) :: force_norm(Natom,Nbeads)

    Fout='force_max.out'
    allocate(force(3,Natom,Nbeads,TNstep))

    open(Uout,file=Fout,status='replace')
    Istep = 0
    do Ifile = 1, Nfile
      Finp=trim(DirResult(Ifile))//'/force.dat'

      open(newunit=Uinp, file=Finp, status='old', iostat=ierr)
        if ( ierr > 0 ) then
          print *, 'Check the path : ', Finp
          stop 'ERROR!!: There is no "force.dat"'
        end if

        do i = 1, Ncut(Ifile)
          read(Uinp,'()')
          do j = 1, Nbeads
            read(Uinp,'()')
          end do
        end do

        do k = Ncut(Ifile)+1, Nstep(Ifile)
          Istep = Istep + 1
          read(Uinp,'()')
          do j = 1, Nbeads
            do i = 1, Natom
              read(Uinp,*) force(:,i,j,Istep)
            end do
          end do
          force_norm(:,:) = force(1,:,:,Istep)**2 + &
                            force(2,:,:,Istep)**2 + &
                            force(3,:,:,Istep)**2
          Iloc(:) = maxloc(force_norm)
          !write(Uout,*) Istep, dsqrt(maxval(force_norm)), trim(label(Iloc(1))), Iloc(1)
          write(Uout,'(I5,F10.5,2X,A,I0)') &
            Istep, dsqrt(maxval(force_norm)), trim(label(Iloc(1))), Iloc(1)
        end do
      close(Uinp)
    end do
    close(Uout)

  end subroutine force_analysis
! ++++++++++++++++++++++
! +++ force_analysis +++
! ++++++++++++++++++++++



! +++++++++++++++++++++++
! +++ charge_analysis +++
! +++++++++++++++++++++++
  subroutine charge_analysis
    integer :: Istep, i, j, k, Ifile
    character(len=128) :: name_out
    allocate(charge(Natom,Nbeads,TNstep))

    Istep = 0
    do Ifile = 1, Nfile
      PathFile=trim(DirResult(Ifile))//'/charge.dat'

      open(newunit=Uinp, file=PathFile, status='old', iostat=ierr)
        if ( ierr > 0 ) then
          print *, 'Check the path : ', PathFile
          stop 'ERROR!!: There is no "charge.dat"'
        end if

        read(Uinp,'()')
        do i = 1, Ncut(Ifile)
          read(Uinp,'()')
          do j = 1, Nbeads
            read(Uinp,'()')
          end do
        end do

        do k = Ncut(Ifile)+1, Nstep(Ifile)
          Istep = Istep + 1
          read(Uinp,'()')
          do j = 1, Nbeads
            read(Uinp,*) charge(:,j,Istep)
          end do
        end do
      close(Uinp)
    end do

    if ( jobtype == 62 ) then

    if ( save_beads .eqv. .True. ) then
      open(newunit=Uout, file=FNameBinary1, form='unformatted', access='stream', status='replace')
        do Istep = 1, TNstep
          do j = 1, Nbeads
            write(Uout) charge(atom1,j,Istep)
          end do
        end do
      close(Uout)
    end if
    write(*,*) "***** atomic charge of ", atom1, "is saved *****"
    write(*,*) "***** in ", FNameBinary1, " *****"
    data_beads = charge(atom1,:,:)
    name_out = "hist_charge.out"
    call calc_1Dhist(out_hist=name_out)


    open(newunit=Uout,file="step_charge.out",status='replace')
      do k = 1, TNstep
        if (mod(k,graph_step) == 0 ) write(Uout,'(I7,F10.5)') k, sum(charge(atom1,:,k))/dble(Nbeads)
      end do
    close(Uout)
    end if

    do i = 1, Natom
      print '(I3,F12.7)', i, sum(charge(i,:,:)) / dble(TNstep*Nbeads)
    end do
  end subroutine charge_analysis
! +++++++++++++++++++++++++++
! +++ end charge_analysis +++
! +++++++++++++++++++++++++++

! +++++++++++++++++++++++
! +++ dipole_analysis +++
! +++++++++++++++++++++++
  subroutine dipole_analysis
    integer :: i, j, k, Istep, Ifile
    real(8) :: abs_dipole
    allocate(dipole(3,Nbeads,TNstep))

    Istep = 0
    do Ifile = 1, Nfile
      PathFile=trim(DirResult(Ifile))//'/dipole.dat'

      open(newunit=Uinp, file=PathFile,status='old',iostat=ierr)
        if ( ierr > 0 ) then
          print *, 'Check the path : ',PathFile 
          stop 'ERROR!!: There is no "dipole.dat"'
        end if

        do k = 1, Ncut(Ifile)
          read(Uinp, '()')
          do j = 1, Nbeads
            read(Uinp,*) dipole(:,j,k)
          end do
        end do

        do k = Ncut(Ifile)+1, Nstep(Ifile)
          Istep = Istep + 1
          read(Uinp, '()')
          do j = 1, Nbeads
            read(Uinp,*) dipole(:,j,Istep)
          end do
        end do
      close(Uinp)
    end do

    abs_dipole = 0.0d0
    do Istep = 1, TNstep
      do j = 1, Nbeads
        data_beads(j,Istep) = dsqrt( dot_product(dipole(:,j,Istep),dipole(:,j,Istep)) )
      end do
    end do
    abs_dipole = sum(data_beads)/dble(TNstep*Nbeads)
    print '(a)', "The absolute dipole moment"
    print *, abs_dipole, " D "

    if ( save_beads .eqv. .True. ) then
      open(newunit=Uout,file=FNameBinary1, form='unformatted', access='stream', status='replace')
        do Istep = 1, TNstep
          do i = 1, Nbeads
            write(Uout) data_beads(i,Istep)
          end do
        end do
      close(Uout)
    end if

    return
  end subroutine dipole_analysis
! +++++++++++++++++++++++++++
! +++ end dipole_analysis +++
! +++++++++++++++++++++++++++

! +++++++++++++++++++++
! +++ hfcc_analysis +++
! +++++++++++++++++++++
  subroutine hfcc_analysis
    integer :: i, j, k, Istep, Ifile
    character(len=128) :: name_out
    allocate(hfcc(Natom,Nbeads,TNstep))

    print '(" ***** START HFCC analysis *****")'

    Istep = 0
    do Ifile = 1, Nfile
      PathFile=trim(DirResult(Ifile))//'/hfcc.dat'

      open(newunit=Uinp, file=PathFile, status='old', iostat=ierr)
        if ( ierr > 0 ) then
          print *, 'Check the path : ',PathFile 
          stop 'ERROR!!: There is no "hfcc.dat"'
        end if

        read(Uinp,'()')   ! Skip Header
        do i = 1, Ncut(Ifile)
          read(Uinp,'()',end=900)
          do j = 1, Nbeads
            read(Uinp,'()')
          end do
        end do

        do k = Ncut(Ifile)+1, Nstep(Ifile)
          Istep = Istep + 1
          read(Uinp,'()',end=900)
          do j = 1, Nbeads
            read(Uinp,*) hfcc(:,j,Istep)
          end do
        end do
      close(Uinp)
    end do

    if ( save_beads .eqv. .True. ) then
      open(newunit=Uout, file=FNameBinary1, form='unformatted', access='stream', status='replace')
        do Istep = 1, TNstep
          do j = 1, Nbeads
            write(Uout) hfcc(atom1,j,Istep)
          end do
        end do
      close(Uout)
      write(*,*) '    Binary data is saved in "', FNameBinary1,'"    '
    end if
    write(*,'("     HFCC of atom",I0, " is saved    ")')  atom1

    data_beads = hfcc(atom1,:,:)
    name_out = "hist_hfcc.out"
    do k = 1, TNstep
      data_step(k) = sum(data_beads(:,k))/dble(Nbeads)
    end do

    call calc_1Dhist(out_hist=name_out)

    open(newunit=Uout,file="step_hfcc.out",status='replace')
      do k = 1, TNstep
        if (mod(k,graph_step) == 0 ) write(Uout,'(I7,F10.5)') k, data_step(k) !sum(hfcc(atom1,:,k))/dble(Nbeads)
      end do
    close(Uout)

    write(*,*) "*** All the HFCCs are as follows ***"
    write(*,*) "   Num. HFCC"
    do i = 1, Natom
      print '("   ",I4,F12.7)', i, sum(hfcc(i,:,:)) / dble(TNstep*Nbeads)
    end do

    call reblock_step()

    return
900 call err_line_exceed(k*(Nbeads+1) + j, k)
  end subroutine hfcc_analysis
! +++++++++++++++++++++++++
! +++ end hfcc_analysis +++
! +++++++++++++++++++++++++

subroutine err_line_exceed(Iline, Istep)
  integer, intent(in) :: Iline, Istep
  print *, 'The line is ', Iline
  print *, 'The step is ', Istep
  stop 'ERROR!!: Reading line exceed the coor lines'
end subroutine err_line_exceed

end module mod_other_quantities
! you can change text to line




