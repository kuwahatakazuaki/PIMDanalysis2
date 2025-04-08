!
! The weight of hydrogen atoms are change to zero 
subroutine rotation
  use input_parameter,  only: label, TNstep, save_beads, Nbeads, Natom, &
      FNameBinary1, graph_step, weight, r_ref, jobtype, label, Ndiv, &
      muon => atom_cube, Nhyd, hyd, r, AngtoAU, AUtoAng
  use utility,          only: calc_deviation, calc_cumulative, get_rot_mat, atom2num, save_cube
  implicit none
  real(8), parameter :: pi = 4.0d0*atan(1.0d0)
  integer :: i, j, k, xyz, info, Uout
  real(8), allocatable :: rnew(:,:,:,:), eigenArray(:)
  real(8) :: rot(3,3), qua(4), rave(3), Tweight
  real(8) :: matA(4,4), matB(4,4), eigen(4)

  allocate(rnew(3,Natom,Nbeads,TNstep))
  allocate(eigenArray(TNstep))

  print '(a)',  " ***** START Removing rotation freedom *****"
  do j = 1, Nhyd
    weight(hyd(j)) = 0.0d0
  end do
  Tweight = sum(weight(:))
  do xyz = 1, 3
    rave(xyz) = dot_product(r_ref(xyz,:), weight(:)) / (Tweight)
  end do
  r_ref(:,:) = r_ref(:,:) - spread(rave(:), dim=2, ncopies=Natom)

  select case(jobtype)
    case(71)
      call remove_rotation_allbeads ! older version
    case(72:73)
      call remove_rotation_eachbeads ! newer version
  end select

  !call save_cube_sub
  if ( jobtype == 73 ) then
    call save_cube_FeH()
  else
    call save_cube(rnew,[muon],'cube.cube')
  end if

  open(newunit=Uout,file='eigen.out')
    do k = 1, TNstep
      if ( mod(k,graph_step) == 0 ) then
        write(Uout,9999) eigenArray(k)
      end if
    end do
  close(Uout)

  !select case(jobtype)
  !  case(71)
  !    call save_movie
  !  case(72)
  !    call save_cube_sub
  !end select
  print '(a)',  " ***** END Removing rotation freedom *****"

  9999 format(E12.5)
contains

  subroutine save_cube_FeH()
    integer, parameter :: Nf = 6
    real(8) :: rH(3,Nbeads,TNstep*4)
    real(8) :: rave(3,5) ! 1:H, 2:F1, 3,F2,
    real(8) :: center(3), rF(3,Nf)
    real(8) :: Lmin(3), Lmax(3), dL(3)
    !real(8) :: base_vec(3,3)
    real(8) :: grid(Ndiv,Ndiv,Ndiv)
    real(8), parameter :: margine = 1d-1
    integer :: xyz, Iatom, i, Ibead, Istep
    character(len=*), parameter :: Fout = 'oct.cube'

    rnew(:,:,:,:) = rnew(:,:,:,:) * AngtoAU
    do Iatom = 1, Natom
      do xyz = 1, 3
        rave(xyz,Iatom) = sum(rnew(xyz,Iatom,:,:))/dble(TNstep*Nbeads)
      end do
    end do
    center(:) = 0.5d0 * ( rave(:,2)+rave(:,3) )
    do xyz = 1, 3
      rnew(xyz,:,:,:) = rnew(xyz,:,:,:) - center(xyz)
      rave(xyz,:) = rave(xyz,:) - center(xyz)
    end do
    rF(:,:) = 0.0d0
    rF(3,1) = rave(3,2)
    rF(3,2) = (-1.0d0)*rF(3,1)
    rF(1,3) = 0.5d0*(rave(1,4)+rave(2,5))
    rF(2,4) = rF(1,3)
    rF(1,5) = (-1.0d0)*rF(1,3)
    rF(2,6) = (-1.0d0)*rF(1,3)

    do i = 1, 4
      rH(:,:,TNstep*(i-1)+1:TNstep*i) = rnew(:,1,:,:)
    end do
    i = 2; rH(1,:,TNstep*(i-1)+1:TNstep*i) = rH(1,:,TNstep*(i-1)+1:TNstep*i) * (-1.0d0)
    i = 3; rH(2,:,TNstep*(i-1)+1:TNstep*i) = rH(2,:,TNstep*(i-1)+1:TNstep*i) * (-1.0d0)
    i = 4; rH(1,:,TNstep*(i-1)+1:TNstep*i) = rH(1,:,TNstep*(i-1)+1:TNstep*i) * (-1.0d0); &
           rH(2,:,TNstep*(i-1)+1:TNstep*i) = rH(2,:,TNstep*(i-1)+1:TNstep*i) * (-1.0d0)
    do xyz = 1, 3
      Lmin(xyz) = minval(rH(xyz,:,:)) - margine
      Lmax(xyz) = maxval(rH(xyz,:,:)) + margine
    end do
    Lmax(:) = 0.5d0*(Lmax(:)-Lmin(:))
    Lmin(:) = -Lmax(:)
    dL(:) = (Lmax(:)-Lmin(:)) / dble(Ndiv)

    block
      integer :: coun(3,Nbeads,TNstep*4)
      integer :: cx, cy, cz
      grid(:,:,:) = 0.0d0
      do Istep = 1, TNstep*4
      do Ibead = 1, Nbeads
        !coun(:,Ibead,Istep) = int( (rH(:,Ibead,Istep)-Lmin(:))/dL(:) ) + 1
        coun(:,Ibead,Istep) = int( (rH(:,Ibead,Istep)-Lmin(:))/dL(:) ) + 1
      end do
      end do
      !do Istep = 1, TNstep*2
      do Istep = 1, TNstep*4
      do Ibead = 1, Nbeads
        if ( all(coun(:,Ibead,Istep) >= 1 .and. coun(:,Ibead,Istep) <= Ndiv )  ) then
        cx = coun(1,Ibead,Istep)
        cy = coun(2,Ibead,Istep)
        cz = coun(3,Ibead,Istep)
        grid(cx,cy,cz) = grid(cx,cy,cz) + 1.0d0
        else
          print *, Ibead, Istep, ":", coun(:,Ibead,Istep)
        end if
      end do
      end do
      grid(:,:,:) = grid(:,:,:) / dble(4*TNstep*Nbeads)
    end block

    open(newunit=Uout,file=Fout,status='replace')
      write(Uout,'(a,3F9.4,a)') "# average of H ", rave(:,1)*AUtoAng, " in Angstrom"
      write(Uout,'(a,3F9.4,a)') "# average of Fe", rave(:,3)*AUtoAng, " in Angstrom"
      write(Uout,9999) Nf, Lmin(:) + dL(:)*0.5d0
      write(Uout,9999) Ndiv, dL(1), 0.0,    0.0
      write(Uout,9999) Ndiv, 0.0,   dL(2),  0.0
      write(Uout,9999) Ndiv, 0.0,   0.0,    dL(3)
      do Iatom = 1, Nf
        write(Uout,9999) 26, dble(Iatom), rF(:,Iatom)
      end do

      do i = 1, Ndiv
        do j = 1, Ndiv
          do k = 1, Ndiv
            write(Uout,'(E13.5)',advance='no') grid(i,j,k)
            if ( mod(k,6) == 0 ) write(Uout,*)
          end do
          write(Uout,*)
        end do
      end do
    close(Uout)

    block
    real(8) :: rx, ry
    real(8) :: grid_cut(Ndiv,Ndiv)
    integer :: mid0, Imax(1)
    mid0 = Ndiv/2 + 1
    do i = 1, Ndiv
      do j = mid0, Ndiv
        grid_cut(i,j) = sum(grid(i,j,mid0-1:mid0+1) + grid(Ndiv-i+1,j,mid0-1:mid0+1))/6.0d0
      end do
    end do

    open(newunit=Uout,file='plot_xy.dat',status='replace')
      write(Uout,'(a,3F9.4,a)') "# average of H ", rave(:,1)*AUtoAng, " in Angstrom"
      write(Uout,'(a,3F9.4,a)') "# average of Fe", rave(:,3)*AUtoAng, " in Angstrom"
      do i = 1, Ndiv
        write(Uout,*) "# index of x", i
        rx = ( Lmin(1)+dL(1)*(dble(i)-0.5d0) ) * AUtoAng
        do j = mid0, Ndiv
          ry = ( Lmin(2)+dL(2)*(dble(j)-0.5d0) ) * AUtoAng
          write(Uout,*) rx, ry, grid_cut(i,j)
          !write(Uout,*) rx, ry, sum(grid(i,j,mid0-1:mid0+1) + grid(Ndiv-i+1,j,mid0-1:mid0+1))/6.0d0
        end do
        write(Uout,*) ""
      end do
    close(Uout)

    open(newunit=Uout,file='max_xy.dat',status='replace')
      write(Uout,'(a,3F9.4,a)') "# average of H ", rave(:,1)*AUtoAng, " in Angstrom"
      write(Uout,'(a,3F9.4,a)') "# average of Fe", rave(:,3)*AUtoAng, " in Angstrom"
      do i = 1, Ndiv
        rx = ( Lmin(1)+dL(1)*(dble(i)-0.5d0) ) * AUtoAng
        !Imax(:) = maxloc( grid(i,mid0:,mid0) ) + mid0-1
        Imax(:) = maxloc( grid_cut(i,mid0:) ) + mid0-1
        ry = ( Lmin(2)+dL(2)*(dble(Imax(1))-0.5d0) ) * AUtoAng
        write(Uout,*) rx, ry, Imax(1)
      end do
    close(Uout)

    end block

  9998  format(I5,4F12.6)
  9999  format(I5,4F12.6)
  end subroutine save_cube_FeH

  subroutine remove_rotation_eachbeads
  integer :: Istep, Ibead
  integer, parameter :: N = 4
  integer, parameter :: liwork = 5*N+3,  lwork = 2*N*N+6*N+1
  real(8) :: work(lwork)
  integer :: iwork(liwork)
! --- Start Remove center of mass ---
  Loop_step1:do Istep = 1, TNstep
    rave(:) = 0.0d0
    do j = 1, Nbeads
      do xyz = 1, 3
        rave(xyz) = rave(xyz) + dot_product(r(xyz,:,j,Istep), weight(:))
      end do
      rave(:) = rave(:) / Tweight
      r(:,:,j,Istep) = r(:,:,j,Istep) - spread(rave(:),dim=2,ncopies=Natom)
    end do
  end do Loop_step1
! --- End Remove center of mass ---

! --- Start Rotation to r_ref ---
  Loop_step2:do Istep = 1, TNstep
    matB(:,:) = 0.0d0
    do j = 1, Nbeads
      do i = 1, Natom
        matA(:,:) = make_matA(r_ref(:,i)+r(:,i,j,Istep),r_ref(:,i)-r(:,i,j,Istep))
        matB(:,:) = matB(:,:) + weight(i)*matmul(transpose(matA),matA)
      end do
      matB(:,:) = matB(:,:) / Tweight

      call dsyevd('V', 'U', N, matB, N, eigen, work, lwork, iwork, liwork, info)
      qua(:) = matB(:,1)
      eigenArray(Istep) = eigen(1)
      rot(:,:) = get_rot_mat(qua(:))

      do i = 1, Natom
        rnew(:,i,j,Istep) = matmul(rot(:,:), r(:,i,j,Istep))
      end do
    end do
  end do Loop_step2
! --- End Rotation to r_ref ---
  end subroutine remove_rotation_eachbeads

! =========================================================
! === Older version, this will be deleted in the future ===
! =========================================================
  subroutine remove_rotation_allbeads
  integer, parameter :: N = 4
  integer, parameter :: liwork = 5*N+3,  lwork = 2*N*N+6*N+1
  real(8) :: work(lwork)
  integer :: iwork(liwork)
! --- Start Remove center of mass ---
  step_loop1:do k = 1, TNstep
    rave(:) = 0.0d0
    do j = 1, Nbeads
      do xyz = 1, 3
        rave(xyz) = rave(xyz) + dot_product(r(xyz,:,j,k), weight(:))
      end do
    end do
    rave(:) = rave(:) / (Tweight * dble(Nbeads))
    r(:,:,:,k) = r(:,:,:,k) - spread( spread(rave(:),dim=2,ncopies=Natom),dim=3,ncopies=Nbeads )
  end do step_loop1
! --- End Remove center of mass ---

! --- Start Rotation to r_ref ---
  step_loop2:do k = 1, TNstep
    matB(:,:) = 0.0d0
    do j = 1, Nbeads
      do i = 1, Natom
        matA(:,:) = make_matA(r_ref(:,i)+r(:,i,j,k),r_ref(:,i)-r(:,i,j,k))
        matB(:,:) = matB(:,:) + weight(i)*matmul(transpose(matA),matA)
      end do
    end do
    matB(:,:) = matB(:,:) / (Tweight*dble(Nbeads))

    call dsyevd('V', 'U', N, matB, N, eigen, work, lwork, iwork, liwork, info)
    qua(:) = matB(:,1)
    eigenArray(k) = eigen(1)
    rot(:,:) = get_rot_mat(qua(:))
    do j = 1, Nbeads
      do i = 1, Natom
        rnew(:,i,j,k) = matmul(rot(:,:), r(:,i,j,k))
      end do
    end do
  end do step_loop2
! --- End Rotation to r_ref ---
  end subroutine remove_rotation_allbeads

!  ! Please reffer the save_cube in utility, in the futhure 
!  subroutine save_cube_sub
!    integer :: Uout
!    real(8), parameter :: Ledge = 10.0d0
!    real(8), parameter :: Angs2Bohr = 1.8897259886d0
!    real(8), parameter :: margine = 1d-1
!    real(8) :: grid(Ndiv,Ndiv,Ndiv)
!    real(8) :: Lmin(3), Lmax(3)
!    real(8) :: dL(3), base_vec(3,3)
!    integer, allocatable :: coun(:,:,:)
!    integer :: cx,cy,cz
!
!    rnew(:,:,:,:) = rnew(:,:,:,:) * Angs2Bohr
!    do i = 1, 3
!      Lmin(i) = minval(rnew(i,muon,:,:)) - margine
!      Lmax(i) = maxval(rnew(i,muon,:,:)) + margine
!    end do
!    dL(:) = (Lmax(:) - Lmin(:)) / dble(Ndiv)
!
!    print '(a)',        '   Constructing the cube file '
!    print '(a,I4)',     '     Ndiv     = ', Ndiv
!    print '(a,I4)',     '     cube atom= ', muon
!    print '(a,1pe11.3)','     margine  = ', margine
!    print '(a,3F10.5)', '     dL       = ', dL(:)
!    print '(a,3F10.5)', '     Lmin     = ', Lmin(:)
!
!    base_vec(:,:) = 0.0d0
!    do i = 1, 3
!      base_vec(i,i) = dL(i)
!    end do
!
!    allocate(coun(3,Nbeads,TNstep))
!    do k = 1, TNstep
!      do j = 1, Nbeads
!        coun(:,j,k) = int( ( rnew(:,muon,j,k)-Lmin(:) ) / dL(:) ) + 1
!      end do
!    end do
!
!    grid(:,:,:) = 0.0d0
!    do k = 1, TNstep
!      do j = 1, Nbeads
!          cx = coun(1,j,k)
!          cy = coun(2,j,k)
!          cz = coun(3,j,k)
!          grid(cx,cy,cz) = grid(cx,cy,cz) + 1.0d0
!      end do
!    end do
!
!    grid(:,:,:) = grid(:,:,:) / dble(TNstep*Nbeads)
!    open(newunit=Uout,file='hyd.cube',status='replace')
!      write(Uout,*) "commnet"
!      write(Uout,*) "commnet"
!      write(Uout,9999) Natom-Nhyd, Lmin(:)
!      do i = 1, 3
!        write(Uout,9999) Ndiv, base_vec(i,:)
!      end do
!      j = 1
!      do i = 1, Natom
!        if ( i == hyd(j) ) then
!          j = j + 1
!          cycle
!        end if
!        write(Uout,9999) atom2num(trim(label(i))), dble(i), &
!                         [sum(rnew(1,i,:,:)),sum(rnew(2,i,:,:)),sum(rnew(3,i,:,:))]/dble(TNstep*Nbeads)
!      end do
!      do i = 1, Ndiv
!        do j = 1, Ndiv
!          do k = 1, Ndiv
!            write(Uout,'(E13.5)',advance='no') grid(i,j,k)
!            if ( mod(k,6) == 0 ) write(Uout,*)
!          end do
!          write(Uout,*)
!        end do
!      end do
!    close(Uout)
!    print '(a)', '   Cube file is saved in "hyd.cube"'
!
!  9998  format(I5,4F12.6)
!  9999  format(I5,4F12.6)
!  end subroutine save_cube_sub

  subroutine save_movie
    integer :: Uout
    open(newunit=Uout,file='movie.xyz',status='replace')
      do k = 1, TNstep
        if ( mod(k,10) == 0) then
          write(Uout,'(I10)') Natom*Nbeads
          write(Uout,'(I10)') k
          do j = 1, Nbeads
            do i = 1, Natom
              write(Uout,9999) label(i), rnew(:,i,j,k)
            end do
          end do
        end if
      end do
    close(Uout)
    9999 format(a,4F11.7)
  end subroutine save_movie

  function make_matA(x,y) result(mat)
    real(8) :: x(3), y(3)
    real(8) :: mat(4,4)
    mat(1,:) = [0.d0,  -y(1), -y(2), -y(3)]
    mat(2,:) = [y(1),   0.d0, -x(3),  x(2)]
    mat(3,:) = [y(2),   x(3),  0.d0, -x(1)]
    mat(4,:) = [y(3),  -x(2),  x(1),  0.d0]
  end function make_matA
end subroutine rotation

