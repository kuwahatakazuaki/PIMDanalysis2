!
! The weight of hydrogen atoms are change to zero 
subroutine rotation
  use input_parameter,  only: label, TNstep, save_beads, Nbeads, Natom, &
      FNameBinary1, graph_step, weight, r_ref, jobtype, label, Ndiv, &
      muon => atom_cube, Nhyd, hyd, r
  use utility,          only: calc_deviation, calc_cumulative, get_rot_mat, atom2num
  implicit none
  real(8), parameter :: pi = 4.0d0*atan(1.0d0)
  integer, parameter :: N = 4
  integer, parameter :: liwork = 5*N+3,  lwork = 2*N*N+6*N+1
  real(8) :: work(lwork)
  integer :: iwork(liwork)
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

  open(newunit=Uout,file='eigen.out')
    do k = 1, TNstep
      if ( mod(k,graph_step) == 0 ) then
        write(Uout,9999) eigenArray(k)
      end if
    end do
  close(Uout)

  select case(jobtype)
    case(71)
      call save_movie
    case(72)
      call save_cube_sub
  end select
  print '(a)',  " ***** END Removing rotation freedom *****"

  9999 format(E12.5)
contains

  ! Please reffer the save_cube in utility, in the futhure 
  subroutine save_cube_sub
    integer :: Uout
!    integer, parameter :: Ndiv = 20
    real(8), parameter :: Ledge = 10.0d0
    real(8), parameter :: Bohr2Angs = 0.529177249
    real(8), parameter :: Angs2Bohr = 1.8897259886
    real(8), parameter :: margine = 1d-1
    real(8) :: grid(Ndiv,Ndiv,Ndiv)
    real(8) :: Lmin(3), Lmax(3)
    real(8) :: dL(3), base_vec(3,3)
    integer, allocatable :: coun(:,:,:)
    integer :: cx,cy,cz

    rnew(:,:,:,:) = rnew(:,:,:,:) * Angs2Bohr
    do i = 1, 3
      Lmin(i) = minval(rnew(i,muon,:,:)) - margine
      Lmax(i) = maxval(rnew(i,muon,:,:)) + margine
    end do
    dL(:) = (Lmax(:) - Lmin(:)) / dble(Ndiv)

    print '(a)',        '   Constructing the cube file '
    print '(a,I4)',     '     Ndiv     = ', Ndiv
    print '(a,I4)',     '     cube atom= ', muon
    print '(a,1pe11.3)','     margine  = ', margine
    print '(a,3F10.5)', '     dL       = ', dL(:)
    print '(a,3F10.5)', '     Lmin     = ', Lmin(:)
    base_vec(:,:) = 0.0d0
    do i = 1, 3
      base_vec(i,i) = dL(i)
    end do

    allocate(coun(3,Nbeads,TNstep))
    do k = 1, TNstep
      do j = 1, Nbeads
        coun(:,j,k) = int( ( rnew(:,muon,j,k)-Lmin(:) ) / dL(:) ) + 1
      end do
    end do

    grid(:,:,:) = 0.0d0
    do k = 1, TNstep
      do j = 1, Nbeads
          cx = coun(1,j,k)
          cy = coun(2,j,k)
          cz = coun(3,j,k)
          grid(cx,cy,cz) = grid(cx,cy,cz) + 1.0d0
      end do
    end do

    grid(:,:,:) = grid(:,:,:) / dble(TNstep*Nbeads)
    open(newunit=Uout,file='hyd.cube',status='replace')
      write(Uout,*) "commnet"
      write(Uout,*) "commnet"
      write(Uout,9999) Natom-Nhyd, Lmin(:)
      do i = 1, 3
        write(Uout,9999) Ndiv, base_vec(i,:)
      end do
      j = 1
      do i = 1, Natom
        if ( i == hyd(j) ) then
          j = j + 1
          cycle
        end if
        write(Uout,9999) atom2num(trim(label(i))), dble(i), &
                         [sum(rnew(1,i,:,:)),sum(rnew(2,i,:,:)),sum(rnew(3,i,:,:))]/dble(TNstep*Nbeads)
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
    print '(a)', '   Cube file is saved in "hyd.cube"'

  9998  format(I5,4F12.6)
  9999  format(I5,4F12.6)
  end subroutine save_cube_sub

  subroutine save_movie
    integer :: Uout
    open(newunit=Uout,file='vmd.xyz',status='replace')
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

