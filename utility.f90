module utility
  use input_parameter,  only: data_beads, data_step, TNstep, graph_step, Ndiv, label
  implicit none
  private
  real(8), parameter :: pi = atan(1.0d0)*4.0d0
  real(8), parameter :: real_max = huge(real(8))
  real(8), parameter :: real_min = tiny(real(8))

  public :: reblock_step, get_rot_mat, calc_cumulative, calc_deviation, get_inv_mat, pi, &
            rand3, random_seed_ini, get_volume, get_qua_theta, norm, lowerchr, outer_product, &
            atom2num, save_cube, real_max, real_min, sort_real
contains

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++ Start save_cube ++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine save_cube(rcub,Iatoms,Fout)
    integer, intent(in) :: Iatoms(:)
    character(len=*), optional :: Fout
    real(8), intent(inout) :: rcub(:,:,:,:) ! rcub(3,Natom,Nbeads,TNstep)
    real(8), parameter :: Ledge = 10.0d0
    real(8), parameter :: Bohr2Angs = 0.529177249
    real(8), parameter :: Angs2Bohr = 1.8897259886
    real(8), parameter :: margine = 1d-1
    real(8) :: grid(Ndiv,Ndiv,Ndiv)
    real(8) :: Lmin(3), Lmax(3)
    real(8) :: dL(3), base_vec(3,3)
    integer, allocatable :: coun(:,:,:,:)
    integer :: Uout,i,j,k
    integer :: uboun(4), Natom, Nbeads, TNstep, Ncube

    uboun(:) = ubound(rcub)
    Natom  = uboun(2)
    Nbeads = uboun(3)
    TNstep = uboun(4)
    Ncube  = size(Iatoms)

    rcub(:,:,:,:) = rcub(:,:,:,:) * Angs2Bohr

    block
      real(8), allocatable :: temp1(:), temp2(:)
      allocate(temp1(Ncube), temp2(Ncube))
      do i = 1, 3
        do j = 1, Ncube
          temp1(j) = minval(rcub(i,Iatoms(j),:,:))
          temp2(j) = maxval(rcub(i,Iatoms(j),:,:))
        end do
        Lmin(i) = minval(temp1(:)) - margine
        Lmax(i) = maxval(temp2(:)) + margine
      end do
    end block

    dL(:) = (Lmax(:) - Lmin(:)) / dble(Ndiv)

    print '(a)',        '  *** Constructing the cube file ***'
    print '(a,I4)',     '     Ndiv      = ', Ndiv
    print '(a,1pe11.3)','     margine   = ', margine
    print '(a,*(I4))',  '     cube atom = ', Iatoms(:)
    print '(a,3F10.5)', '     dL        = ', dL(:)
    print '(a,3F10.5)', '     Lmin      = ', Lmin(:)
    base_vec(:,:) = 0.0d0
    do i = 1, 3
      base_vec(i,i) = dL(i)
    end do

    allocate(coun(3,Ncube,Nbeads,TNstep))
    do i = 1, Ncube
      do k = 1, TNstep
        do j = 1, Nbeads
          coun(:,i,j,k) = int( ( rcub(:,Iatoms(i),j,k)-Lmin(:) ) / dL(:) ) + 1
        end do
      end do
    end do

    block
      integer :: cx,cy,cz
      grid(:,:,:) = 0.0d0
      do i = 1, Ncube
        do k = 1, TNstep
          do j = 1, Nbeads
              cx = coun(1,i,j,k)
              cy = coun(2,i,j,k)
              cz = coun(3,i,j,k)
              grid(cx,cy,cz) = grid(cx,cy,cz) + 1.0d0
          end do
        end do
      end do
    end block

    grid(:,:,:) = grid(:,:,:) / dble(TNstep*Nbeads)
    open(newunit=Uout,file='cube.cube',status='replace')
      write(Uout,*) "commnet"
      write(Uout,*) "commnet"
      write(Uout,9999) Natom-Ncube, Lmin(:)
      do i = 1, 3
        write(Uout,9999) Ndiv, base_vec(i,:)
      end do
      j = 1
      do i = 1, Natom
        if ( i == Iatoms(j) ) then
          j = j + 1
          cycle
        end if
        !write(Uout,9999) 8, dble(i), &  ! Only for Tetrahedra O-H4
        write(Uout,9999) atom2num(trim(label(i))), dble(i), &
                         [sum(rcub(1,i,:,:)),sum(rcub(2,i,:,:)),sum(rcub(3,i,:,:))]/dble(TNstep*Nbeads)
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
    print '(a)', '  *** Cube file is saved in "cube.cube" ***'

  9998  format(I5,4F12.6)
  9999  format(I5,4F12.6)

  end subroutine save_cube
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++ End save_cube ++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++ Start Reblocking  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine reblock_step()
    integer :: Istep, j
    integer :: Neach, Ndata, Uout
    integer :: N1!, N2
    real(8) :: ave, var, dev, err
    real(8), allocatable :: data1(:), data2(:)

    N1 = size(data_step)
    allocate(data1(N1),data2(N1))
    data1(:) = data_step(:)
    Neach = N1
    Istep = 0

    print '(a)', " ***** START Reblocking methods *****"
    open(newunit=Uout,file='reblock.out')
      do
        Istep = Istep + 1
        print *, '   Neach is :', Neach
        if ( mod(Neach,2) == 0 ) then
          Neach = Neach/2
        else
          Neach = (Neach-1)/2
        end if
        if ( Neach <= 1) exit
        do j = 1, Neach
          data2(j) = 0.5d0 * ( data1(2*j-1)+data1(2*j) )
        end do

        ave = sum(data2(1:Neach))/dble(Neach)
        var = sum( data2(1:Neach)*data2(1:Neach) ) / dble(Neach) - ave**2
        dev = dsqrt(var/dble(Neach-1))
        err = dev / dsqrt(2.0d0*dble(Neach-1))
        write(Uout,'(I10,2E13.5)') Istep, dev, err

        data1(1:Neach) = data2(1:Neach)
      end do
    close(Uout)
    print '(a)',
    print '(a)', '   Reblocking data is saved in "reblock.out"'
    print '(a)', '   Type "plot "reblock.out" with errorbar" in Gnuplot'
    print '(a)', " ***** End Reblocking methods*****"
    print '(a)',
  end subroutine reblock_step
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++ End Reblocking  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++ Start Quaternion  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  function get_rot_mat(qua) result(rot)
    real(8) :: qua(4)
    real(8) :: rot(3,3)
    real(8) :: q0, q1, q2, q3
    q0 = qua(1)
    q1 = qua(2)
    q2 = qua(3)
    q3 = qua(4)
    rot(1,1) = 1.0d0 - 2.0d0 * ( q2**2 + q3**2)
    rot(1,2) = 2.0d0 * ( q1*q2 - q0*q3 )
    rot(1,3) = 2.0d0 * ( q1*q3 + q0*q2 )

    rot(2,1) = 2.0d0 * ( q2*q1 + q0*q3 )
    rot(2,2) = 1.0d0 - 2.0d0 * ( q3**2 + q1**2)
    rot(2,3) = 2.0d0 * ( q2*q3 - q0*q1 )
!    rot(2,3) = 2.0d0 * ( q1*q2 - q0*q1 )

    rot(3,1) = 2.0d0 * ( q3*q1 - q0*q2 )
    rot(3,2) = 2.0d0 * ( q3*q2 + q0*q1 )
    rot(3,3) = 1.0d0 - 2.0d0 * ( q1**2 + q2**2)
  end function get_rot_mat

  function get_qua_theta(theta, direc)  result(qua)
    real(8) :: theta, direc(3)
    real(8) :: qua(4)
    real(8) :: half
    half = 0.5d0 * theta * pi / 180.0d0
    qua(1) = dcos(half)
    qua(2:4) = direc(:) * dsin(half)
  end function
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++ End Quaternion  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++ Start calc_cumulative NEW verions ++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine calc_cumulative()
    integer :: i, cumu_step = 100
    real(8) :: data_dev, data_err
    character(len=:), allocatable :: out_cumulative
    out_cumulative="cumu_new.out"

    open(20, file=out_cumulative, status='replace')
    do i = 1, TNstep
      if (mod(i,graph_step*cumu_step) == 0) then
        call calc_deviation(data_dev, data_err, end_step=i)
        write(20,'(I6,3F13.6)') i, sum(data_step(1:i))/dble(i), data_dev, data_err
      end if
    end do
    close(20)
  end subroutine calc_cumulative
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++ End calc_cumulative  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++ Start calc_cumulative Older verions ++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine calc_cumulative_old()
!    use input_parameter, only: TNstep, graph_step, data_step
    integer :: i, cumu_step = 100
    real(8) :: data_dev, data_err
    character(len=:), allocatable :: out_cumulative
    out_cumulative="cumu_old.out"

    open(20, file=out_cumulative, status='replace')
    do i = 1, TNstep
      if (mod(i,graph_step*cumu_step) == 0) then
        call calc_deviation(data_dev, data_err, end_step=i)
        write(20,'(I6,3F13.6)') i, sum(data_step(1:i))/dble(i), data_dev, data_err
      end if
    end do
    close(20)
  end subroutine calc_cumulative_old
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++ End calc_cumulative  +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++ Start calc_deviation  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine calc_deviation(data_dev, data_err, end_step)
    use input_parameter, only: Nbeads, TNstep, data_beads
    integer i, j, Nstep
    real(8) :: data_ave
    integer, intent(in), optional :: end_step
    real(8), intent(out) :: data_dev
    real(8), intent(out), optional :: data_err

    if (present(end_step)) then
      Nstep = end_step
    else
      Nstep = ubound(data_beads,2)
    end if

    data_ave = sum(data_beads)/size(data_beads)
    data_dev = 0.0d0
    do i = lbound(data_beads,2), Nstep
      do j = 1, Nbeads
        data_dev = data_dev + (data_beads(j,i) - data_ave)**2
      end do
    end do
    data_dev = data_dev / dble(Nstep*Nbeads)
    data_dev = dsqrt(data_dev)
    if (present(data_err)) data_err = data_dev / sqrt(dble(Nstep))
  end subroutine calc_deviation
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++ End calc_deviation  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_inv_mat(mat,inv,n)
    integer :: n
    real(8), intent(in)  :: mat(n,n)
    real(8), intent(out) :: inv(n,n)
    integer :: lwork, lda, info
    real(8), allocatable :: work(:)
    integer, allocatable :: ipiv(:)
    inv(:,:) = mat(:,:)
    lda = n
    lwork = 64*n
    allocate(work(lwork),ipiv(n))
    call dgetrf(N, N, inv, lda, ipiv, info)
    call dgetri(N, inv, lda, ipiv, work, lwork, info)
  end subroutine get_inv_mat

  function rand3() result(r3)
    real(8) :: r3(3)
    call random_number(r3(:))
    r3(:) = r3(:) - 0.5d0
  end function rand3

  subroutine random_seed_ini
    integer :: seedsize
    integer, allocatable :: seed(:)

    call random_seed(size=seedsize)
    allocate(seed(seedsize))
    seed(:) = 123456
    call random_seed(put=seed)
  end subroutine random_seed_ini

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++ Start norm +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  real(8) function norm(x)
    implicit none
    real(8), intent(in) :: x(:)
    norm = dsqrt( sum( x(:)*x(:)) )
  end function
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++ End norm +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++ Start get_volume +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  function get_volume(vec) result(vol)
    real(8) :: vec(3,3)
    real(8) :: vol
    vol = dot_product(vec(1,:), outer_product(vec(2,:),vec(3,:)))
  end function get_volume
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++ End get_volume +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  function lowerchr(str)
    character(*), intent(in) :: str
    character(len(str)) :: lowerchr
    integer :: i
    do i = 1, len_trim(str)
      if ( str(i:i) >= 'A' .and. str(i:i) <= 'Z' ) then
        lowerchr(i:i) = char(ichar(str(i:i))+32)
      else
        lowerchr(i:i) = str(i:i)
      end if
    end do
  end function lowerchr

  function atom2num(cha) result(num)
    character(*) :: cha
    integer :: num
    cha = lowerchr(cha)
    num = 0
    if     ( trim(cha) == 'h' ) then
      num = 1
    elseif ( trim(cha) == 'he' ) then
      num = 2
    elseif ( trim(cha) == 'li' ) then
      num = 3
    elseif ( trim(cha) == 'be' ) then
      num = 4
    elseif ( trim(cha) == 'b' ) then
      num = 5
    elseif ( trim(cha) == 'c' ) then
      num = 6
    elseif ( trim(cha) == 'n' ) then
      num = 7
    elseif ( trim(cha) == 'o' ) then
      num = 8
    elseif ( trim(cha) == 'f' ) then
      num = 9
    elseif ( trim(cha) == 'ne' ) then
      num = 10
    elseif ( trim(cha) == 'na' ) then
      num = 11
    elseif ( trim(cha) == 'mg' ) then
      num = 12
    elseif ( trim(cha) == 'al' ) then
      num = 13
    elseif ( trim(cha) == 'si' ) then
      num = 14
    elseif ( trim(cha) == 'p' ) then
      num = 15
    elseif ( trim(cha) == 's' ) then
      num = 16
    elseif ( trim(cha) == 'cl' ) then
      num = 17
    elseif ( trim(cha) == 'v' ) then
      num = 23
    elseif ( trim(cha) == 'cr' ) then
      num = 24
    elseif ( trim(cha) == 'mn' ) then
      num = 25
    elseif ( trim(cha) == 'fe' ) then
      num = 26
    elseif ( trim(cha) == 'co' ) then
      num = 27
    elseif ( trim(cha) == 'ni' ) then
      num = 28
    elseif ( trim(cha) == 'cu' ) then
      num = 29
    elseif ( trim(cha) == 'zn' ) then
      num = 30
    elseif ( trim(cha) == 'ga' ) then
      num = 31
    elseif ( trim(cha) == 'ge' ) then
      num = 32
    elseif ( trim(cha) == 'as' ) then
      num = 33
    else
      print *, cha, 'is not exist in atom2num'
      stop 'ERROR!! "atom2num" cannot chage '
    end if
  end function

  character(len=2) function itoc(i)
    integer :: i
    select case(i)
      case(1)
        itoc = 'H'
      case(3)
        itoc = 'Li'
      case(5)
        itoc = 'B'
      case(6)
        itoc = 'C'
      case(7)
        itoc = 'N'
      case(8)
        itoc = 'O'
      case(9)
        itoc = 'F'
    end select
  end function itoc

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++ Start outer_product ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  function outer_product(a,b) result(vec)
    real(8), intent(in) :: a(3), b(3)
    real(8) :: vec(3)

    vec(1) = a(2) * b(3) - a(3) * b(2)
    vec(2) = a(3) * b(1) - a(1) * b(3)
    vec(3) = a(1) * b(2) - a(2) * b(1)
  end function outer_product
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++ End outer_product ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++ Start sort_real ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine sort_real(num, idx)
    real(8), intent(inout) :: num(:)
    integer, intent(out), optional :: idx(:)
    integer :: Nele, i, j, Itemp
    real(8) :: temp
    logical :: Lidx
    Nele = size(num)

    Lidx = present(idx)
    if (Lidx) idx = [(i, i = 1, Nele)]

    do i = 1, Nele
      do j = i+1, Nele
        if (num(i) > num(j) ) then
          temp = num(i)
          num(i) = num(j)
          num(j) = temp

          if ( Lidx ) then
            Itemp  = idx(i)
            idx(i) = idx(j)
            idx(j) = Itemp
          end if
        end if
      end do
    end do
  end subroutine sort_real
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++ End!! sort_real ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module utility

