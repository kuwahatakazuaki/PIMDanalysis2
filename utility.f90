module utility
  use input_parameter,  only: data_beads, data_step, TNstep, graph_step
  implicit none
  real(8) :: pi = atan(1.0d0)*4.0d0
contains

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
    integer :: i, seedsize
    integer, allocatable :: seed(:)

    call random_seed(size=seedsize)
    allocate(seed(seedsize))
    seed(:) = 123
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

end module utility

