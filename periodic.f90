module mod_periodic
  use input_parameter, &
      only: jobtype, Natom, Nbeads, TNstep, label, save_beads, &
            hist_max1, hist_max2, hist_min1, hist_min2, &
            Nhist, &
            Ielement1, Ielement2, Felement1, Felement2, Noho, Lbox, label_oho, &
            r
!  use calc_parameter, only: data_beads, data_step
  use calc_histogram1D
  use utility
  implicit none
  private
  integer :: Uout
  real(8) :: Dhist
  public periodic

contains

  subroutine periodic

    select case(jobtype)
      case(81)
        call RDF1
      case(82)
        call RDF2
      case(83)
        call oho_distribution
      case(89)
        call rms_oho
      case default
        stop 'ERROR!!! wrong "Job type" option'
    end select
  end subroutine periodic

  subroutine rms_oho
    integer :: i, j, k, xyz
    integer :: Nh, No, Uout
    real(8) :: r3(3), d3
    real(8), allocatable :: rms(:,:)
    Nh = Felement1 - Ielement1 + 1
    No = Felement2 - Ielement2 + 1

    allocate(rms(Nh,TNstep))

    do k = 1, TNstep
      do i = 1, Nh
        do xyz = 1, 3
          r3(xyz) = sum( r(xyz,i,:,k)-r(xyz,i,:,1) ) / dble(Nbeads)
        end do
        d3 = dsqrt( sum(r3(:)*r3(:)) )
        rms(i,k) = d3
      end do
    end do

    open(newunit=Uout,file='rms.out')
      do k = 1, TNstep
        write(Uout,9999) rms(:,k)
      end do
    close(Uout)
    9999 format(32F10.5)
  end subroutine rms_oho

  subroutine oho_distribution
    integer :: i, j, k, Ioho, Uout
    real(8) :: r12(3), r23(3), r13(3), d12, d23, d13
    real(8), allocatable :: oho_bead(:,:,:), oho_step(:,:), oo_bead(:,:,:), oo_step(:,:), oh_bead(:,:,:)
    integer :: Nh, No
    Nh = Felement1 - Ielement1 + 1
    No = Felement2 - Ielement2 + 1

    allocate(oho_bead(Noho,Nbeads,TNstep))
    allocate(oho_step(Noho,TNstep))
    allocate(oo_bead(Noho,Nbeads,TNstep))
    allocate(oo_step(Noho,TNstep))
    allocate(oh_bead(Noho,Nbeads,TNstep))

    deallocate(data_beads,data_step)
    allocate(data_beads(Nbeads,TNstep*Noho))

    do k = 1, TNstep
      do j = 1, Nbeads
        do Ioho = 1, Noho
          r12(:) = r(:,label_oho(1,Ioho),j,k) - r(:,label_oho(2,Ioho),j,k)
          r23(:) = r(:,label_oho(2,Ioho),j,k) - r(:,label_oho(3,Ioho),j,k)
          r13(:) = r(:,label_oho(1,Ioho),j,k) - r(:,label_oho(3,Ioho),j,k)
          r12(:) = r12(:) - Lbox(:) * nint(r12(:)/Lbox(:))
          r23(:) = r23(:) - Lbox(:) * nint(r23(:)/Lbox(:))
          r13(:) = r13(:) - Lbox(:) * nint(r13(:)/Lbox(:))
          d12 = sqrt( sum(r12(:)*r12(:)) )
          d23 = sqrt( sum(r23(:)*r23(:)) )
          d13 = sqrt( sum(r13(:)*r13(:)) )
          oho_bead(Ioho,j,k) = d12 - d23
          oo_bead(Ioho,j,k) = d13
          oh_bead(Ioho,j,k) = min(d12,d23)
        end do
      end do
    end do

    open(newunit=Uout,file='step_oho.out')
      do k = 1, TNstep
        if (mod(k,graph_step)==0) then
          do i =1, Noho
            write(Uout,9999,advance='no') sum(oho_bead(i,:,k))/dble(Nbeads)
          end do
          write(Uout,*) ""
        end if
      end do
    close(Uout)


    do i = 1, Noho
      data_beads(:,TNstep*(i-1)+1:TNstep*i) = oho_bead(i,:,:)
    end do
    call calc_1Dhist(out_hist="hist1D_oho.out")

    hist_max1 = 0.0d0
    hist_min1 = 0.0d0
    do i = 1, Noho
      data_beads(:,TNstep*(i-1)+1:TNstep*i) = oo_bead(i,:,:)
    end do
    call calc_1Dhist(out_hist="hist1D_oo.out")

    hist_max1 = 0.0d0
    hist_min1 = 0.0d0
    do i = 1, Noho
      data_beads(:,TNstep*(i-1)+1:TNstep*i) = oh_bead(i,:,:)
    end do
    call calc_1Dhist(out_hist="hist1D_oh.out")

    if ( save_beads .eqv. .True. ) then
      open(newunit=Uout,file="oho.bin", form='unformatted', access='stream', status='replace')
        write(Uout) oho_bead(:,:,:)
      close(Uout)
      open(newunit=Uout,file="oo.bin", form='unformatted', access='stream', status='replace')
        write(Uout) oo_bead(:,:,:)
      close(Uout)
    end if


  9999 format(F10.5)
  end subroutine oho_distribution

! ++++++++++++++++++
! +++ Start RDF1 +++
! ++++++++++++++++++
  subroutine RDF1
    integer :: Nelement
    integer :: i, j, k, l, Ihist
    real(8) :: r12(3), minbox, d12, rho
    character(len=:), allocatable :: out_hist
    allocate(histogram(Nhist,2))
    histogram(:,:) = 0.0d0

    if ( trim(out_hist) == "0" ) then
      write(out_hist, '(a,I0,a,I0,a)') "rdf1_", Ielement1, "-", Felement1, ".out"
    end if

    minbox = minval(Lbox(:))
    Nelement = Felement1 - Ielement1 + 1
    rho = dble(Nelement*(Nelement-1)/2) / (Lbox(1)*Lbox(2)*Lbox(3))
    Dhist = minbox / dble(Nhist)
    histogram(:,:) = 0.0d0
    do Ihist = 1, Nhist
      histogram(Ihist,1) = Dhist * dble(Ihist)  ! not dble(Ihist-1)
    end do

    do i = 1, TNstep
      do j = 1, Nbeads
        do k = Ielement1, Felement1
          do l = k+1, Felement1
            r12(:) = r(:,k,j,i) - r(:,l,j,i)
            r12(:) = r12(:) - Lbox(:) * nint(r12(:)/Lbox(:))
            d12 = dsqrt( sum( r12(:)*r12(:) ) )
            do Ihist = 1, Nhist
              if ( d12 <= histogram(Ihist,1) ) then
                histogram(Ihist,2) = histogram(Ihist,2) + 1.0d0
                goto 100
              end if
            end do
            100 continue
          end do
        end do
      end do
    end do
    histogram(:,1) = histogram(:,1) - 0.5d0 * Dhist
    histogram(:,2) = histogram(:,2) / (4*pi*rho*Dhist*TNstep*Nbeads)
    do Ihist = 1, Nhist
      histogram(Ihist,2) = histogram(Ihist,2) / (histogram(Ihist,1)*histogram(Ihist,1))
    end do

    open(Uout, file=trim(out_hist), status='replace')
      do Ihist = 1, Nhist
        write(Uout,'(F13.6, E13.4)') histogram(Ihist,:)
      end do
    close(Uout)

  end subroutine RDF1
! ++++++++++++++++
! +++ End RDF1 +++
! ++++++++++++++++

! ++++++++++++++++++
! +++ Start RDF2 +++
! ++++++++++++++++++
  subroutine RDF2
    integer :: Nelement
    integer :: i, j, k, l, Ihist
    real(8) :: r12(3), minbox, d12, rho
    character(len=:), allocatable :: out_hist
    allocate(histogram(Nhist,2))

    if ( trim(out_hist) == "0" ) then
      write(out_hist, '(a,I0,a,I0,a,I0,a,I0,a)') & 
         "rdf2_", Ielement1, "-", Felement1, "_",Ielement2, "-",Felement2, ".out"
    end if

    minbox = minval(Lbox(:))
    Nelement = (Felement1 - Ielement1 + 1) + (Felement2 - Ielement2 + 1)
    rho = dble(Nelement) / (Lbox(1)*Lbox(2)*Lbox(3))
    Dhist = minbox / dble(Nhist)
    histogram(:,:) = 0.0d0
    do Ihist = 1, Nhist
      histogram(Ihist,1) = Dhist * dble(Ihist)  ! not dble(Ihist-1)
    end do

    do i = 1, TNstep
      do j = 1, Nbeads
        do k = Ielement1, Felement1
          do l = Ielement2, Felement2
            r12(:) = r(:,k,j,i) - r(:,l,j,i)
            r12(:) = r12(:) - Lbox(:) * nint(r12(:)/Lbox(:))
            d12 = dsqrt( sum( r12(:)*r12(:) ) )
            do Ihist = 1, Nhist
              if ( d12 <= histogram(Ihist,1) ) then
                histogram(Ihist,2) = histogram(Ihist,2) + 1.0d0
                goto 100
              end if
            end do
            100 continue
          end do
        end do
      end do
    end do
    histogram(:,1) = histogram(:,1) - 0.5d0 * Dhist
    histogram(:,2) = histogram(:,2) / (4*pi*rho*Dhist*TNstep*Nbeads)
    do Ihist = 1, Nhist
      histogram(Ihist,2) = histogram(Ihist,2) / (histogram(Ihist,1)*histogram(Ihist,1))
    end do

    open(Uout, file=trim(out_hist), status='replace')
      do Ihist = 1, Nhist
        write(Uout,'(F13.6, E13.4)') histogram(Ihist,:)
      end do
    close(Uout)

  end subroutine RDF2
! ++++++++++++++++
! +++ End RDF2 +++
! ++++++++++++++++
end module mod_periodic




