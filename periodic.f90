module mod_periodic
  use input_parameter, &
      only: jobtype, Natom, Nbeads, TNstep, label, save_beads, &
            hist_max1, hist_max2, hist_min1, hist_min2, Nhist, lattice, &
            Ielement1, Ielement2, Felement1, Felement2, Noho, Lbox, label_oho, &
            r, data_beads, data_step, graph_step, atom1, atom2, atom3, atom4, &
            Ntetra, Itetra, Ndiv
  use calc_histogram1D, only: calc_1Dhist
  use utility, only: get_volume, pi, get_inv_mat! , save_cube
  implicit none
  private
  integer :: Uout
  real(8) :: Dhist, lat_inv(3,3)
  real(8), allocatable :: s(:,:,:,:), hist(:,:)
  public periodic

contains

  subroutine periodic
    integer :: i
    print *, "***** START periodic condition *****"
    do i = 1, 3
      print '(a,3F13.6)', '    lattice = ', lattice(i,:)
    end do

    select case(jobtype)
      case(81)
        call RDF1
      case(82)
        call RDF2
      case(83)
        call bond_perio
      case(84)
        call bond_diff_perio
      case(85)
        call RMSDatoms
      !case(88)
      !  call Tetrahedron
      case(89)
        call oho_distribution
      !case(89)
      !  call rms_oho
      case default
        stop 'ERROR!!! wrong "Job type" option'
    end select
    print *, "***** END periodic condition *****"
  end subroutine periodic

  subroutine bond_diff_perio
    integer :: Nelement
    integer :: i, j, k, l, Ihist
    real(8) :: r12(3), r34(3), s12(3), s34(3), d12, d34
    character(len=128) :: out_hist
    allocate(hist(Nhist,2), source=0.0d0)

    write(out_hist, '(a,I0,a,I0,a,I0,a,I0,a)') &
               "hist_", atom1, "_", atom2,"-",atom3,"_",atom4, ".out"

    call get_inv_mat(lattice,lat_inv,3)

    allocate(s(3,Natom,Nbeads,TNstep))
    do k = 1, TNstep
      do j = 1, Nbeads
        do i = 1, Natom
          s(:,i,j,k) = matmul(r(:,i,j,k), lat_inv(:,:))
        end do
      end do
    end do

    do i = 1, TNstep
      do j = 1, Nbeads
        s12(:) = s(:,atom1,j,i) - s(:,atom2,j,i)
        s12(:) = s12(:) - nint(s12(:))
        r12(:) = matmul(s12(:),lattice(:,:))
        d12 = norm2(r12(:))

        s34(:) = s(:,atom3,j,i) - s(:,atom4,j,i)
        s34(:) = s34(:) - nint(s34(:))
        r34(:) = matmul(s34(:),lattice(:,:))
        d34 = norm2(r34(:))
        data_beads(j,i) = d12 - d34
      end do
    end do
    call calc_1Dhist(out_hist=trim(out_hist))
stop 'Not Update'
  end subroutine bond_diff_perio

  subroutine bond_perio
    integer :: Nelement
    integer :: i, j, k, l, Ihist
    real(8) :: r12(3), s12(3), d12
    character(len=128) :: out_hist
    allocate(hist(Nhist,2), source=0.0d0)

    write(out_hist, '(a,I0,a,I0,a)') "hist_", atom1, "-", atom2, ".out"

    call get_inv_mat(lattice,lat_inv,3)

    allocate(s(3,Natom,Nbeads,TNstep))
    do k = 1, TNstep
      do j = 1, Nbeads
        do i = 1, Natom
          s(:,i,j,k) = matmul(r(:,i,j,k), lat_inv(:,:))
        end do
      end do
    end do

    do i = 1, TNstep
      do j = 1, Nbeads
        s12(:) = s(:,atom1,j,i) - s(:,atom2,j,i)
        s12(:) = s12(:) - nint(s12(:))
        r12(:) = matmul(s12(:),lattice(:,:))
        d12 = norm2(r12(:))
        data_beads(j,i) = d12
      end do
    end do
    call calc_1Dhist(out_hist=trim(out_hist))
  end subroutine bond_perio

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++ Start save_cube ++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine save_cube_sub(rcub,Iatoms)
    integer, intent(in) :: Iatoms(:)
    ! character(len=*), optional :: Fout
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

    allocate(coun(3,5,Nbeads,TNstep))
    ! allocate(coun(3,Ncube,Nbeads,TNstep))
    do k = 1, TNstep
      do j = 1, Nbeads
        do i = 1, 5
          coun(:,i,j,k) = int( ( rcub(:,i,j,k)-Lmin(:) ) / dL(:) ) + 1
          !coun(:,i,j,k) = int( ( rcub(:,Iatoms(i),j,k)-Lmin(:) ) / dL(:) ) + 1
        end do
      end do
    end do

    block
      integer :: c1,c2,c3
      grid(:,:,:) = 0.0d0
      do k = 1, TNstep
        do j = 1, Nbeads
          c1 = coun(1,1,j,k); c2 = coun(2,1,j,k); c3 = coun(3,1,j,k)
          grid(c1,c2,c3) = grid(c1,c2,c3) + 1.0d0
          do i = 2, 5
            c1 = coun(1,i,j,k); c2 = coun(2,i,j,k); c3 = coun(3,i,j,k)
            grid(c1,c2,c3) = grid(c1,c2,c3) - 1.0d0
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
        write(Uout,9999) -1, dble(i), &  ! Only for Tetrahedra O-H4
        !write(Uout,9999) atom2num(trim(label(i))), dble(i), &
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

  end subroutine save_cube_sub
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +++++ End save_cube ++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



! ++++++++++++++++++++++
! +++++ Start RMSD +++++
! ++++++++++++++++++++++
  subroutine RMSDatoms
    real(8), allocatable :: rc(:,:,:), rmsd(:,:)
    real(8) :: dis2, rij(3)
    integer :: i, j, k, Uout
    allocate(rc(3,Natom,TNstep), source=0.0d0)
    allocate(rmsd(Natom,TNstep))
    do j = 1, Nbeads
      rc(:,:,:) = rc(:,:,:) + r(:,:,j,:)
    end do
    rc(:,:,:) = rc(:,:,:) / dble(Nbeads)

    do k = 1, TNstep
      do i = 1, Natom
        rij(:) = rc(:,i,k) - rc(:,i,1)
        rmsd(i,k) = dot_product(rij(:),rij(:))
      end do
    end do
    rmsd(:,:) = dsqrt(rmsd(:,:))

    open(newunit=Uout,file='rmsd.out')
      do k = 1, TNstep
        if (mod(k,graph_step) == 0) then
          write(Uout,9999) rmsd(atom1:atom2,k)
        end if
      end do
    close(Uout)
  9999 format(48F8.4)
  end subroutine RMSDatoms
! ++++++++++++++++++++
! +++++ End RMSD +++++
! ++++++++++++++++++++


! ++++++++++++++++++++++
! +++++ Start RDF1 +++++
! ++++++++++++++++++++++
  subroutine RDF1
    integer :: Nelement
    integer :: i, j, k, l, Ihist
    real(8) :: r12(3), s12(3), minedge, rho, d12
    character(len=128) :: out_hist
    allocate(hist(Nhist,2), source=0.0d0)
    !hist(:,:) = 0.0d0

    write(out_hist, '(a,I0,a,I0,a)') "rdf1_", Ielement1, "-", Felement1, ".out"

    minedge = get_min_edge(lattice(:,:))
    Dhist = minedge / dble(Nhist)
    Nelement = Felement1 - Ielement1 + 1
    rho = dble(Nelement*(Nelement-1)/2) / (get_volume(lattice(:,:)))
    !hist(:,:) = 0.0d0

    print '(a,I0,"-"I0)', '    Radial distribution of ',Ielement1,Felement1
    print '(a,I0)',       '    Nelement =  ', Nelement
    print '(a,F13.6)',    '    minedge  =  ', minedge

    do Ihist = 1, Nhist
      hist(Ihist,1) = Dhist * dble(Ihist)  ! not dble(Ihist-1)
    end do

    call get_inv_mat(lattice,lat_inv,3)

    allocate(s(3,Natom,Nbeads,TNstep))
    do k = 1, TNstep
      do j = 1, Nbeads
        do i = 1, Natom
          s(:,i,j,k) = matmul(r(:,i,j,k), lat_inv(:,:))
        end do
      end do
    end do

    do i = 1, TNstep
      do j = 1, Nbeads
        do k = Ielement1, Felement1
          do l = k+1, Felement1
            s12(:) = s(:,k,j,i) - s(:,l,j,i)
            s12(:) = s12(:) - nint(s12(:))
            r12(:) = matmul(s12(:),lattice(:,:))
            d12 = dsqrt( sum( r12(:)*r12(:) ) )

            Ihist = int( (d12-0.0d0)/Dhist ) + 1
            if ( Ihist <= Nhist ) then
              hist(Ihist,2) = hist(Ihist,2) + 1.0d0
            end if

          end do
        end do
      end do
    end do
    hist(:,1) = hist(:,1) - 0.5d0 * Dhist
    hist(:,2) = hist(:,2) / (4*pi*rho*Dhist*TNstep*Nbeads)
    do Ihist = 1, Nhist
      hist(Ihist,2) = hist(Ihist,2) / (hist(Ihist,1)*hist(Ihist,1))
    end do

    open(newunit=Uout, file=trim(out_hist), status='replace')
      write(Uout,'(" # Maximum hist =", F13.6)') maxval(hist(:,2))
      write(Uout,'(" # Max loc step =", I8)')    maxloc(hist(:,2))
      write(Uout,'(" # Max location =", F13.6)') hist(maxloc(hist(:,2)),1)
      do Ihist = 1, Nhist
        write(Uout,'(F13.6, E13.4)') hist(Ihist,:)
      end do
    close(Uout)

  end subroutine RDF1
! ++++++++++++++++++++
! +++++ End RDF1 +++++
! ++++++++++++++++++++

! ++++++++++++++++++++++
! +++++ Start RDF2 +++++
! ++++++++++++++++++++++
  subroutine RDF2
    integer :: Nelement
    integer :: i, j, k, l, Ihist
    real(8) :: r12(3), s12(3), minedge, d12, rho
    character(len=128) :: out_hist
    allocate(hist(Nhist,2))

    write(out_hist, '(a,I0,a,I0,a,I0,a,I0,a)') & 
       "rdf2_", Ielement1, "-", Felement1, "_",Ielement2, "-",Felement2, ".out"

    minedge = get_min_edge(lattice(:,:))
    Dhist = minedge / dble(Nhist)
    Nelement = (Felement1 - Ielement1 + 1) * (Felement2 - Ielement2 + 1)
    rho = dble(Nelement) / get_volume(lattice(:,:))
    hist(:,:) = 0.0d0

    print '(a,I0,"-",I0," to ",I0,"-",I0)', '    Radial distribution of ',Ielement1,Felement1,Ielement2,Felement2
    print '(a,I10)',      '    Nelement =  ', Nelement
    print '(a,F10.6)',    '    minedge  =  ', minedge


    do Ihist = 1, Nhist
      hist(Ihist,1) = Dhist * dble(Ihist)  ! not dble(Ihist-1)
    end do

    call get_inv_mat(lattice,lat_inv,3)

    allocate(s(3,Natom,Nbeads,TNstep))
    do k = 1, TNstep
      do j = 1, Nbeads
        do i = 1, Natom
          s(:,i,j,k) = matmul(r(:,i,j,k), lat_inv(:,:))
        end do
      end do
    end do

    do i = 1, TNstep
      do j = 1, Nbeads
        do k = Ielement1, Felement1
          do l = Ielement2, Felement2
            s12(:) = s(:,k,j,i) - s(:,l,j,i)
            s12(:) = s12(:) - nint(s12(:))
            r12(:) = matmul(s12(:),lattice(:,:))
            d12 = dsqrt( sum( r12(:)*r12(:) ) )

            Ihist = int( (d12-0.0d0)/Dhist ) + 1
            if ( Ihist <= Nhist ) then
              hist(Ihist,2) = hist(Ihist,2) + 1.0d0
            end if

          end do
        end do
      end do
    end do
    hist(:,1) = hist(:,1) - 0.5d0 * Dhist
    hist(:,2) = hist(:,2) / (4*pi*rho*Dhist*TNstep*Nbeads)
    do Ihist = 1, Nhist
      hist(Ihist,2) = hist(Ihist,2) / (hist(Ihist,1)*hist(Ihist,1))
    end do

    open(newunit=Uout, file=trim(out_hist), status='replace')
      write(Uout,'(" # Maximum hist =", F13.6)') maxval(hist(:,2))
      write(Uout,'(" # Max loc step =", I8)')    maxloc(hist(:,2))
      write(Uout,'(" # Max location =", F13.6)') hist(maxloc(hist(:,2)),1)
      do Ihist = 1, Nhist
        write(Uout,'(F13.6, E13.4)') hist(Ihist,:)
      end do
    close(Uout)

  end subroutine RDF2
! ++++++++++++++++++++
! +++++ End RDF2 +++++
! ++++++++++++++++++++

! ++++++++++++++++++++++++++++++++++
! +++++ Start oho_distribution +++++
! ++++++++++++++++++++++++++++++++++
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
! ++++++++++++++++++++++++++++++++
! +++++ End oho_distribution +++++
! ++++++++++++++++++++++++++++++++

  function get_min_edge(vec) result(mini)
    real(8) :: vec(3,3), mini
    real(8) :: edge(3)
    integer :: i
    do i = 1, 3
      edge(i) = dsqrt(dot_product(vec(i,:), vec(i,:)))
    end do
    mini = minval(edge)
  end function get_min_edge

end module mod_periodic

!! +++++++++++++++++++++++++++++
!! +++++ Start Tetrahedron +++++
!! +++++++++++++++++++++++++++++
!  subroutine Tetrahedron
!    integer :: i, j, k, xyz, Utem
!    real(8), allocatable :: rt(:,:,:,:,:), st(:,:,:,:)
!    real(8), allocatable :: rcub(:,:,:,:)
!    real(8) :: rc(3), Lbox(3)
!    integer :: Iatoms(4) = [2,3,4,5]
!
!    do i = 1, 3
!      Lbox(i) = lattice(i,i)
!    end do
!
!    !call get_inv_mat(lattice,lat_inv,3)
!
!    allocate(rt(3,5,Ntetra,Nbeads,TNstep))
!    allocate(st(3,5,Ntetra,Nbeads))
!
!    do k = 1, TNstep
!      do xyz = 1, 3
!        rc(xyz) = sum(r(xyz,:,:,k)) / dble(Natom*Nbeads)
!        r(xyz,:,:,k) = r(xyz,:,:,k) - rc(xyz)
!      end do
!    end do
!
!    do i = 1, Ntetra
!      do j = 1, 5
!        rt(:,j,i,:,:) = r(:,Itetra(i,j),:,:)
!      end do
!      do xyz = 1, 3
!        rc(xyz) = sum(r(xyz,Itetra(i,1),:,:)) / dble(Nbeads*TNstep)
!        rt(xyz,:,i,:,:) = rt(xyz,:,i,:,:) - rc(xyz)
!      end do
!    end do
!
!    do k = 1, TNstep
!      do i = 1, Ntetra
!        do xyz = 1, 3
!          st(xyz,:,i,:) = rt(xyz,:,i,:,k) / Lbox(xyz)
!        end do
!        st(:,:,:,:) = st(:,:,:,:) - nint(st(:,:,:,:))
!        do xyz = 1, 3
!          rt(xyz,:,i,:,k) = st(xyz,:,i,:) * Lbox(xyz)
!        end do
!      end do
!    end do
!
!    allocate(rcub(3,5,Nbeads,TNstep*Ntetra))
!    do i = 1, Ntetra
!      rcub(:,:,:,TNstep*(i-1)+1:TNstep*i) = rt(:,:,i,:,:)
!    end do
!
!    call save_cube_sub(rcub,Iatoms)
!  end subroutine Tetrahedron
!! +++++++++++++++++++++++++++
!! +++++ End Tetrahedron +++++
!! +++++++++++++++++++++++++++

