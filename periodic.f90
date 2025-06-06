module mod_periodic
  use input_parameter, &
      only: jobtype, Natom, Nbeads, TNstep, label, save_beads, &
            hist_max1, hist_max2, hist_min1, hist_min2, Nhist, lattice, &
            Ielement1, Ielement2, Felement1, Felement2, Noho, Lbox, label_oho, &
            r, data_beads, data_step, graph_step, atom1, atom2, atom3, atom4, &
            Ntetra, Itetra, Ndiv, Natom_peri, FNameBinary1
  use calc_histogram1D, only: calc_1Dhist
  use utility, only: get_volume, pi, get_inv_mat, real_max, sort !, sort_real
  implicit none
  private
  integer :: Uout
  real(8) :: Dhist, lat_inv(3,3)
  real(8), allocatable :: s(:,:,:,:), hist(:,:)
  integer :: Nelement1, Nelement2
  public periodic

contains

  subroutine periodic
    integer :: i, j, k

    Nelement1 = Felement1-Ielement1+1
    Nelement2 = Felement2-Ielement2+1

    print *, "***** START periodic condition *****"
    do i = 1, 3
      print '(a,3F13.6)', '    lattice = ', lattice(i,:)
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
      case(86)
        call Minimum_bond
      case(87)
        call Near_str1
      case(88)
        call Near_str2
      case(89)
        call Near_atoms1
      case(98)
        call Fealloy_near4
      case(99)
        call oho_distribution
      case default
        stop 'ERROR!!! wrong "Job type" option'
    end select
    print *, "***** END periodic condition *****"
  end subroutine periodic
! +++++++++++++++++++++++++++++++
! +++++ Start Fealloy_near4 +++++
! +++++++++++++++++++++++++++++++
  subroutine Fealloy_near4
    use utility, only: count_letter
    integer, parameter :: Nsite = 2592
    integer :: Iatom, Istep, Ibead, i, j, k, Idis, Uinp, ios, Isite
    !integer :: near_idx(Natom_peri),   temp_idx(Nelement2)
    integer :: near_idx(Natom_peri+1),   temp_idx(Nelement2)
    real(8) :: near_str(3,Natom_peri), temp_dis(Nelement2), temp_str(3,Nelement2)
    real(8) :: s12(3), r12(3), dis2
    integer :: site_dat(4,2,Nsite)
    logical :: Lmatch
    character(len=*), parameter :: Fsite="site.dat", Fout = 'near_atoms.xyz'
    character :: Cdummy

    open(newunit=Uinp,file=Fsite,status='old',iostat=ios)
      if (ios /= 0) then
        print *, "ERROR!! opening file : ", Fout
        stop
      end if
      do i = 1, Nsite
        read(Uinp,*) site_dat(:,1,i), Cdummy, site_dat(:,2,i)
      end do
    close(Uinp)

    open(newunit=Uout,file=Fout,status='replace')
    !LoopStep : do Istep = 1, 10
    LoopStep : do Istep = 1, TNstep
    LoopAtom : do Iatom = Ielement1, Felement1
      write(Uout,*) (Natom_peri+1)*Nbeads
      write(Uout,*) Istep, Iatom
      do Ibead = 1, Nbeads
        temp_dis(:) = real_max

        Idis = 1
        do i = Ielement2, Felement2
          s12(:) = s(:,i,Ibead,Istep) - s(:,Iatom,Ibead,Istep)
          s12(:) = s12(:) - anint(s12(:))
          r12(:) = matmul(s12(:),lattice(:,:))
          dis2 = dot_product(r12(:),r12(:))
          temp_dis(Idis) = dis2
          temp_str(:,Idis) = r12(:)
          Idis = Idis + 1
        end do
        call sort(temp_dis,temp_idx)
        near_idx(1:Natom_peri+1) = temp_idx(1:Natom_peri+1)
        !near_idx(1:Natom_peri) = temp_idx(1:Natom_peri)

        call compare_array(Isite,Lmatch)
        if ( Lmatch ) then
          do i = 1, Natom_peri
            near_str(:,i) = temp_str(:,site_dat(i,2,Isite))
          end do
        else
          do i = 1, Natom_peri
            near_str(:,i) = temp_str(:,near_idx(i))
print *, near_str(:,i)
          end do
stop 'HERE'
        end if

        write(Uout,*) label(Iatom), 0.0d0, 0.0d0, 0.0d0
        do i = 1, Natom_peri
          write(Uout,*) label(near_idx(i)+Felement1), near_str(:,i), site_dat(i,2,Isite)
          !write(Uout,*) label(near_idx(i)+Felement1), near_str(:,i), near_idx(i)
        end do
      end do

    end do LoopAtom
    end do LoopStep
    close(Uout)

    !open(newunit=Uout,file='near_label.dat',status='replace')
    !  do Iatom = atom1, atom2
    !    do Istep = 1, TNstep
    !      do Ibead = 1, Nbeads
    !        !write(Uout,*) Alllabel(:,Ibead,Istep,Iatom), count_letter(Alllabel(:,Ibead,Istep,Iatom),'Fe')
    !        write(Uout,*) Istep, Ibead, count_letter(Alllabel(:,Ibead,Istep,Iatom),'Fe')
    !      end do
    !    end do
    !  end do
    !close(Uout)
contains
    subroutine compare_array(Isite, Lmatch)
      integer :: temp_arry(Natom_peri+1)
      !integer :: temp_arry(Natom_peri)
      integer, intent(out) :: Isite
      logical, intent(out) :: Lmatch
      integer :: i, j
      logical :: Lenddo
      temp_arry(:) = near_idx(:)
      call sort(temp_arry(1:4))
      Lenddo = .True.
      do Isite = 1, Nsite
        if ( all(temp_arry(1:4) == site_dat(:,1,Isite)) ) then
          Lenddo = .False.
          exit
        end if
      end do
!print *, Isite, ":", near_idx(1:4)
!print *, site_dat(:,1,Isite)
!print *, site_dat(:,2,Isite)

      if ( Lenddo ) then
        Lenddo = .True.
!print *,Istep, Iatom, Ibead, ":", temp_arry(:)
        temp_arry(:) = near_idx(:)
        temp_arry(4) = near_idx(5)
        call sort(temp_arry(1:4))
        do Isite = 1, Nsite
          if ( all(temp_arry(1:4) == site_dat(:,1,Isite)) ) then
            !print *, "match!! ",i, ":", temp_arry(:)
            Lenddo = .False.
            exit
          end if
        end do
      end if
      if ( Lenddo ) then
        Lmatch = .True.
      end if
    end subroutine compare_array
  end subroutine Fealloy_near4
! +++++++++++++++++++++++++++++++
! +++++ End!! Fealloy_near4 +++++
! +++++++++++++++++++++++++++++++

! +++++++++++++++++++++++++++
! +++++ Start Near_str1 +++++
! +++++++++++++++++++++++++++
  subroutine Near_atoms1
    use utility, only: count_letter
    integer :: Iatom, Istep, Ibead, i, j, k!, Ulab
    integer :: near_idx(Natom_peri), temp_idx(Natom)
    real(8) :: near_str(3,Natom_peri), temp_dis(Natom), temp_str(3,Natom)
    real(8) :: s12(3), r12(3), dis2
    character(:), allocatable :: Fout!, Flab
    character(len=2) :: Alllabel(Natom_peri,Nbeads,TNstep,atom2-atom1+1)
    !integer :: count_letter

    Fout = 'near_atoms.xyz'
    open(newunit=Uout,file=Fout,status='replace')
    LoopStep : do Istep = 1, TNstep
    LoopAtom : do Iatom = atom1, atom2
      write(Uout,*) (Natom_peri+1)*Nbeads
      write(Uout,*) Istep, Iatom
      do Ibead = 1, Nbeads
        temp_dis(:) = real_max

        do i = 1, Natom
          if ( i == Iatom ) cycle
          s12(:) = s(:,i,Ibead,Istep) - s(:,Iatom,Ibead,Istep)
          s12(:) = s12(:) - anint(s12(:))
          r12(:) = matmul(s12(:),lattice(:,:))
          dis2 = dot_product(r12(:),r12(:))
          temp_dis(i) = dis2
          temp_str(:,i) = r12(:)
        end do
        call sort(temp_dis,temp_idx)
        near_idx(1:Natom_peri) = temp_idx(1:Natom_peri)
        do j = 1, Natom_peri
          near_str(:,j) = temp_str(:,near_idx(j))
        end do

        write(Uout,*) label(Iatom), 0.0d0, 0.0d0, 0.0d0
        do i = 1, Natom_peri
          write(Uout,*) label(near_idx(i)), near_str(:,i)
        end do
        do i = 1, Natom_peri
          Alllabel(i,Ibead,Istep,Iatom) = label(near_idx(i))
        end do
      end do

    end do LoopAtom
    end do LoopStep
    close(Uout)

    open(newunit=Uout,file='near_label.dat',status='replace')
      do Iatom = atom1, atom2
        do Istep = 1, TNstep
          do Ibead = 1, Nbeads
            !write(Uout,*) Alllabel(:,Ibead,Istep,Iatom), count_letter(Alllabel(:,Ibead,Istep,Iatom),'Fe')
            write(Uout,*) Istep, Ibead, count_letter(Alllabel(:,Ibead,Istep,Iatom),'Fe')
          end do
        end do
      end do
    close(Uout)
  end subroutine Near_atoms1
! +++++++++++++++++++++++++++
! +++++ End!! Near_str1 +++++
! +++++++++++++++++++++++++++


! +++++++++++++++++++++++++++
! +++++ Start Near_str2 +++++
! +++++++++++++++++++++++++++
  subroutine Near_str2
    integer, parameter :: max_atom = 20
    integer :: Istep, Ibead, i, j, k
    integer :: Nnear, near_idx(max_atom)
    real(8) :: si1(3), ri1(3), si2(3), ri2(3), dis1, dis2
    real(8) :: diff_r(3), diff_s(3)
    real(8) :: near_str(3,max_atom)
    real(8), parameter :: cut_dis = 3.3d0
    character(:), allocatable :: Fout

    Fout = 'near_str2.xyz'
    open(newunit=Uout,file=Fout,status='replace')
    do Istep = 1, TNstep
      do Ibead = 1, Nbeads
        Nnear = 0
        do i = 1, Natom
          if ( i == atom1 .or. i == atom2) cycle
          diff_r(:) = r(:,atom2,Ibead,Istep) - r(:,atom1,Ibead,Istep)
          diff_s(:) = matmul(diff_r(:),lat_inv(:,:))
          diff_s(:) = diff_s(:) - anint(diff_s(:))
          diff_r(:) = matmul(diff_s(:),lattice(:,:))

          si1(:) = s(:,i,Ibead,Istep) - s(:,atom1,Ibead,Istep)
          si2(:) = s(:,i,Ibead,Istep) - s(:,atom2,Ibead,Istep)

          si1(:) = si1(:) - anint(si1(:))
          si2(:) = si2(:) - anint(si2(:))

          ri1(:) = matmul(si1(:),lattice(:,:))
          ri2(:) = matmul(si2(:),lattice(:,:))

          dis1 = dot_product(ri1(:),ri1(:))
          dis2 = dot_product(ri2(:),ri2(:))

          if ( dis1 <= cut_dis**2 ) then
            Nnear = Nnear + 1
            if ( Nnear > max_atom ) stop 'ERROR!! Nnear exceed max_atom'
            near_idx(Nnear) = i
            near_str(:,Nnear) = ri1(:)
          else if ( dis2 <= cut_dis**2  ) then
            Nnear = Nnear + 1
            if ( Nnear > max_atom ) stop 'ERROR!! Nnear exceed max_atom'
            near_idx(Nnear) = i
            near_str(:,Nnear) = ri2(:) + diff_r(:)
          end if
        end do
      end do

      if ( 4000 <= Istep .and. Istep <= 4400 ) then
      write(Uout,*) Nnear + 2
      write(Uout,*) Istep
      write(Uout,*) label(atom1), 0.0d0, 0.0d0, 0.0d0
      write(Uout,*) label(atom2), diff_r(:)
      do i = 1, Nnear
        write(Uout,*)  label(near_idx(i)), near_str(:,i)
      end do
      end if
    end do
    close(Uout)
  end subroutine Near_str2
! +++++++++++++++++++++++++++
! +++++ End!! Near_str2 +++++
! +++++++++++++++++++++++++++

! +++++++++++++++++++++++++++
! +++++ Start Near_str1 +++++
! +++++++++++++++++++++++++++
  subroutine Near_str1
    integer, parameter :: max_atom = 20
    integer :: Istep, Ibead, i, j, k
    integer :: Nnear, near_idx(max_atom)
    real(8) :: s12(3), r12(3), dis2
    real(8) :: near_str(3,max_atom)
    real(8), parameter :: cut_dis = 4.0d0
    character(:), allocatable :: Fout

    Fout = 'near_str1.xyz'
    open(newunit=Uout,file=Fout,status='replace')
    do Istep = 1, TNstep
      do Ibead = 1, Nbeads
        Nnear = 0
        do i = 1, Natom
          if ( i == atom1 ) cycle
          s12(:) = s(:,i,Ibead,Istep) - s(:,atom1,Ibead,Istep)
          s12(:) = s12(:) - anint(s12(:))
          r12(:) = matmul(s12(:),lattice(:,:))
          dis2 = dot_product(r12(:),r12(:))
          if ( dis2 <= cut_dis**2 ) then
            Nnear = Nnear + 1
            if ( Nnear > max_atom ) stop 'ERROR!! Nnear exceed max_atom'
            near_idx(Nnear) = i
            near_str(:,Nnear) = r12(:)
          end if
        end do
      end do

      write(Uout,*) Nnear + 1
      write(Uout,*) Istep
      write(Uout,*) label(atom1), 0.0d0, 0.0d0, 0.0d0
      do i = 1, Nnear
        write(Uout,*) label(near_idx(i)), near_str(:,i)
      end do
    end do
    close(Uout)
  end subroutine Near_str1
! +++++++++++++++++++++++++++
! +++++ End!! Near_str1 +++++
! +++++++++++++++++++++++++++

  subroutine Minimum_bond
    integer :: Istep, i, j, k
    integer :: Iloc(3)
    real(8) :: s12(3), r12(3), d12, dis2(Natom,Natom,Nbeads)
    character(:), allocatable :: Fout

    Fout='step_mini_bond.out'

    open(newunit=Uout,file=Fout,status='replace')
    write(Uout,*) "# step, bond, atom1, atom2, Beads"

    !dis2(:,:,:) = 1.0d+100
    dis2(:,:,:) = real_max
    do Istep = 1, TNstep
      do k = 1, Nbeads
        do i = 1, Natom-1
          do j = i+1, Natom
            s12(:) = s(:,i,k,Istep) - s(:,j,k,Istep)
            s12(:) = s12(:) - anint(s12(:))
            r12(:) = matmul(s12(:),lattice(:,:))
            dis2(i,j,k) = dot_product(r12(:),r12(:))
          end do
        end do
      end do

      Iloc = minloc(dis2)
      write(Uout,'(I5,F10.5,3X,A,I0,"-",A,I0)', advance='no') &
          Istep, dsqrt(minval(dis2)), trim(label(Iloc(1))),Iloc(1),  trim(label(Iloc(2))),Iloc(2)! , Iloc(3)

      !dis2(Iloc(1), Iloc(2), :) = real_max ??
      dis2(Iloc(1), Iloc(2), Iloc(3)) = real_max
      Iloc = minloc(dis2)
      write(Uout,'(F10.5,3X,A,I0,"-",A,I0)') &
           dsqrt(minval(dis2)), trim(label(Iloc(1))),Iloc(1),  trim(label(Iloc(2))),Iloc(2)! , Iloc(3)
    end do
    close(Uout)
  end subroutine Minimum_bond

  subroutine bond_diff_perio
    integer :: i, j, k, l, Ihist
    real(8) :: r12(3), r34(3), s12(3), s34(3), d12, d34
    character(len=128) :: out_hist
    allocate(hist(Nhist,2), source=0.0d0)

    write(out_hist, '(a,I0,a,I0,a,I0,a,I0,a)') &
               "hist_", atom1, "_", atom2,"-",atom3,"_",atom4, ".out"

    do i = 1, TNstep
      do j = 1, Nbeads
        s12(:) = s(:,atom1,j,i) - s(:,atom2,j,i)
        s12(:) = s12(:) - anint(s12(:))
        r12(:) = matmul(s12(:),lattice(:,:))
        d12 = norm2(r12(:))

        s34(:) = s(:,atom3,j,i) - s(:,atom4,j,i)
        s34(:) = s34(:) - anint(s34(:))
        r34(:) = matmul(s34(:),lattice(:,:))
        d34 = norm2(r34(:))
        data_beads(j,i) = d12 - d34
      end do
    end do
    call calc_1Dhist(out_hist=trim(out_hist))
stop 'Not Update'
  end subroutine bond_diff_perio

  subroutine bond_perio
    integer :: i, j, k, l, Ihist
    real(8) :: r12(3), s12(3), d12
    character(len=128) :: out_hist
    allocate(hist(Nhist,2), source=0.0d0)

    write(out_hist, '(a,I0,a,I0,a)') "hist_", atom1, "-", atom2, ".out"

    do i = 1, TNstep
      do j = 1, Nbeads
        s12(:) = s(:,atom1,j,i) - s(:,atom2,j,i)
        s12(:) = s12(:) - anint(s12(:))
        r12(:) = matmul(s12(:),lattice(:,:))
        d12 = norm2(r12(:))
        data_beads(j,i) = d12
      end do
    end do
    call calc_1Dhist(out_hist=trim(out_hist))
  end subroutine bond_perio

! ++++++++++++++++++++++
! +++++ Start RMSD +++++
! ++++++++++++++++++++++
  subroutine RMSDatoms
    real(8), allocatable :: rc(:,:,:), rave(:,:)
    real(8), allocatable :: rmsd_atom(:,:), rmsd(:), msd2(:)
    real(8) :: dis2, rij(3)
    integer :: i, j, k, Uout, Nrmsd
    character(len=32)  :: fmt1
    character(len=256) :: fmt2

    Nrmsd = atom2 - atom1 + 1
    allocate(rc(3,Natom,TNstep), source=0.0d0)
    allocate(rave(3,TNstep), source=0.0d0)
    allocate(rmsd_atom(Nrmsd,TNstep))
    allocate(rmsd(TNstep))
    allocate(msd2(TNstep))

    ! +++ Removing the center of mass +++
    do j = 1, Nbeads
      rc(:,:,:) = rc(:,:,:) + r(:,:,j,:)
    end do
    rc(:,:,:) = rc(:,:,:) / dble(Nbeads)
    do i = 1, Natom
      rave(:,:) = rave(:,:) + rc(:,i,:)
    end do
    rave(:,:) = rave(:,:) / dble(Natom)
    do i = 1, Natom
      rc(:,i,:) = rc(:,i,:) - rave(:,:)
    end do
    ! +++ Removing the center of mass +++

    do k = 1, TNstep
      do i = atom1, atom2
        rij(:) = rc(:,i,k) - rc(:,i,1)
        rmsd_atom(i,k) = dot_product(rij(:),rij(:))
      end do
    end do
    do k = 1, TNstep
      rmsd(k) = sum(rmsd_atom(:,k))
    end do
    rmsd_atom(:,:) = dsqrt(rmsd_atom(:,:))
    msd2(:) = rmsd(:)/dble(Nrmsd)
    rmsd(:) = dsqrt(msd2(:))
    !rmsd(:) = dsqrt(rmsd(:)/dble(Nrmsd))

    write(fmt1,'(A,I0,A)') '(',Nrmsd+2,'G10.3)'
    fmt2 = "# Ave       d2"
    do i = atom1, atom2
      write(fmt2,'(a, "       ",a,I0)') trim(fmt2), label(i), i
    end do

    open(newunit=Uout,file='rmsd.out')
      write(Uout,'(a)') trim(fmt2)
      do k = 1, TNstep
        !if (mod(k,graph_step) == 0) then
          write(Uout,fmt1) rmsd(k), msd2(k), rmsd_atom(:,k)
        !end if
      end do
    close(Uout)
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

    write(out_hist, '(a,I0,a,I0,a)') "rdf1_", Ielement1, "-", Felement1, ".out"

    minedge = get_min_edge(lattice(:,:))
    Dhist = minedge / dble(Nhist)
    Nelement = Felement1 - Ielement1 + 1
    rho = dble(Nelement*(Nelement-1)/2) / (get_volume(lattice(:,:)))

    print '(a,I0,"-"I0)', '    Radial distribution of ',Ielement1,Felement1
    print '(a,I0)',       '    Nelement =  ', Nelement
    print '(a,F13.6)',    '    minedge  =  ', minedge

    do Ihist = 1, Nhist
      hist(Ihist,1) = Dhist * dble(Ihist)  ! not dble(Ihist-1)
    end do

    do i = 1, TNstep
      do j = 1, Nbeads
        do k = Ielement1, Felement1
          do l = k+1, Felement1
            s12(:) = s(:,k,j,i) - s(:,l,j,i)
            s12(:) = s12(:) - anint(s12(:))
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

    do i = 1, TNstep
      do j = 1, Nbeads
        do k = Ielement1, Felement1
          do l = Ielement2, Felement2
            s12(:) = s(:,k,j,i) - s(:,l,j,i)
            s12(:) = s12(:) - anint(s12(:))
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
          r12(:) = r12(:) - Lbox(:) * anint(r12(:)/Lbox(:))
          r23(:) = r23(:) - Lbox(:) * anint(r23(:)/Lbox(:))
          r13(:) = r13(:) - Lbox(:) * anint(r13(:)/Lbox(:))
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

  !subroutine sort_real_label(num, idx)
  !  real(8), intent(inout) :: num(:)
  !  integer, intent(inout) :: idx(:)
  !  integer :: Nele, i, j, Itemp, Nidx
  !  real(8) :: temp
  !  logical :: Lidx
  !  Nele = size(num)
  !  Nidx = size(idx)
  !  if ( Nele /= Nidx ) then
  !    print *, 'Nele and Nidx is different'
  !    stop "ERROR!!"
  !  end if

  !  do i = 1, Nele
  !    do j = i+1, Nele
  !      if (num(i) > num(j) ) then
  !        temp = num(i)
  !        num(i) = num(j)
  !        num(j) = temp

  !        if ( Lidx ) then
  !          Itemp  = idx(i)
  !          idx(i) = idx(j)
  !          idx(j) = Itemp
  !        end if
  !      end if
  !    end do
  !  end do
  !end subroutine sort_real_label
end module mod_periodic

