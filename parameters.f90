module input_parameter
  implicit none

  integer :: Natom, Nbeads, TNstep, Nhist, Nbond, Nfile, Natom_peri
  integer :: atom1, atom2, atom3, atom4, atom5
  integer :: jobtype
  integer :: type_dummy, atom_dummy1, atom_dummy2
  integer :: graph_step = 10
  integer, allocatable :: Nstep(:), Ncut(:)
  real(8) :: hist_min1=0.0d0, hist_max1=0.0d0, hist_min2=0.0d0, hist_max2=0.0d0, hist_margin !, hist_max
  character(len=128), allocatable :: FileName(:), DirResult(:)
  integer, allocatable :: Imulti(:,:)
  logical :: Lfirst  = .False.
  logical :: Lfolding = .False.
  logical :: save_beads = .False.
  character(len=:), allocatable :: FNameBinary1, FNameBinary2
  integer :: umbrella_type, umbrella_atom1, umbrella_atom2, umbrella_atom3
  real(8) :: umbrella_force, temperature, Lbox(3)
  integer :: Ielement1, Felement1, Ielement2, Felement2, Nunit, Nhyd, Ndiv = 30
  character(len=2), allocatable :: label(:)
  integer, allocatable :: hyd(:), label_oho(:,:), Itetra(:,:)
  integer :: atom_cube, Noho, Ntetra
  real(8), allocatable :: r_ref(:,:), weight(:)  ! r_ref(xyz,Natom), weight(Natom)
  real(8) :: lattice(3,3)
  real(8) :: bin_min, bin_max

  real(8), save, allocatable :: r(:,:,:,:) ! r(xyz,i=atom,j=beads,k=step)
  real(8), save, allocatable :: data_step(:), data_beads(:,:) ! step(TNstep), beads(Nbeads,TNstep)

  real(8), parameter :: KtoAU     = 1.98624d-3/627.51d0  ! Boltzmann constant K to AU (kcal/K hartree/kcal)
  real(8), parameter :: AUtoK     = 627.51d0/1.98624d-3  ! Boltzmann constant K to AU
  real(8), parameter :: AngtoAU   = 1/0.529177249d0 ! 1.8897259886d0
  real(8), parameter :: AUtoAng   = 0.529177249d0
  real(8), parameter :: Angs2Bohr = 1.8897259886d0

!  character(len=:), allocatable :: out_hist
end module input_parameter



