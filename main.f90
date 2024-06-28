! Analysing PIMD data
! last modified 2024/06/25
! Program structure
! All of the data are stored in the 'data_beads', and this is analyzed according to the 'jobtype' option

! Choose "job type" as follows options:
!   1 : Bond length                 (atom1-atom2)
!   2 : Angle histgram              (atom1-atom2-atom3)
!   3 : Dihedral angle              (atom1-atom2-atom3-atom4)
!   4 : Bond diff                   (atom1-atom2  -  atom3-atom4)
!   5 : Bond add                    (atom1-atom2  +  atom3-atom4)
!   9 : 1D histgram from External
!  11 : Multi bond for all
!  12 : Multi bond sort
!  13 : Multi bond sum
!  21 : 2D histogram_bond            (atom1-atom2 and atom3-atom4)
!  22 : 2D histogram_angle           (atom1-atom2 and atom3-atom4-atom5)
!  29 : 2D histogram from External   (Arbitrary number line)
!  31 : binary mask average    (coor.xyz if (bin1) )
!  32 : binary mask each beads (coor.xyz if (bin1) )
!  33 : binary append
!  34 : binary add             (bin1 + bin2)
!  35 : binary diff            (bin1 - bin2)
!  41 : Dummy atom (X) for bond      (atom1-atomX)
!  42 : Dummy atom (X) for angle     (atom1-atom2-atomX)
!  43 : Dummy atom (X) for dihedral  (atom1-atom2-atomX-atom4)
!  51 : Beads expansion      (all atoms)
!  52 : Beads expansion      (export binary of atom1)
!  53 : Beads expansion      (atom1 projected to atom2-atom3)
!  61 : charge analysis      (all atoms)
!  62 : charge analysis      (atom1)
!  63 : dipole analysis
!  64 : hfcc analysis        (atom1)
!  65 : Force analysis
!  71 : Rotation          (movie)
!  72 : Rotation          (cube file)
!  8* : === Periodic boundary condition ===
!  81 : Radial distribution (element1)
!  82 : Radial distribution (element1 to element2)
!  83 : Bond length with periodic         (atom1-atom2)
!  84 : Bond diff with periodic           (atom1-atom2  -  atom3-atom4)
!  85 : RMSD (Root mean square deviation)
!  86 : Minimum bond length               (from atom1)
!  89 : OHO distribution
!  91 : Out of plane         (atom2-atom1-atom3 -> atom1-atom4)
! 191 : PbHPO4  (O-O distribution)
! 192 : PbHPO4  (dleta OH distribution)
!!!  29 : 2D histogram from External   (Old, use the 28 mode)
!!!  91 : projection           (atom1-atom2  T  atom3-atom4)

program analysis
use input_parameter, &
    only: Natom, Nbeads, TNstep, label, jobtype, r, data_step, data_beads, Lfirst
!use calc_centoroid
!use calc_histogram1D
use calc_histogram2D
use mod_other_quantities
use mod_special_case
use mod_periodic
implicit none


! +++ Reading the input file +++
call read_input

! r(:,i,j,k) = r(xyz,atom,beads,step)
if (Lfirst .eqv. .True.) then
  allocate(r(3,Natom,Nbeads,TNstep))
  if ( .not. allocated(label) ) allocate(label(Natom))
end if

! +++ Reading coordinate +++
call read_coor
allocate(data_step(TNstep), source=0.0d0)
allocate(data_beads(Nbeads,TNstep), source=0.0d0)

select case(jobtype)
  case(1)
    call calc_bond
  case(2)
    call calc_angle
  case(3)
    call calc_dihedral
  case(4:5)
    call calc_bond
  case(11:14)
    call multi_bond
!  case(21:22)
!    call calc_2Dhist
  case(29)
    call external_2Dhits_arbitrary
  case(31:33)
    call binary_calc
  case(41:49)
    call dummy_atom
!  case(51:54)
!    call beads_expansion
  case(61:65)
    call other_quantities
  case(71:73)
    call rotation
  case(81:89)
    call periodic
  case(91)
    call special_case
  case(191:195)
    call pbhpo4
  case default
    stop 'ERROR!!! wrong "Job type" option'
end select

!deallocate(r)
!deallocate(data_step)
!deallocate(data_beads)
print '(a,/)', " Normal termination"
stop
contains

end program analysis
