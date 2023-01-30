subroutine binary_calc
  use input_parameter,  &
      only: jobtype, Natom, Nbeads, TNstep, label, &
            atom1, atom2, bin_min, bin_max, &
            r, data_beads, data_step, save_beads, FNameBinary1, FNameBinary2
!  use calc_histogram1D, only: calc_1Dhist
!  use utility,          only: calc_deviation, calc_cumulative, reblock_step
  implicit none
  integer :: Uin, Uout, ios, i, j, k
  real(8), allocatable :: binary1(:,:), binary2(:,:)
!  logical :: L2file = .True.

  allocate(binary1(Nbeads,TNstep),binary2(Nbeads,TNstep))

  open(newunit=Uin, file=FNameBinary1, form='unformatted', access='stream', status='old', iostat=ios)
    if ( ios /= 0 ) then
      print '(a,a)', 'ERROR!!: There is no input file of ', FNameBinary1; stop
    end if

    do k = 1, TNstep
      do j = 1, Nbeads
        read(Uin) binary1(j,k)
      end do
    end do
  close(Uin)

  select case(jobtype)
    case(31)
      call binary_mask
  end select

contains
  subroutine binary_mask
    real(8) :: ave
    real(8), allocatable :: rnew(:,:,:,:)
    character(len=:), allocatable :: ftemp
    integer :: Imask, Nmask

    allocate(rnew(3,Natom,Nbeads,TNstep))
    ftemp = 'coor_mask.bin'

    Nmask = 0
    do k = 1, TNstep
      ave = sum(binary1(:,k))/dble(Nbeads)
      if ( bin_min < ave .and. ave < bin_max ) then
        Nmask = Nmask + 1
        rnew(:,:,:,Nmask) = r(:,:,:,k)
      end if
    end do

    print '(" ***** START compression *****")'
    open(newunit=Uout,file=ftemp, form='unformatted', access='stream', status='replace')
      write(Uout) Natom, Nbeads, 1, Nmask, Nmask
      write(Uout) label(:)
      do Imask = 1, Nmask
        do j = 1, Nbeads
          do i = 1, Natom
            write(Uout) rnew(:,i,j,Imask)
          end do
        end do
      end do
    close(Uout)

    print *, 'Nmask is ', Nmask

    print '(a,a,a)', '    Binary data is saved as "',ftemp,'"'
!    print '(a)', '    in "'//trim(DirResult(1))//'" directory '
    print '(" ***** END compression *****")'
    print *, ""
  end subroutine binary_mask

end subroutine binary_calc

