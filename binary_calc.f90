subroutine binary_calc
  use input_parameter,  &
      only: jobtype, Natom, Nbeads, TNstep, label, &
            atom1, atom2, bin_min, bin_max, &
            r, data_beads, data_step, save_beads, FNameBinary1, FNameBinary2
!  use utility,          only: calc_deviation, calc_cumulative, reblock_step
  implicit none
  integer :: Uin, Uout, ios, i, j, k
  real(8), allocatable :: binary1(:,:), binary2(:,:)

  allocate(binary1(Nbeads,TNstep),binary2(Nbeads,TNstep))


  select case(jobtype)
    case(31)
      call read_binary1
      call binary_mask_ave
    case(32)
      call read_binary1
      call binary_mask_each
    case(33)
      call read_binary2
      call binary_append
  end select

contains

  subroutine binary_append
    deallocate(data_beads)
    allocate(data_beads(Nbeads,2*TNstep))
    data_beads(:,1:TNstep) = binary1(:,1:TNstep)
  end subroutine binary_append

  subroutine read_binary1
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
  end subroutine read_binary1

  subroutine read_binary2
    call read_binary1
    open(newunit=Uin, file=FNameBinary2, form='unformatted', access='stream', status='old', iostat=ios)
      if ( ios /= 0 ) then
        print '(a,a)', 'ERROR!!: There is no input file of ', FNameBinary2; stop
      end if

      do k = 1, TNstep
        do j = 1, Nbeads
          read(Uin) binary1(j,k)
        end do
      end do
    close(Uin)
  end subroutine read_binary2

  subroutine binary_mask_each
    real(8), allocatable :: rnew(:,:,:)
    character(len=:), allocatable :: ftemp
    integer :: Imask, Nmask

    allocate(rnew(3,Natom,Nbeads*TNstep))
    ftemp = 'coor_mask.bin'

    Nmask = 0
    do k = 1, TNstep
      do j = 1, Nbeads
      if ( bin_min < binary1(j,k) .and. binary1(j,k) < bin_max ) then
        Nmask = Nmask + 1
        rnew(:,:,Nmask) = r(:,:,j,k)
      end if
      end do
    end do

    print '(" ***** START compression *****")'
    open(newunit=Uout,file=ftemp, form='unformatted', access='stream', status='replace')
      write(Uout) Natom, 1, 0, Nmask, Nmask
      write(Uout) label(:)
      do Imask = 1, Nmask
        do i = 1, Natom
          write(Uout) rnew(:,i,Imask)
        end do
      end do
    close(Uout)

    print *, 'Nmask is ', Nmask

    print '(a,a,a)', '    Binary data is saved as "',ftemp,'"'
    print '(" ***** END compression *****")'
    print *, ""
  end subroutine binary_mask_each

  subroutine binary_mask_ave
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
      write(Uout) Natom, Nbeads, 0, Nmask, Nmask
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
    print '(" ***** END compression *****")'
    print *, ""
  end subroutine binary_mask_ave

end subroutine binary_calc

