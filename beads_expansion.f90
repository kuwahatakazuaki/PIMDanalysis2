
subroutine beads_expansion
use input_parameter
use calc_parameter, only: r, data_beads !, data_step
use calc_histogram1D
use utility
implicit none
  integer :: i, j, k

  select case(jobtype)
    case(51)
      call beads_all
    case(52)
      call beads_bin
    case(53)
      call beads_projection
  end select


contains

  subroutine beads_bin
    real(8) :: beads_ave1(3,TNstep), beads_dev1(TNstep)
    real(8) :: beads_temp(Nbeads,TNstep)
    integer :: atom1, Ounit, step

    atom1 = atom_num(1,1)  ! atom_num(Iatom,Ifile)

    beads_ave1(:,:) = 0.0d0
    beads_dev1(:) = 0.0d0
    do j = 1, Nbeads
      beads_ave1(:,:) = beads_ave1(:,:) + r(:,atom1,j,:)
    end do
    beads_ave1(:,:) = beads_ave1(:,:) / Nbeads

    do k = 1, TNstep
      do j = 1, Nbeads
        beads_dev1(k) = beads_dev1(k) + dot_product( r(:,atom1,j,k)-beads_ave1(:,k) , r(:,atom1,j,k)-beads_ave1(:,k) )
      end do
    end do
    beads_dev1(:) = dsqrt( beads_dev1(:) / Nbeads )
    do step = 1, TNstep
      beads_temp(:,step) = beads_dev1(step)
    end do

    open(newunit=Ounit,file=FNameBinary1, form='unformatted', access='stream', status='replace')
      do step = 1, TNstep
        do j = 1, Nbeads
 !         write(Ounit) beads_dev1(step)
          write(Ounit) beads_temp(j,step)
        end do
      end do
    close(Ounit)
  end subroutine beads_bin


  subroutine beads_all
    integer :: Nbeads_back
    character(len=19) :: out_expan
!    character(len=2) :: dummy = "C"
    real(8) :: beads_ave(3,Natom,TNstep), beads_dev(Natom,TNstep)
    real(8) :: data_dev

    write(out_expan, '("beads_expansion.out")')

    beads_ave(:,:,:) = 0.0d0  ! beads_ave(3,i,k)
    beads_dev(:,:) = 0.0d0
    do j = 1, Nbeads
      beads_ave(:,:,:) = beads_ave(:,:,:) + r(:,:,j,:)
    end do
    beads_ave(:,:,:) = beads_ave(:,:,:) / Nbeads

    do i = 1, Natom
      do k = 1, TNstep
        do j = 1, Nbeads
!          beads_dev(i,k) = beads_dev(i,k) + norm(r(:,i,j,k)-beads_ave(:,i,k))**2
          beads_dev(i,k) = beads_dev(i,k) + dot_product( r(:,i,j,k)-beads_ave(:,i,k) , r(:,i,j,k)-beads_ave(:,i,k) )
        end do
      end do
    end do
    beads_dev(:,:) = dsqrt( beads_dev(:,:) / Nbeads )

    deallocate(data_beads)
    allocate(data_beads(1,TNstep))

    open(25, file=out_expan, status='replace')
      write(25,*) "     # Number of Atom = ", Natom
      write(25,*) "     # Total Step     = ", TNstep
      write(25, '(" # Ave")', advance='no')
      do i = 1, Natom
        write(25,'(F12.6)', advance='no') sum(beads_dev(i,:))/TNstep
      end do
      write(25, '(/," #    ")', advance='no')
    !  write(25,'(" #    ", 1000A12)') atom(:)
      do i = 1, Natom
        write(25,'(a10,I2)',advance='no') trim(atom(i)),i
      end do
      write(25,*) ""
      do k = 1, TNstep
        write(25,'(I6, 1000F12.6)') k, beads_dev(:,k)
      end do
    close(25)

    ! back up Nbeads_back
    Nbeads_back = Nbeads
    Nbeads = 1
    do i = 1, Natom
    !  hist_min(1) = 0.0d0
    !  hist_max(1) = 0.0d0
      data_beads(1,:) = beads_dev(i,:)
    !  data_max = maxval(data_beads)
    !  data_min = minval(data_beads)
    !  data_ave = sum(data_beads)/size(data_beads)
      call calc_deviation(data_dev)
      write(out_hist,'("hist_beads_",a,I0,".out")') trim(atom(i)), i
      call calc_1Dhist(0.0d0,0.0d0)
    end do
    Nbeads = Nbeads_back
  end subroutine beads_all

  subroutine beads_projection
  end subroutine beads_projection

end subroutine beads_expansion


