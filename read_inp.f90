subroutine read_input
use,intrinsic :: iso_fortran_env
use input_parameter
implicit none
integer :: i, j, Uin, ios, Ifile
character(:), allocatable :: input_file
character(len=128) :: line
character(len=128) :: FNtemp1 = "bin1.bin", FNtemp2 = "bin2.bin"

Ifile = 0

print '(" ***** START reading parameters *****")'
block
  integer :: leng
  if ( command_argument_count() == 0) then
    print '(a)',   "   There is no argument"
    print '(a,/)', '   Reading from "input.inp"'
    allocate(character(9) :: input_file)
    write(input_file,'(a)') "input.inp"
  else
    call get_command_argument(1, length=leng)
      allocate(character(leng) :: input_file)
      call get_command_argument(1, input_file)
    print '(a,a,/)', "   Reading from ", '"'//input_file//'"'
  end if
end block

open(newunit=Uin,file=input_file,status='Old',iostat=ios)
  if ( ios /= 0) then
    print '(a,a)', 'ERROR!!: There is no input file of ', input_file; stop
  end if

  InputFile:do
    read(Uin,'(a)',iostat=ios) line
!    print *, ios, line
    if ( ios == IOSTAT_END ) exit

! --- Start job type ---
    if (index(line,"# job type") == 1 )  then
      do
        read(Uin,'(a)',iostat=ios) line
        if ( ios == IOSTAT_END ) then
          print *, 'ERROR!!: There is no "# end "'
          stop
        elseif (index(line,"$Lfirst" )        == 1) then; read(Uin,*,err=901) Lfirst
        elseif (index(line, "$Nfile")         == 1) then; read(Uin,*) Nfile
        elseif (index(line, "$Job type")      == 1) then; read(Uin,*) jobtype
        elseif (index(line,"$Natom" )    == 1) then; read(Uin,*) Natom
        elseif (index(line,"$Nbeads")    == 1) then; read(Uin,*) Nbeads
        elseif (index(line,"$atom1" )    == 1) then; read(Uin,*) atom1
        elseif (index(line,"$atom2" )    == 1) then; read(Uin,*) atom2
        elseif (index(line,"$atom3" )    == 1) then; read(Uin,*) atom3
        elseif (index(line,"$atom4" )    == 1) then; read(Uin,*) atom4
        elseif (index(line,"$atom5" )    == 1) then; read(Uin,*) atom5
        elseif (index(line, "$graph_step")    == 1) then; read(Uin,*) graph_step
        elseif (index(line, "$save_beads")    == 1) then; read(Uin,*) save_beads
        elseif (index(line, "$name_binary1")  == 1) then; read(Uin,*) FNtemp1
        elseif (index(line, "$name_binary2")  == 1) then; read(Uin,*) FNtemp2
        elseif (index(line, "# end job type") == 1) then; exit
        end if
      end do
      FNameBinary1 = trim(FNtemp1)
      FNameBinary2 = trim(FNtemp2)
      allocate(DirResult(Nfile))
      allocate(FileName(Nfile))
      allocate(Ncut(Nfile))
      allocate(Nstep(Nfile))
      FileName="0"
! --- End job type ---

! --- Start input file ---
    else if (index(line,"# input file") == 1) then
      Ifile = Ifile + 1
      do
        read(Uin,'(a)',iostat=ios) line
        if ( ios == IOSTAT_END ) then
          print *, 'ERROR!!: There is no "# end "'
          stop
        elseif (index(line,"$DirResult") == 1) then; read(Uin,'(a)') DirResult(Ifile)
        elseif (index(line,"$FileName")  == 1) then; read(Uin,'(a)') FileName(Ifile)
        elseif (index(line,"$Ncut")      == 1) then; read(Uin,*) Ncut(Ifile)
        elseif (index(line,"$Nstep" )    == 1) then; read(Uin,*) Nstep(Ifile)
! I will move these parameters to 'jobtype'
          elseif (index(line,"$Natom" )    == 1) then; read(Uin,*) Natom
          elseif (index(line,"$Nbeads")    == 1) then; read(Uin,*) Nbeads
          elseif (index(line,"$atom1" )    == 1) then; read(Uin,*) atom1
          elseif (index(line,"$atom2" )    == 1) then; read(Uin,*) atom2
          elseif (index(line,"$atom3" )    == 1) then; read(Uin,*) atom3
          elseif (index(line,"$atom4" )    == 1) then; read(Uin,*) atom4
          elseif (index(line,"$atom5" )    == 1) then; read(Uin,*) atom5
! I will move these parameters to 'jobtype'
          elseif (index(line,"# end file") == 1) then; exit
        end if
      end do
! --- End input file ---

! --- Start histogram ---
    else if (index(line,"# histogram parameters") == 1 ) then
      do
        read(Uin,'(a)',iostat=ios) line
        if ( ios == IOSTAT_END ) then
          print *, 'ERROR!!: There is no "# end "'
          stop
        elseif (index(line,"$Nhist")       == 1) then; read(Uin,*) Nhist
        elseif (index(line,"$Xrange_min" ) == 1) then; read(Uin,*) hist_min1
        elseif (index(line,"$Xrange_max" ) == 1) then; read(Uin,*) hist_max1
        elseif (index(line,"$Yrange_min" ) == 1) then; read(Uin,*) hist_min2
        elseif (index(line,"$Yrange_max" ) == 1) then; read(Uin,*) hist_max2
        elseif (index(line,"$hist_margin") == 1) then; read(Uin,*) hist_margin
        elseif (index(line,"$folding")     == 1) then; read(Uin,*) Lfolding
!        elseif (index(line,"$Output_name") == 1) then; read(Uin,*) out_hist
        elseif (index(line,"# end"       ) == 1) then; exit
        end if
      end do
! --- End histogram ---

! --- Start binary calculation ---
    else if (index(line,"# binary calculation") == 1 ) then
      do
        read(Uin,'(a)',iostat=ios) line
        if ( ios == IOSTAT_END ) then
          print *, 'ERROR!!: There is no "# end "'
          stop
        elseif (index(line,"$bin_min") == 1 ) then; read(Uin,*) bin_min
        elseif (index(line,"$bin_max") == 1 ) then; read(Uin,*) bin_max
        elseif (index(line,"# end"   ) == 1 ) then; exit
        end if
      end do
! --- End binary calculation ---

! --- Start multi bond ---
    else if (index(line,"# multi bond") == 1 ) then
      do
        read(Uin,'(a)',iostat=ios) line
        if ( ios == IOSTAT_END ) then
          print *, 'ERROR!!: There is no "# end "'
          stop
        elseif (index(line,"$Nbond") == 1 ) then; read(Uin,*) Nbond
          allocate(Imulti(2,Nbond))
          read(Uin,'()')
          do i = 1, Nbond
            read(Uin,*) Imulti(:,i)
          end do
        else if (index(line,"# end"   ) == 1) then; exit
        end if
      end do
! --- End multi bond ---

! --- Start dummy atom ---
    else if (index(line,"# dummy atom") == 1 ) then
      do
        read(Uin,'(a)',iostat=ios) line
        if ( ios == IOSTAT_END ) then
          print *, 'ERROR!!: There is no "# end "'
          stop
        elseif (index(line,"$type of dummy") == 1 ) then; read(Uin,*) type_dummy
        elseif (index(line,"$atom_temp1"   ) == 1 ) then; read(Uin,*) atom_dummy1
        elseif (index(line,"$atom_temp2"   ) == 1 ) then; read(Uin,*) atom_dummy2
        elseif (index(line,"# end"         ) == 1 ) then; exit
        end if
      end do
! --- End dummy atom ---

! --- Start rotation ---
    else if (index(line,"# Rotation") == 1 ) then
      do
        read(Uin,'(a)',iostat=ios) line
        if ( ios == IOSTAT_END ) then
          print *, 'ERROR!!: There is no "# end "'
          stop
        elseif (index(line,"$Nhyd")  == 1 ) then
          read(Uin,*) Nhyd
! HERE is bug in setting the "label"
          allocate(hyd(Nhyd), r_ref(3,Natom), weight(Natom), label(Natom))
        elseif (index(line,"$Ndiv")   == 1 ) then; read(Uin,*) Ndiv
        elseif (index(line,"$Hatom")  == 1)  then
          do i = 1, Nhyd
            read(Uin,*) hyd(i)
          end do
        elseif (index(line,"$atom for cube")  == 1)  then
          read(Uin,*) atom_cube
        elseif (index(line,"$coord") == 1)  then
          do i = 1, Natom
            read(Uin,*) label(i), weight(i), r_ref(:,i)
          end do
        elseif (index(line,"$end coord") == 1)  then
          exit
        end if
      end do
! --- End rotation ---

! --- Start periodic ---
    else if (index(line,"# periodic") == 1 ) then
      do
        read(Uin,'(a)',iostat=ios) line
        if ( ios == IOSTAT_END ) then
          print *, 'ERROR!!: There is no "# end "'
          stop
        elseif (index(line,"$Lbox") == 1 ) then
          do j = 1, 3
            read(Uin,*) Lbox(j)
          end do
        elseif (index(line,"$Lattice") == 1 ) then
          do j = 1, 3
            read(Uin,*) lattice(j,:)
          end do
        elseif (index(line,"$Ielement1")        == 1 ) then; read(Uin,*) Ielement1
        elseif (index(line,"$Felement1")        == 1 ) then; read(Uin,*) Felement1
        elseif (index(line,"$Ielement2")        == 1 ) then; read(Uin,*) Ielement2
        elseif (index(line,"$Felement2")        == 1 ) then; read(Uin,*) Felement2
        elseif (index(line,"$OHO distribution") == 1 ) then
          read(Uin,'(a)') line
            read(Uin,*) Noho
            allocate(label_oho(3,Noho))
          read(Uin,'(a)') line
            do j = 1, Noho
              read(Uin,*) label_oho(:,j)
            end do
        elseif (index(trim(line) ,"# end periodic") == 1)  then
          exit
        end if
      end do
! --- End periodic ---

! --- Start umbrella ---
    else if (index(line,"# umbrella sampling") == 1 ) then
      do
        read(Uin,'(a)',iostat=ios) line
        if ( ios == IOSTAT_END ) then
          print *, 'ERROR!!: There is no "# end "'
          stop
        elseif (index(line,"$type")        == 1 )  then; read(Uin,*) umbrella_type
        elseif (index(line,"$temperature") == 1 )  then; read(Uin,*) temperature
        elseif (index(line,"$atom1")       == 1 )  then; read(Uin,*) umbrella_atom1
        elseif (index(line,"$atom2")       == 1 )  then; read(Uin,*) umbrella_atom2
        elseif (index(line,"$atom3")       == 1 )  then; read(Uin,*) umbrella_atom3
        elseif (index(line,"$force")       == 1 )  then; read(Uin,*) umbrella_force
        elseif (index(line,"# end")        == 1 )  then; exit
        end if
      end do
! --- End umbrella ---


! --- Start PbHPO4 ---
    else if (index(line,"# PbHPO4") == 1 ) then
      do
        read(Uin,'(a)',iostat=ios) line
        if ( ios == IOSTAT_END ) then
          print *, 'ERROR!!: There is no "# end "'
          stop
        elseif (index(line,"$Lattice") == 1 ) then
          do j = 1, 3
            read(Uin,*) lattice(j,:)
          end do
        elseif (index(line,"$Nunit") == 1 ) then; read(Uin,*) Nunit
        elseif (index(line,"# end")  == 1 ) then; exit
        end if
      end do

    end if
  end do InputFile
close(Uin)

! --- Print input parameters ---
  print '(" *** Input parameters as follows ***")'
  print '("   jobtype = ",I0)', jobtype
  print '("   Natom   = ",I0)', Natom
  print '("   Nbeads  = ",I0)', Nbeads
  print '("   LFirst  = ",L)',  Lfirst
  print '("   atom",I0,"  = ", I0)', 1, atom1
  print '("   atom",I0,"  = ", I0)', 2, atom2
  print '("   atom",I0,"  = ", I0)', 3, atom3
  print '("   atom",I0,"  = ", I0)', 4, atom4
  print '("   atom",I0,"  = ", I0)', 5, atom5

  do Ifile = 1, Nfile
    print '(a,I0,a)',    " *** Input from the file # ",Ifile," ***"
    print '(a,a)',  "   DirResult = ", trim(DirResult(Ifile))
    print '(a,a)',  "   FileName = ", trim(FileName(Ifile))
    print '(a,I0)', "   Ncut     = ", Ncut(Ifile)
    print '(a,I0)', "   Nstep    = ", Nstep(Ifile)
  end do
  print *, ""

  TNstep = 0
  do j = 1, Nfile
    TNstep = TNstep + Nstep(j) - Ncut(j)
  end do
  print '(a, i0)', "   The total number of step = ", TNstep
! --- End Print input parameters ---

close(Uin)
print '(a,/)', " ***** END reading parameters *****"

return
!  100 print *, 'ERROR!!: There is no "# End ~~"'; stop
!  101 print *, 'ERROR!!: There is no "# job type", check -name_binary'; stop
!  102 print *, 'ERROR!!: There is no "# input file"'; stop
!  103 print *, 'ERROR!!: There is no "# histogram parameters"'; stop
!  120 print *, 'ERROR!!: There is no "# end histogram parameters"'; stop
!  900 print *, 'ERROR!!: There is no "input.dat"'; stop
  901 print *, 'ERROR!!: "$Lbinary" must be T or F'; stop
end subroutine read_input

!! --- Erro Check !! ---
!  if     ( Natom < 0) then; print *, "ERROR!!: Write Natom!!";  stop
!  elseif (Nbeads < 0) then; print *, "ERROR!!: Write Nbeads!!"; stop
!  endif
!  do j = 1, Nfile
!    if     (Nstart(j) < 0) then; print *, "ERROR!!: Write Nstart!!"; stop
!    elseif ( Nstep(j) < 0) then; print *, "ERROR!!: Write Nstep!!";  stop
!    end if
!    do k = 1, 5
!      if (atom_num(k,j) < 0) then
!        print *, "ERROR!!: Write atom of ",k; stop
!      end if
!    end do
!  end do
!! --- End Erro Check !! ---
!
