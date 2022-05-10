!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module contains subroutine to read in INPUT_PARAMETERS.txt file

MODULE Read_input_file
use Universal_Constants
use Objects
use Dealing_with_arrays
use Integration
implicit none

contains

subroutine read_input_parameters_file(path_sep, Matter, Num_th, Error_message, read_well)
    character(1), intent(in) :: path_sep
    type(Solid), intent(out) ::  Matter   ! name
    integer, intent(out) :: Num_th  ! number of threads for OpenMP
    type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
    logical, intent(inout) :: read_well ! did the file read well?
    integer i, FN, Reason
    logical file_exist    ! to check where file to be open exists
    logical file_opened   ! to check if a file is still opened
    character(100) :: Input_paramters, Err_data
    Input_paramters = 'INPUT_PARAMETERS.txt'    ! this name is fixed
    
    FN = 200
    inquire(file=trim(adjustl(Input_paramters)),exist=file_exist)    ! check if input file excists
    if (file_exist) then    ! read this file
        open(FN, file=trim(adjustl(Input_paramters)), action='read')
        print*, 'Material parameters are in the file '//trim(adjustl(Input_paramters))
        i = 0
        READ(FN,*,IOSTAT=Reason) Matter%Name    ! name
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2013
        
        READ(FN,*,IOSTAT=Reason) Matter%Mass    ! mass [a.m.u]
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2013
        
        READ(FN,*,IOSTAT=Reason) Matter%Ne    ! number of VB or CB electrons per atom
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2013
        
        READ(FN,*,IOSTAT=Reason) Matter%Dens    ! [g/cm^3] atomic density
        Matter%At_Dens = Matter%Dens*1d-3/(Matter%Mass*g_Mp)  ! [1/cm^3] atomic material density
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2013
        
        READ(FN,*,IOSTAT=Reason) Matter%Ef    ! [eV] Fermi-energy
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2013
       
        READ(FN,*,IOSTAT=Reason) Matter%Vs    ! speed of sound [m/s]
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2013
        !Matter%q_debye = (6.0d0*g_Pi*g_Pi*Matter%At_dens*1d6*Matter%Ne)**(1.0d0/3.0d0)  ! [1/m] Debye vector
        Matter%q_debye = (6.0d0*g_Pi*g_Pi*Matter%At_dens*1d6)**(1.0d0/3.0d0)  ! [1/m] Debye vector
        Matter%E_debye = g_h*Matter%Vs*Matter%q_debye/g_e   ! [eV] Debye energy
        
        READ(FN,*,IOSTAT=Reason) Matter%Ta    ! lattice temperature [K]
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2013
        
        READ(FN,*,IOSTAT=Reason) Num_th    ! number of openmp threads
        call read_file(Reason, i, read_well) ! reports if everything read well
        if (.not. read_well) goto 2013
    else
        Err_data = 'File '//trim(adjustl(Input_paramters))//' is not found.' ! no input data found
        call Save_error_details(Error_message, 1, Err_data)
        read_well = .false. ! failed to read the input-file
    endif

2013    continue
    inquire(unit=FN,opened=file_opened)    ! check if this file is opened
    if (file_opened) close(FN)             ! and if it is, close it
end subroutine


subroutine read_DOS_file(path_sep, Matter, Mat_DOS, Error_message, read_well)
    character(1), intent(in) :: path_sep
    type(Solid), intent(out) ::  Matter   ! name
    type(Density_of_states), intent(inout) :: Mat_DOS  ! materail DOS
    type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
    logical, intent(inout) :: read_well ! did the file read well?
    
    real(8), dimension(:,:), allocatable :: Temp_DOS
    real(8) dE, E, Sum_DOS, loc_DOS
    integer i, FN2, Reason, N, M, N_tot
    logical file_exist    ! to check where file to be open exists
    logical file_opened   ! to check if a file is still opened
    character(100) :: File_DOS, Err_data
    File_DOS = 'DOS'//path_sep//trim(adjustl(Matter%Name))//'.txt'
    FN2 = 200
    inquire(file=trim(adjustl(File_DOS)),exist=file_exist)    ! check if input file excists
    if (file_exist) then    ! read this file
       print*, 'Corresponding DOS is in the file '//trim(adjustl(File_DOS))
       open(FN2, file=trim(adjustl(File_DOS)), action='read')
       call Count_lines_in_file(FN2, N)
       if (.not. allocated(Temp_DOS)) allocate(Temp_DOS(2,N)) ! to read from file
       do i = 1, N  ! read all lines
           read(FN2, *, IOSTAT=Reason) Temp_DOS(1,i), Temp_DOS(2,i)
           call read_file_here(Reason, i, read_well)
           if (.not. read_well) goto 2014
           !print*, Temp_DOS(1,i), Temp_DOS(2,i)
       enddo
       ! Do we shift the energy start to the bottom of the CB? Comment the next line out, if not!
       !Temp_DOS(1,:) = Temp_DOS(1,:) - Temp_DOS(1,1)    ! shift the bottom level to 'zero'
       !CBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCB
       
       N_tot = 100	! number of grid points per eV
       dE = 1.0d0/dble(N_tot)	! [eV] energy grid
       M = FLOOR(100.0d0/dE)	! 100 eV is the maximal value we reach
       !print*, 'M', M, Temp_DOS(1,N), Temp_DOS(1,1)
       if (.not. allocated(Mat_DOS%E)) allocate(Mat_DOS%E(M))
       if (.not. allocated(Mat_DOS%DOS)) allocate(Mat_DOS%DOS(M))
       if (.not. allocated(Mat_DOS%Int_DOS)) allocate(Mat_DOS%Int_DOS(M))
       if (.not. allocated(Mat_DOS%k)) allocate(Mat_DOS%k(M))
       E = Temp_DOS(1,1)-dE    ! start from the start [eV]
       do i = 1, M  ! fill all points in the array:
           Mat_DOS%E(i) = E ! [eV] energy
           if (i .EQ. 1) then
                Mat_DOS%DOS(i) = 0.0d0  !Temp_DOS(2,1) ! [1/eV] DOS
           !else if (i .LT. size(Temp_DOS,2)) then
           else if (E .LT. Temp_DOS(1,size(Temp_DOS,2))) then
                call Linear_approx(Temp_DOS, E, loc_DOS, (Temp_DOS(1,1)-dE), 0.0d0)
                Mat_DOS%DOS(i) = loc_DOS ! [1/eV] DOS
           else
                Mat_DOS%DOS(i) = Temp_DOS(2,size(Temp_DOS,2))/sqrt(Temp_DOS(1,size(Temp_DOS,2)))*SQRT(E) ! [1/eV] Free-electron DOS
           endif
           E = E + dE   ! energy grid [eV]
           !write(*,'(i,f,f,f,f)') i, Mat_DOS%E(i), Mat_DOS%DOS(i), Temp_DOS(2,size(Temp_DOS,2)), Temp_DOS(2,size(Temp_DOS,2))/sqrt(Temp_DOS(1,size(Temp_DOS,2)))*SQRT(E)
       enddo
       ! Do we shift the energy start to the bottom of the CB? Comment the next line out, if not!
       !Mat_DOS%E = Mat_DOS%E - Mat_DOS%E(1) ! shift it to 'zero', and make it positive [eV]
       !CBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCBCB
       
       !Mat_DOS%E = Mat_DOS%E(size(Mat_DOS%E):1:-1) ! make the array increasing
       !Mat_DOS%DOS = Mat_DOS%DOS(size(Mat_DOS%DOS):1:-1) ! make the array according to the increasing energy
       
       call Integrate_function(0, Mat_DOS%E, Mat_DOS%DOS, Mat_DOS%E(1), Matter%Ef, Sum_DOS, Error_message)
       Mat_DOS%DOS = Mat_DOS%DOS/Sum_DOS*Matter%Ne  ! normalize to the number of electrons in the CB
       
!        print*, Sum_DOS, Matter%Ne
!        pause 'TEST DOS'
       
       call Integrate_function(0, Mat_DOS%E, Mat_DOS%DOS, Mat_DOS%E(1), 100.0d0, Mat_DOS%Int_DOS, Error_message)
       Mat_DOS%k = (3.0d0*2.0d0*g_Pi*g_Pi/2.0d0*Mat_DOS%Int_DOS*Matter%At_Dens*1d6)**(1.0d0/3.0d0)  ! [1/m]
!        do i = 1, size(Mat_DOS%E)
!         write(*,'(a,f,f,f,es)') 'DOS:', Mat_DOS%E(i), Mat_DOS%DOS(i), Mat_DOS%Int_DOS(i), Mat_DOS%k(i)
!        enddo
    else
        Err_data = 'File '//trim(adjustl(File_DOS))//' is not found.' ! no input data found
        call Save_error_details(Error_message, 2, Err_data)
        read_well = .false. ! failed to read the input-file
    endif

2014    continue
    inquire(unit=FN2,opened=file_opened)    ! check if this file is opened
    if (file_opened) close(FN2)             ! and if it is, close it
end subroutine


subroutine Count_lines_in_file(File_num, N)
    integer, INTENT(in) :: File_num     ! number of file to be opened
    integer, INTENT(out) :: N           ! number of lines in this file
    integer i
    i = 0
    do
        read(File_num, *, end=603)
        i = i + 1
    enddo
603 continue
    rewind (File_num) ! to read next time from the beginning, not continue from the line we ended now.
    N = i
end subroutine Count_lines_in_file


subroutine read_file(Reason, i, read_well)
   integer, intent(in) :: Reason    ! file number where to read from
   integer, intent(inout) :: i      ! line number
   logical, intent(inout) :: read_well  ! did we read ok?
   i = i + 1    ! it's next line
   IF (Reason .GT. 0)  THEN ! ... something wrong ...
       write(*,'(a,i3,a)') 'Problem reading input file in line ', i, ', wrong type of variable'
       read_well = .false.
   ELSE IF (Reason .LT. 0) THEN ! ... end of file reached ...
       write(*,'(a,i3,a)') 'Problem reading input file in line ', i, ', unexpected END of file'
       read_well = .false.
   ELSE   ! normal reading
       read_well = .true.  ! it read well, nothing to report
   END IF
end subroutine read_file

subroutine read_file_here(Reason, i, read_well)
   integer, intent(in) :: Reason    ! file number where to read from
   integer, intent(in) :: i      ! line number
   logical, intent(inout) :: read_well  ! did we read ok?
   IF (Reason .GT. 0)  THEN ! ... something wrong ...
       write(*,'(a,i3,a)') 'Problem reading input file in line ', i, ', wrong type of variable'
       read_well = .false.
   ELSE IF (Reason .LT. 0) THEN ! ... end of file reached ...
       write(*,'(a,i3,a)') 'Problem reading input file in line ', i, ', unexpected END of file'
       read_well = .false.
   ELSE   ! normal reading
       read_well = .true.  ! it read well, nothing to report
   END IF
end subroutine read_file_here

END MODULE Read_input_file
