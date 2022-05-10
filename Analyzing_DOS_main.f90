!***************************************************************
!* This code was developed in 2014 by
!* Dr. Nikita Medvedev, CFEL at DESY, Hamburg, Germany, 
!* The code calculates electronic parameters of a matreial
!* for given DOS vs electron temperature: 
!* 1) electron-phonon coupling factor
!* 2) electron heat capacity
!* 3) chemical potential
!* 4) electron effective mass in one-parabolic-band approximation
!* Theoretical background of this work relies on the paper:
!* B. Y. Mueller and B. Rethfeld
!* "Relaxation dynamics in laser-excited metals under nonequilibrium conditions"
!* Phys. Rev. B 87, 035139 (2013)
!* Should you have any questions, address them to the author:
!* nikita.medvedev@desy.de
!***************************************************************

! ifort.exe /F9999999999 /O3 /Qipo /Qvec-report1 /fpp /Qopenmp /D OMP_inside /heap-arrays Analyzing_DOS_main.f90 -o Analyzing_DOS.exe /link /stack:9999999999 

!***************************************************************

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Include all the separate files with modules to use in the main program:
include 'List_of_modules.f90'   ! in this file, all the modules to be included are listed
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


PROGRAM Analyzing_DOS
!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
! Initiate modules:
use Universal_Constants
use Objects
use Dealing_with_arrays
use Path_separator
use Read_input_file
use Integration
use Produce_output
use Fermi_Bose
use Temperature_stuff
use Electron_phonon_coupling
!MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
implicit none
!VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
type(Error_handling) Error_message	! error messages are dealed with as objects
type(Solid) :: Matter   ! all material parameters, will be read from a file
type(Density_of_states) :: Mat_DOS  ! DOS of the given material, will be read from a file
real(8), dimension(:), allocatable :: Te_grid   ! [K]
real(8), dimension(:), allocatable :: mu_grid   ! [eV]
real(8), dimension(:), allocatable :: Ce   ! electron heat capacity
real(8), dimension(:), allocatable :: ksc   ! [1/m] screening length
real(8), dimension(:), allocatable :: gfactor   ! electron-phonon coupling factor
real(8) as1, Etot, Te, mu
integer i, j, k, ctim(8), c1(8)
logical read_well
character(100) File_name
CHARACTER(len=1) :: path_sep
! For OPENMP:
integer Num_th, my_id, OMP_GET_THREAD_NUM, OMP_GET_NUM_THREADS
read_well = .true.  ! assuming all files read well
! Error log file:
Error_message%Err = .false. ! no error at the start
File_name = 'OUTPUT_Error_log.dat'
open(unit = 100, FILE = trim(adjustl(File_name)))
Error_message%File_Num = 100	! file number with error-log
   
call get_path_separator(path_sep, Error_message, read_well)   ! module: Path_separator
if (.not. read_well) goto 2012  ! if there was an error, there is nothing else to do, go to end
call date_and_time(values=c1) ! standard FORTRAN time and date
 ctim=c1
write(*, 1005) ctim(5), ctim(6), ctim(7), ctim(3), ctim(2), ctim(1)

! Read input files:
call read_input_parameters_file(path_sep, Matter, Num_th, Error_message, read_well) ! input data
if (.not. read_well) goto 2012
call read_DOS_file(path_sep, Matter, Mat_DOS, Error_message, read_well)  ! DOS for the material, specifiedin the input data
if (.not. read_well) goto 2012
!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
call OMP_SET_DYNAMIC(0) 	        ! standard openmp subroutine
call OMP_SET_NUM_THREADS(Num_th)    ! start using threads with openmp: Num_th is the number of threads, defined in the input file
!IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII



call Temperature_grid(Te_grid)  ! create temperature grid [K]
call date_and_time(values=c1)	    ! For calculation of the time of execution of the program
write(*, 1010) 'Grid is done at 	', c1(5),c1(6),c1(7),c1(8), c1(3), c1(2), c1(1)


!------
! TEST
!  call get_E_from_T_array(Matter, Mat_DOS, Te_grid, mu_grid, Error_message)
!  pause 'Calculated Ee, mu'
! !------
! !------
! ! TEST
! call get_T_from_E_array(Matter, Mat_DOS, Te_grid, mu_grid, Error_message)
! pause 'Calculated T, mu'
!------


call Chemical_potential(Te_grid, Matter, Mat_DOS, Error_message, mu_grid)   ! calculate chemical potential
call date_and_time(values=c1)	    ! For calculation of the time of execution of the program
write(*, 1010) 'Chemical potential at	', c1(5),c1(6),c1(7),c1(8), c1(3), c1(2), c1(1)

call heat_capacity(Te_grid, mu_grid, Matter, Mat_DOS, Error_message, Ce, ksc)    ! calculate electron heat capacity
call date_and_time(values=c1)	    ! For calculation of the time of execution of the program
write(*, 1010) 'Heat capacity at	', c1(5),c1(6),c1(7),c1(8), c1(3), c1(2), c1(1)

!OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
! Saving already existing output data into files (without coupling parameter yet):
call save_output(path_sep, Matter, Mat_DOS, Te_grid, mu_grid, Ce, ksc, gfactor, .true.)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



call date_and_time(values=c1)	    ! For calculation of the time of execution of the program
write(*, 1010) 'Coupling starts at	', c1(5),c1(6),c1(7),c1(8), c1(3), c1(2), c1(1)
call Electron_phonon(Te_grid, mu_grid, Matter, Mat_DOS, ksc, Error_message, gfactor)  ! coupling parameter
call date_and_time(values=c1)	    ! For calculation of the time of execution of the program
write(*, 1010) 'Coupling is done at	', c1(5),c1(6),c1(7),c1(8), c1(3), c1(2), c1(1)

!OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
! Saving full output data into files:
call save_output(path_sep, Matter, Mat_DOS, Te_grid, mu_grid, Ce, ksc, gfactor)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End of the program:
! Printing out the duration of the program, starting and ending time and date:
2012 call date_and_time(values=c1)	    ! For calculation of the time of execution of the program
as1=REAL(24*60*60*(c1(3)-ctim(3))+3600*(c1(5)-ctim(5))+60*(c1(6)-ctim(6))+(c1(7)-ctim(7))+(c1(8)-ctim(8))*0.001)	! sec
print*, '   '
write(*,'(a, es16.8, a)') 'Duration of execution of the program: ', as1, ' [sec]'
print*, '   '
write(*, 1001) ctim(5),ctim(6),ctim(7), ctim(3), ctim(2), ctim(1)
write(*, 1002) c1(5), c1(6), c1(7), c1(3),c1(2), c1(1)		
print*, '   '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Closing the remaining opened files:
if (Error_message%Err) then ! if error occured (thus it's = "true") then save the error log file:
   close(100)
else ! if there was no error, no need to keep the file, delete it:
   close(100, status='delete')
endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Formats defined for printing out on the screen:
1001 format ('Begining: ', i2.2, ':', i2.2, ':', i2.2, '  ', i2.2, '/', i2.2, '/', i4.4)
1002 format ('The end:  ', i2.2, ':', i2.2, ':', i2.2, '  ', i2.2, '/', i2.2, '/', i4.4)
1005 format ('Start at: ', i2.2, ':', i2.2, ':', i2.2, '  ', i2.2, '/', i2.2, '/', i4.4)
1006 format ('Step at: ', i2.2, ':', i2.2, ':', i2.2, '  ', i2.2, '/', i2.2, '/', i4.4)
1007 format ('Step at: ', i2.2, ':', i2.2, ':', i2.2, ':', i3.3, '  ', i2.2, '/', i2.2, '/', i4.4)
1010 format (a, i2.2, ':', i2.2, ':', i2.2, ':', i3.3, '  ', i2.2, '/', i2.2, '/', i4.4)
1008 format (a, i4, a, i6, a, i2.2, ':', i2.2, ':', i2.2)

if (path_sep .EQ. '\') then
!     PAUSE 'The program is finished, press RETURN to go out...'
endif    
STOP
END PROGRAM Analyzing_DOS