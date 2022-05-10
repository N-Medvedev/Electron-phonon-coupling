! 1111111111111111111111111111111111111111111111111111111111111
! This module contains subroutines to deal with objects/types,
! in the framework of the object oriented programming:

MODULE Objects
  implicit none 

!==============================================
! Error handling as an object:
type Error_handling
   LOGICAL Err		! indicates that an error occured
   integer Err_Num	! assign a number to an error
   integer File_Num		! number of the file with error log
   character(100) Err_descript	! describes more details about the error
end type
!==============================================

!==============================================
! Material parameters are all in here:
type :: Solid
    character(100) Name
    real(8) Mass    ! average mass of atom [a.m.u.]
    integer Ne      ! number of electrons in the valence band per molecule
    real(8) Dens    ! [g/cm^3] material density
    real(8) Ef      ! [eV] fermi-energy
    real(8) Vs      ! [m/s] speed of sound in the material
    real(8) At_Dens ! [1/cm^3] atomic material density
    real(8) q_debye ! [1/m] Debye-vector
    real(8) E_debye ! [eV] Debye energy
    real(8) Ta      ! [K] lattice temperature
end type solid
!==============================================

!==============================================
! Density of states of the material:
type :: Density_of_states
   real(8), dimension(:), allocatable :: E  ! [eV] energy
   real(8), dimension(:), allocatable :: DOS    ! [1/eV] DOS itself
   real(8), dimension(:), allocatable :: Int_DOS  ! integral of DOS
   real(8), dimension(:), allocatable :: k  ! [1/m] k-vector
end type Density_of_states
!==============================================

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! if contains subroutines, make them here:
contains


! How to write the log about an error:
subroutine Save_error_details(Err_name, Err_num, Err_data, FN)
   class(Error_handling) :: Err_name    ! object containing all details
   integer, intent(in) :: Err_num       ! number of error asigned
   character(100), intent(in) :: Err_data   ! description of the error
   integer, intent(in), optional :: FN   ! number of file where to write error log
   integer FN_err
   logical file_exist
   character Error_log_file
   Err_name%Err = .true.    ! error occured, mark it as "true"
   Err_name%Err_Num = Err_num   ! number of error we asign to it
   Err_name%Err_descript = Err_data ! descriptino of an error
   
   FN_err = Err_name%File_Num
   
   if (present(FN)) then
    write(FN_err, '(a,i2,1x,a,1x,i3)') 'Error #', Err_name%Err_Num, trim(adjustl(Err_name%Err_descript)), FN   ! write it all into the file
    write(*, '(a,i2,1x,a,1x,i3)') 'Error #', Err_name%Err_Num, trim(adjustl(Err_name%Err_descript)), FN   ! write it all into the file
   else
    write(FN_err, '(a,i2,1x,a)') 'Error #', Err_name%Err_Num, trim(adjustl(Err_name%Err_descript))   ! write it all into the file
    write(*, '(a,i2,1x,a)') 'Error #', Err_name%Err_Num, trim(adjustl(Err_name%Err_descript))   ! write it all into the file
   endif
end subroutine Save_error_details

end MODULE Objects