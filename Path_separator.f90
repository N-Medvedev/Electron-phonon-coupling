!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module contains subroutine to get the environmental variable
! for the path separator, which is different for windows and linux

MODULE Path_separator
use Objects
implicit none

contains
subroutine get_path_separator(path_sep, Error_message, read_well)
    type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
    CHARACTER(len=1), intent(out) :: path_sep   ! path separator
    logical, intent(inout) :: read_well ! did the data read well?
    CHARACTER(len = 100) :: path, Err_data
    CALL get_environment_variable("PATH",path)
    if (path(1:1) .EQ. '/') then        !unix based OS
        path_sep = '/'
    else if (path(3:3) .EQ. '\') then   !Windows OS
        path_sep = '\'
    else
        Err_data = 'Path separator is not defined' ! unknown OS
        call Save_error_details(Error_message, 0, Err_data)
        read_well = .false. ! didn't read well this data
    endif
end subroutine


END MODULE Path_separator