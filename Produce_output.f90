!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module contains subroutine to print out output file(s)

MODULE Produce_output
use Universal_Constants
use Objects
implicit none

contains

subroutine save_output(path_sep, Matter, Mat_DOS, Te, mu, Ce, ksc, gfactor, part)
    character(1), intent(in), optional :: path_sep
    type(Solid), intent(in) ::  Matter   ! name
    type(Density_of_states), intent(in) :: Mat_DOS  ! materail DOS
    real(8), dimension(:), intent(in) :: Te, mu, Ce, ksc, gfactor
    logical, intent(in), optional :: part ! save only part of the output?
    
    integer FN, i, FN2
    character(100) Output_file, File_name, ch1
    logical file_exist, file_opened
    
!     print*, 'save_output 0'
    FN = 300
    write(Output_file,'(a,a,a)') 'OUTPUT_', trim(adjustl(Matter%Name)), '_DOS'
    if (present(part)) then
       if (part) write(Output_file,'(a,a)') trim(adjustl(Output_file)), '_part'
    endif
    write(File_name,'(a,a)') trim(adjustl(Output_file)), '.dat'
    inquire(file=trim(adjustl(File_name)),exist=file_exist)    ! check if input file excists
    i = 0
    do while (file_exist)
        i = i + 1
        write(ch1,'(i6)') i
        write(File_name,'(a,a,a,a)') trim(adjustl(Output_file)), '_', trim(adjustl(ch1)), '.dat'
        inquire(file=trim(adjustl(File_name)),exist=file_exist)    ! check if input file excists
    enddo
    open(unit = FN, FILE = trim(adjustl(File_name)))
    inquire(file=trim(adjustl(File_name)),opened=file_opened)    ! check if input file excists
    
    print*, 'Output DOS-data are saved in the file: ', trim(adjustl(File_name))
    
    write(FN,'(a)') '#Energy    DOS k_vector'
    write(FN,'(a)') '#eV    eV^-1    m^-1'
    do i = 1, size(Mat_DOS%E)
        write(FN, '(f,f,es)') Mat_DOS%E(i), Mat_DOS%DOS(i), Mat_DOS%k(i)
    enddo
    
    
    FN2 = 301
    write(Output_file,'(a,a,a)') 'OUTPUT_', trim(adjustl(Matter%Name)), '_electron_properties'
    if (present(part)) then   
       if (part) write(Output_file,'(a,a)') trim(adjustl(Output_file)), '_part'
    endif
    write(File_name,'(a,a)') trim(adjustl(Output_file)), '.dat'
    inquire(file=trim(adjustl(File_name)),exist=file_exist)    ! check if input file excists
    i = 0
    do while (file_exist)
        i = i + 1
        write(ch1,'(i6)') i
        write(File_name,'(a,a,a,a)') trim(adjustl(Output_file)), '_', trim(adjustl(ch1)), '.dat'
        inquire(file=trim(adjustl(File_name)),exist=file_exist)    ! check if input file excists
    enddo
    open(unit = FN2, FILE = trim(adjustl(File_name)))
    inquire(file=trim(adjustl(File_name)),opened=file_opened)    ! check if input file excists
    
    print*, 'Output electron properties are saved in the file: ', trim(adjustl(File_name))
    
    if (present(part)) then
     if (part) then
       write(FN2,'(a)') '#Temperature    Chem.Potential     Heat_capacity   Screening   Debye_screening'
       write(FN2,'(a)') '#K    eV   J/m^3/K     1/m	1/m'
       do i = 1, size(Te) 
           write(FN2, '(es,es,es,es,es, es)') Te(i), mu(i), Ce(i), ksc(i), sqrt( Matter%At_dens*Matter%Ne*1d6*g_e*g_kb/(g_e0*Te(i)) )
       enddo
     else ! full output including coupling parameter:
       write(FN2,'(a)') '#Temperature    Chem.Potential     Heat_capacity   Screening   G_factor    Debye_screening'
       write(FN2,'(a)') '#K    eV   J/m^3/K     1/m     J/(s*K*m^3)     1/m'
       do i = 1, size(Te) 
           write(FN2, '(es,es,es,es,es, es)') Te(i), mu(i), Ce(i), ksc(i), gfactor(i), sqrt( Matter%At_dens*Matter%Ne*1d6*g_e*g_kb/(g_e0*Te(i)) )
       enddo
     endif
    else
      write(FN2,'(a)') '#Temperature    Chem.Potential     Heat_capacity   Screening   G_factor    Debye_screening'
      write(FN2,'(a)') '#K    eV   J/m^3/K     1/m     J/(s*K*m^3)     1/m'
      do i = 1, size(Te) 
         write(FN2, '(es,es,es,es,es, es)') Te(i), mu(i), Ce(i), ksc(i), gfactor(i), sqrt( Matter%At_dens*Matter%Ne*1d6*g_e*g_kb/(g_e0*Te(i)) )
      enddo 
    endif
    
2015    continue
    inquire(unit=FN,opened=file_opened)    ! check if this file is opened
    if (file_opened) close(FN)             ! and if it is, close it
    inquire(unit=FN2,opened=file_opened)    ! check if this file is opened
    if (file_opened) close(FN2)             ! and if it is, close it
end subroutine 

END MODULE Produce_output