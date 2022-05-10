!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This module contains subroutines to deal with temperature-related things

MODULE Temperature_stuff
use Universal_Constants
use Objects
use Fermi_Bose
use Dealing_with_arrays
implicit none

 contains

subroutine Temperature_grid(T)
    real(8), dimension(:), allocatable, intent(out) :: T    ! [K]
    real(8) Te, T_limit, T_start
    integer i

    T_start = 100.0d0 ! [K] start from this temperature
    T_limit = 5.0d4 ! until this electronic temperature [Ks]

    if (.not. allocated(T)) then
        i = 0
        Te = T_start    ! [K]
        do while (Te .LT. T_limit)   ! [K]
            i = i + 1
            Te = next_temperature_step(Te, 1)
            !print*, i, Te
        enddo
        allocate(T(i))
    endif
    i = 0
    Te = T_start    ! [K]
    do while (Te .LT. T_limit)   ! [K]
        i = i + 1
        Te = next_temperature_step(Te, 1)
        T(i) = Te
    enddo
end subroutine


function next_temperature_step(Te, ind)
   real(8) :: next_temperature_step ! [K]
   real(8) :: Te ! [K] previous temperature
   integer :: ind ! which scale to use
   real(8) :: Te1
   select case (ind)
   case (0) ! just linear scale
      Te1 = Te + 1.0d0
   case default
      if (Te .LT. 1.0d3) then
         Te1 = Te + 10.0d0
      else if (Te .LT. 1d4) then
         Te1 = Te + 100.0d0
      else if (Te .LT. 1d5) then
         Te1 = Te + 1.0d3
      else 
         Te1 = Te + 1.0d4
      endif
   end select
   next_temperature_step = Te1
end function next_temperature_step


subroutine Chemical_potential(Te, Matter, Mat_DOS, Error_message, mu) ! calculate chemical potential
    real(8), dimension(:), intent(in) :: Te    ! [K]
    real(8), dimension(:), allocatable, intent(out) :: mu    ! [eV]
    type(Solid), intent(in) ::  Matter
    type(Density_of_states), intent(in) :: Mat_DOS  ! materail DOS
    type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
    real(8) T, mu_temp
    real(8), dimension(:), allocatable :: mu_ch
    integer i, N, my_id, OMP_GET_THREAD_NUM
    N = size(Te)
    if (.not. allocated(mu)) allocate(mu(N))
    mu = 0.0d0
    if (.not. allocated(mu_ch)) allocate(mu_ch(N))
    mu_ch = 0.0d0
    !print*, 'Start of the parallel region'
    !$omp parallel &
    !$omp private (i, my_id, T, mu_temp)
    !$omp do reduction( + : mu_ch)
    do i = 1, N ! all gridpoints
        T = Te(i)
        call get_chem_pot(T, Matter, Mat_DOS, mu_temp, Error_message)
        !mu(i) = mu_temp
        mu_ch(i) = mu_temp
        my_id = 1 + OMP_GET_THREAD_NUM() ! identify which thread it is
         write(*,'(a,i4,i2, f, f)') 'Mu', i, my_id, T, mu_temp
    enddo
    !$omp end do
    !$omp end parallel
    mu = mu_ch
    !print*, mu
    if (allocated(mu_ch)) deallocate(mu_ch)
    !pause 'CHECKING'
end subroutine


subroutine get_chem_pot(T, Matter, Mat_DOS, mu, Error_message, ne_out)
    real(8), intent(in)  :: T    ! [K]
    real(8), intent(out) :: mu    ! [eV]
    type(Solid), intent(in) ::  Matter
    type(Density_of_states), intent(in) :: Mat_DOS  ! materail DOS
    type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
    real(8), intent(out), optional :: ne_out ! electron density
    real(8) ne, a, b
    integer i
    ! start to look for it in this interval:
    a = Matter%Ef - 10.0d0*T/g_kb
    b = Matter%Ef + 10.0d0*T/g_kb
    i = 0
    ne = 10.0d20 ! to start
    do while (abs(ne-dble(Matter%Ne))/dble(Matter%Ne) .GT. 1d-10)
        i = i + 1   ! count interations
        mu = (a+b)/2.0d0
        call Chem_pot_for_Te(T, Matter, Mat_DOS, mu, ne, Error_message)
        if (ne .GT. dble(Matter%Ne)) then
            b = mu
        else
            a = mu
        endif
        if (abs(a-b) .LT. 1d-10) exit
!         print*, 'mu=', a, b, ne
    enddo
    if (abs(ne-dble(Matter%Ne))/dble(Matter%Ne) .GT. 1d-5) then
       write(*,'(a,f,f)') 'get_chem_pot : calculations did not converge!', T, mu
       write(*,'(a,f,f,e,a)') 'Ne=', ne, dble(Matter%Ne), abs(ne-dble(Matter%Ne))/dble(Matter%Ne)*100.0d0, '%'
    endif
    if (present(ne_out)) ne_out = ne
end subroutine


subroutine Chem_pot_for_Te(T, Matter, Mat_DOS, mu, ne, Error_message)
    real(8), intent(in)  :: T    ! [K]
    real(8), intent(in) :: mu    ! [eV]
    type(Solid), intent(in) ::  Matter
    type(Density_of_states), intent(in) :: Mat_DOS  ! materail DOS
    real(8), intent(out) :: ne  ! calculated ne for given Te and mu
    type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
    real(8) Ferm, f, f0, T_cur, mu_cur, E, E0, dE, third, third2
    real(8) f1, f2, f3, f4, g1, g2, g3, g4, fmin, fmax, dEmin, dEmax
    integer i, N
    third = 1.0d0/3.0d0
    third2 = 2.0d0/3.0d0
    T_cur = T/g_kb
    mu_cur = mu
    f0 = 1.0d0
    E = Mat_DOS%E(1)
    f1 = Fermi(T_cur, mu_cur, E)
    ne = 0.0d0
    do while(f1 .GT. 1d-12)
        dE = 0.1d0
        call adjust_dE(T_cur, mu_cur, 1.0d0, E, Mat_DOS, .false., dE, .false.) ! below
        
        f1 = Fermi(T_cur, mu_cur, E)
        f2 = Fermi(T_cur, mu_cur, E+third*dE) 
        f3 = Fermi(T_cur, mu_cur, E+third2*dE)
        f4 = Fermi(T_cur, mu_cur, E+dE)
        f0 = f4
        call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E, g1)
        call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E+third*dE, g2)
        call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E+third2*dE, g3)
        call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E+dE, g4)
        ne = ne + dE/8.0d0*(f1*g1 + 3.0d0*(f2*g2+f3*g3) + f4*g4)    ! Simpsone 3/8 rule
        E = E + dE
!         write(*,'(f,f,es,es,es)') T_cur, E, dE, ne
    enddo
end subroutine Chem_pot_for_Te


!EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE



subroutine get_E_from_T_array(Matter, Mat_DOS, Te_grid, mu_grid, Error_message)
    type(Solid), intent(in) ::  Matter
    type(Density_of_states), intent(in) :: Mat_DOS  ! materail DOS
    real(8), dimension(:), intent(inout), allocatable :: Te_grid, mu_grid    ! [K], [eV]
    type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
    !-----------------------------------------------------
    real(8) :: Te, Te0	! [eV] electron temperature corresponging to the given energy
    real(8) :: mu	! [eV] chemical potential
    real(8) :: Etot, Etot0, dT ! [eV]
    real(8) :: ne, Ce
    integer i, FN
    FN = 100
    open(unit = FN, FILE = 'TEST_T_E_mu.txt')
    ! Starting point: T=0

    if (.not.allocated(mu_grid)) then
       allocate(mu_grid(size(Te_grid)))
       open(unit = 111, FILE = 'OUTPUT_Ni_electron_properties_part.dat')
       read(111,*)
       read(111,*)
       do i = 1, size(Te_grid)
          read(111,'(e,e)') Te_grid(i), mu_grid(i)
!           print*, i, Te_grid(i), mu_grid(i)
       enddo
       close(111)
!        pause 'TEST 000'
    endif

    call E_for_given_Te_mu(100.0d0/g_kb, Matter, Mat_DOS, Matter%Ef, Etot, Error_message)
    call E_for_given_Te_mu(0.0d0/g_kb, Matter, Mat_DOS, Matter%Ef, Etot0, Error_message)
    
    print*, 'Ef=', Matter%Ef, Etot
    
    dT = 0.01d0 ! [eV]
    !Etot = Etot - 10.0d0 ! TEST from smaller energies
    Te0 = 0.0d0

    do i = 1, 1000
       !call get_T_from_E(Matter, Mat_DOS, Etot, Te_grid, mu_grid, Te, mu, ne, Error_message)
       call get_chem_pot(Te*g_kb, Matter, Mat_DOS, mu, Error_message, ne)	! -> mu [eV]
       call E_for_given_Te_mu(Te, Matter, Mat_DOS, mu, Etot, Error_message)
       Ce = (Etot-Etot0)/(Te-Te0)*Matter%At_Dens*1d6*g_e/g_kb
       write(FN,'(f,f,f,f,e)') Te, Etot, mu, ne, Ce
       write(*,'(f,f,f,f,e)') Te, Etot, mu, ne, Ce
       Te0 = Te
       Etot0 = Etot
       Te = Te + dT
    enddo
    
    close(FN)
end subroutine get_E_from_T_array



subroutine get_T_from_E_array(Matter, Mat_DOS, Te_grid, mu_grid, Error_message)
    type(Solid), intent(in) ::  Matter
    type(Density_of_states), intent(in) :: Mat_DOS  ! materail DOS
    real(8), dimension(:), intent(inout), allocatable :: Te_grid, mu_grid    ! [K], [eV]
    type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
    !-----------------------------------------------------
    real(8) :: Te, Te0	! [K] electron temperature corresponging to the given energy
    real(8) :: mu	! [eV] chemical potential
    real(8) :: Etot, Etot0, dE ! [eV]
    real(8) :: ne, Ce
    integer i, FN
    FN = 100
    open(unit = FN, FILE = 'TEST_E_T_mu.txt')
    ! Starting point: T=0

    if (.not.allocated(mu_grid)) then
       allocate(mu_grid(size(Te_grid)))
       open(unit = 111, FILE = 'OUTPUT_Ni_electron_properties_part.dat')
       read(111,*)
       read(111,*)
       do i = 1, size(Te_grid)
          read(111,'(e,e)') Te_grid(i), mu_grid(i)
!           print*, i, Te_grid(i), mu_grid(i)
       enddo
       close(111)
!        pause 'TEST 000'
    endif

    call E_for_given_Te_mu(100.0d0/g_kb, Matter, Mat_DOS, Matter%Ef, Etot, Error_message)
    Te0 = 0.0d0
    Etot0 = Etot
    print*, 'Ef=', Matter%Ef, Etot
    
    dE = 0.1d0 ! [eV]
    !Etot = Etot - 10.0d0 ! TEST from smaller energies
    do i = 1, 1000
       call get_T_from_E(Matter, Mat_DOS, Etot, Te_grid, mu_grid, Te, mu, ne, Error_message)
       Ce = (Etot-Etot0)/(Te-Te0)*g_e/g_kb*(Matter%At_Dens)*1d6
       write(FN,'(f,f,f,f,e)') Te, Etot, mu, ne, Ce
       write(*,'(f,f,f,f,e)') Te, Etot, mu, ne, Ce
       Te0 = Te
       Etot0 = Etot
       Etot = Etot + dE
    enddo
    
    close(FN)
end subroutine get_T_from_E_array


subroutine get_T_from_E(Matter, Mat_DOS, Etot, Te_grid, mu_grid, T, mu, ne, Error_message)
    type(Solid), intent(in) ::  Matter
    type(Density_of_states), intent(in) :: Mat_DOS  ! materail DOS
    real(8), intent(in) :: Etot ! [eV]
    real(8), dimension(:), intent(in) :: Te_grid, mu_grid    ! [K], [eV]
    real(8), intent(out) :: T	! [K] electron temperature corresponging to the given energy
    real(8), intent(out) :: mu	! [eV] chemical potential
    real(8), intent(out) :: ne	! electron density
    type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
    !-----------------------------------------------------
    real(8) :: Ef, Te, mu0, mu1, Te0, Etot0, Etot1, ne_cur, a, b
    integer :: i, j
    mu0 = Matter%Ef    ! [eV] Fermi-energy
    Te0 = 10.0d0/g_kb   ! to start with [eV]
    Te = 1.0d6/g_kb ! just to start  [eV]
    call get_chem_pot(Te0, Matter, Mat_DOS, mu1, Error_message) ! -> mu [eV]
    
    call Find_in_array_monoton(Te_grid, Te*g_kb, j) ! get number of corresponding temperature
    mu1 = mu_grid(j) ! [eV]
       
       
!     print*, 'TEST 0', Te0-Te
    Etot0 = Etot + 10.0d0
    
    i = 0
    a = Te0
    b = Te
    Te = (a+b)/2.0d0
    write(*,'(i,f,f,f,f)') i, Te*g_kb, mu1,  Etot0, Etot
    
    do while( abs(Etot0-Etot)/abs(Etot) >= 1.0d-8 )
       i = i + 1
       !call get_T(Te, Matter, Mat_DOS, mu0, Etot, Error_message)	! -> Te [eV]
       call get_chem_pot(Te*g_kb, Matter, Mat_DOS, mu1, Error_message, ne_cur)	! -> mu [eV]
       !call Find_in_array_monoton(Te_grid, Te*g_kb, j) ! get number of corresponding temperature
       !mu1 = mu_grid(j) ! [eV]
       call E_for_given_Te_mu(Te, Matter, Mat_DOS, mu1, Etot0, Error_message)
        !print*, 'Changing', Etot0, Etot
        write(*,'(i,f,f,f,f)') i, Te, mu1,  Etot0, Etot
       if (Etot0 > Etot) then
          b = Te
       else
          a = Te
       endif
       Te = (a+b)/2.0d0
       
       mu0 = mu1
       if (i > 100) then
          print*, 'get_T_from_E : DID NOT CONVERGE after iteration #', i
          exit
       endif
    enddo
    T = Te*g_kb ! -> T [K]
    mu = mu1	! -> mu [eV]
    ne = ne_cur
end subroutine get_T_from_E


subroutine get_T(T, Matter, Mat_DOS, mu, Etot, Error_message)
    real(8), intent(out)  :: T   ! [eV]
    real(8), intent(in) :: mu    ! [eV]
    real(8), intent(in) :: Etot  ! [eV]
    type(Solid), intent(in) ::  Matter
    type(Density_of_states), intent(in) :: Mat_DOS  ! materail DOS
    type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
    real(8) ne, a, b, Ee, Te
    integer i
    ! start to look for it in this interval:
    a = 1.0d0/g_kb	! starting lower limit
    b = 25.0d0		! starting upper limit
    i = 0
    Ee = 10.0d20 ! to start
    do while ( abs(abs(Ee)-abs(Etot))/abs(Etot) .GT. 1d-8 )
        i = i + 1   ! count interations
        Te = (a+b)/2.0d0
        call E_for_given_Te_mu(Te, Matter, Mat_DOS, mu, Ee, Error_message)
        if (Ee >= Etot) then
            b = Te
        else
            a = Te
        endif
        if (abs(a-b) .LT. 1d-10) exit
!         print*, 'mu=', a, b, Ee, Etot, abs(abs(Ee)-abs(Etot))/abs(Etot)
    enddo
!     pause 'CHECK'
    T = Te
!     print*, 'Chemical potential', mu
    
    if (abs(abs(Ee)-abs(Etot))/abs(Etot) .GT. 1d-5) then
       write(*,'(a,f,f)') 'get_T : calculations did not converge!', T, mu
       write(*,'(a,es,es,es,a)') 'Ee=', Ee, Etot, abs(abs(Ee)-abs(Etot))/abs(Etot)*100.0d0, '%'
    endif
end subroutine get_T



subroutine E_for_given_Te_mu(T, Matter, Mat_DOS, mu, Ee, Error_message)
    real(8), intent(in)  :: T    ! [eV]
    real(8), intent(in) :: mu    ! [eV]
    type(Solid), intent(in) ::  Matter
    type(Density_of_states), intent(in) :: Mat_DOS  ! materail DOS
    real(8), intent(out) :: Ee  ! calculated Ee for given Te and mu
    type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
    real(8) Ferm, f, f0, T_cur, mu_cur, E0, dE, third, third2
    real(8) f1, f2, f3, f4, g1, g2, g3, g4, fmin, fmax, dEmin, dEmax
    real(8) E, E1, E2, E3
    integer i, N
    third = 1.0d0/3.0d0
    third2 = 2.0d0/3.0d0
    T_cur = T
    mu_cur = mu
    E = Mat_DOS%E(1)
    f0 = Fermi(T_cur, mu_cur, E)
    Ee = 0.0d0
    do while(f0 .GT. 1d-12)
        dE = 0.1d0
        call adjust_dE(T_cur, mu_cur, 1.0d0, E, Mat_DOS, .false., dE, .false.) ! below
        E1 = E+third*dE
        E2 = E+third2*dE
        E3 = E+dE
        f0 = Fermi(T_cur, mu_cur, E)
        f1 = f0*E
        f2 = Fermi(T_cur, mu_cur, E1)*E1 
        f3 = Fermi(T_cur, mu_cur, E2)*E2
        f4 = Fermi(T_cur, mu_cur, E3)*E3
        call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E, g1)
        call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E1, g2)
        call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E2, g3)
        call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E3, g4)
        Ee = Ee + dE/8.0d0*(f1*g1 + 3.0d0*(f2*g2+f3*g3) + f4*g4)    ! Simpsone 3/8 rule
        E = E + dE
!         write(*,'(f,f,es,es,es)') T_cur, E, dE, ne
    enddo
end subroutine E_for_given_Te_mu


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
subroutine heat_capacity(Te_grid, mu_grid, Matter, Mat_DOS, Error_message, Ce, ksc)    ! calculate electron heat capacity
    real(8), dimension(:), intent(in) :: Te_grid    ! [K]
    real(8), dimension(:), intent(in) :: mu_grid    ! [eV]
    type(Solid), intent(in) ::  Matter
    type(Density_of_states), intent(in) :: Mat_DOS  ! materail DOS
    real(8), dimension(:), allocatable, intent(out) :: Ce  ! calculated Ce
    real(8), dimension(:), allocatable, intent(out) :: ksc   ! [1/m] screening length
    type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
    real(8), dimension(:), allocatable :: Ce_ch  ! calculated Ce - temporary for openmp
    real(8), dimension(:), allocatable :: ksc_ch  ! calculated ksc - temporary for openmp
    real(8) Ce_cur, ksc_cur
    integer i, N, my_id, OMP_GET_THREAD_NUM
    character(10) :: time, proper_time
    N = size(Te_grid)
    if (.not. allocated(Ce)) allocate(Ce(N))
    if (.not. allocated(ksc)) allocate(ksc(N))
    if (.not. allocated(Ce_ch)) allocate(Ce_ch(N))  ! temporary
    if (.not. allocated(ksc_ch)) allocate(ksc_ch(N))  ! temporary
    Ce_ch = 0.0d0
    ksc_ch = 0.0d0
!     call date_and_time(TIME=time) ! system function
!     print*, 'Start calculations of heat capacity at', time
    !$omp parallel &
    !$omp private (i, my_id, Ce_cur, ksc_cur, time)
    !$omp do reduction( + : Ce_ch, ksc_ch)
    do i = 1, N
        call get_Ce(i, Te_grid, mu_grid, Matter, Mat_DOS, Error_message, Ce_cur, ksc_cur)
        !Ce(i) = Ce_cur*Matter%At_Dens*1d6*g_e/g_kb
        Ce_ch(i) = Ce_cur*Matter%At_Dens*1d6*g_e/g_kb
        !ksc(i) = sqrt(abs(ksc_cur*Matter%At_Dens*1d6*g_e/g_e0))
        ksc_ch(i) = sqrt(abs(ksc_cur*Matter%At_Dens*1d6*g_e/g_e0))
        my_id = 1 + OMP_GET_THREAD_NUM()
!         call date_and_time(TIME=time) ! system function
!         write(*,'(a,i2,a,i3,a,f7.0,a)') 'Node# ', my_id, ' Point# ', i, ' Temperature ', Te_grid(i), ' '//time
    enddo
    !$omp end do
    !$omp end parallel
    ! Correct:
    ksc = ksc_ch
    ! Testing:
    !ksc = ksc_ch*sqrt(2.0d0)	! artificially increased!
!     ksc = ksc_ch/1.62d0	! artificially decreased!
    
    Ce = Ce_ch
    if (allocated(Ce_ch)) deallocate(Ce_ch)  ! temporary
    if (allocated(ksc_ch)) deallocate(ksc_ch)  ! temporary
    !print*, 'Dulong-Petit limit:', 1.5d0/g_kb*g_e*Matter%At_Dens*1d6*dble(Matter%Ne)
    !pause 'KSC'
!     call date_and_time(TIME=time) ! system function
!     print*, 'Heat capacity is done at: ', time
end subroutine



subroutine get_Ce(NT, T, mu, Matter, Mat_DOS, Error_message, Ce_cur, ksc_cur)
    integer, intent(in) :: NT   ! number of temperature point
    real(8), dimension(:), intent(in) :: T    ! [K]
    real(8), dimension(:), intent(in) :: mu   ! [eV]
    type(Solid), intent(in) ::  Matter
    type(Density_of_states), intent(in) :: Mat_DOS  ! materail DOS
    real(8), intent(out) :: Ce_cur, ksc_cur  ! calculated Ce and ksc
    type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
    real(8) Ferm, f, f0, T_cur, E, E0, dE, third, third2, dmu
    real(8) f1, f2, f3, f4, g1, g2, g3, g4, fmin, fmax, dEmin, dEmax, E1, E2, E3, E4
    real(8) df1, df2, df3, df4, mu1
    integer i, N
    logical :: partial
    partial = .true.
    third = 1.0d0/3.0d0
    third2 = 2.0d0/3.0d0
    T_cur = T(NT)/g_kb  ! [eV]
    mu1 = mu(NT)     ! [eV]
    if (NT .EQ. 1) then
        dmu = (mu(2) - mu(1))/(T(2)-T(1))
    else
        dmu = (mu(NT) - mu(NT-1))/(T(NT)-T(NT-1))
    endif

    f = 1.0d0
    E = Mat_DOS%E(1)
    Ce_cur = 0.0d0
    ksc_cur = 0.0d0
    do while ((f .GT. 1d-12) .or. (f1 .GT. 1d-12) )
        dE = 0.1d0
        call adjust_dE(T_cur, mu1, dmu, E, Mat_DOS, partial, dE) ! below
        
        ! Distributions:
        f = Fermi(T_cur, mu1, E)

        !f1 = Diff_Fermi_Te(T_cur, mu1, dmu, E, partial)*(E-mu(1))
        f1 = Diff_Fermi_Te(T_cur, mu1, dmu, E, partial)*(E-mu1)
        !f1 = Diff_Fermi_Te(T_cur, mu1, dmu, E, partial)*E
        E2 = E+third*dE
        !f2 = Diff_Fermi_Te(T_cur, mu1, dmu, E2, partial)*(E2-mu(1))
        f2 = Diff_Fermi_Te(T_cur, mu1, dmu, E2, partial)*(E2-mu1)
        !f2 = Diff_Fermi_Te(T_cur, mu1, dmu, E2, partial)*E2
        E3 = E+third2*dE
        !f3 = Diff_Fermi_Te(T_cur, mu1, dmu, E3, partial)*(E3-mu(1))
        f3 = Diff_Fermi_Te(T_cur, mu1, dmu, E3, partial)*(E3-mu1)
        !f3 = Diff_Fermi_Te(T_cur, mu1, dmu, E3, partial)*E3
        E4 = E+dE
        !f4 = Diff_Fermi_Te(T_cur, mu1, dmu, E4, partial)*(E4-mu(1))
        f4 = Diff_Fermi_Te(T_cur, mu1, dmu, E4, partial)*(E4-mu1)
        !f4 = Diff_Fermi_Te(T_cur, mu1, dmu, E4, partial)*E4

        ! DOS:
        call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E, g1)
        call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E2, g2)
        call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E3, g3)
        call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E4, g4)
        ! Heat capacity:
        Ce_cur = Ce_cur + dE/8.0d0*(f1*g1 + 3.0d0*(f2*g2+f3*g3) + f4*g4)    ! Simpsone 3/8 rule
        
        ! Screening:
        df1 = Diff_Fermi(T_cur, mu1, E)
        df2 = Diff_Fermi(T_cur, mu1, E2)
        df3 = Diff_Fermi(T_cur, mu1, E3)
        df4 = Diff_Fermi(T_cur, mu1, E4)
        ! Correct:
        ksc_cur = ksc_cur + dE/8.0d0*(df1*g1 + 3.0d0*(df2*g2+df3*g3) + df4*g4)
        
        ! Make a step:
        E = E + dE
        !write(*,'(f,f,f,f,f,f)') E, dE, g1, f1, T_cur, dmu
    enddo
    !pause 'Ce'
end subroutine



subroutine adjust_dE(T_cur, mu_cur, dmu, E, Mat_DOS, partial, dE, deriv)
   real(8), intent(in) :: T_cur, mu_cur, dmu, E ! temperature, checmical potential, its step, current energy, its default step
   type(Density_of_states), intent(in) :: Mat_DOS  ! materail DOS
   logical, intent(in) :: partial ! full derivative or only partial
   logical, intent(in), optional :: deriv ! Fermi or its derivative
   real(8), intent(inout) :: dE ! energy step [eV]
   !-----------------------------------------------
   real(8) :: f0, f1, dEmax, dEmin, fmin, fmax, g0, g1, Func0, Func1, F0min
   logical :: do_deriv ! Fermi or its derivative
   if (present(deriv)) then
      do_deriv = deriv
   else ! be default, do derivative
      do_deriv = .true.
   endif
   ! Define the step limits:
    fmin = 1.0d-5
    fmax = 1.0d-2
    F0min = 1.0d-8
    dEmin = 1.0d-4
    dEmax = 0.01d0
   
   ! Distribution function:
   if (do_deriv) then
      f0 = Diff_Fermi_Te(T_cur, mu_cur, dmu, E, partial)*E
   else
      f0 = Fermi(T_cur, mu_cur, E)
   endif
   
   if (do_deriv) then
      f1 = Diff_Fermi_Te(T_cur, mu_cur, dmu, E+dE, partial)*(E+dE)
   else
      f1 = Fermi(T_cur, mu_cur, E+dE)
   endif
   
   ! DOS:
   call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E, g0) ! DOS at point E
   call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E+dE, g1) ! DOS at point E+dE
   ! The function we are integrating:
   Func0 = f0*g0
   Func1 = f1*g1

   if ((Func0 <= 1.0d-12) .and. (Func1 <= 1.0d-12)) goto 2017 ! zero populations, no need to do anything
   if (Func0 <= 1.0d-12) Func0 = Func1
   
   do while ((abs(Func1-Func0)/abs(Func0) .GT. fmax) .AND. (abs(f1) .GT. F0min))
      dE = dE/2.0d0
      if (do_deriv) then
         f1 = Diff_Fermi_Te(T_cur, mu_cur, dmu, E+dE, partial)*(E+dE)
      else
         f1 = Fermi(T_cur, mu_cur, E+dE)
      endif
      call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E+dE, g1) ! DOS at point E+dE
      Func1 = f1*g1
      if (dE .LT. dEmin) dE = dEmin
      if (dE .LE. dEmin) exit
!       write(*,'(a,f,e,e,e,e)') 'a', E, dE, Func1, Func0, abs(Func1-Func0)/abs(Func0)
   enddo
   do while (abs(Func1-Func0)/abs(Func0) .LT. fmin)
      dE = dE*2.0d0
      if (do_deriv) then
         f1 = Diff_Fermi_Te(T_cur, mu_cur, dmu, E+dE, partial)*(E+dE)
      else
         f1 = Fermi(T_cur, mu_cur, E+dE)
      endif
      call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E+dE, g1) ! DOS at point E+dE
      Func1 = f1*g1
      if (dE .GT. dEmax) dE = dEmax
      if (dE .GE. dEmax) exit
!       print*, 'b', E, dE
   enddo
   
2017 continue
!     print*, 'adjust_dE is done', E, dE
end subroutine adjust_dE


END MODULE Temperature_stuff