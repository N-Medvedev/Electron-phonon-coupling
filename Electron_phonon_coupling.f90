!1111111111111111111111111111111111111111111111111111111111111
! This module contains calculations of he electron-phonon coupling

MODULE Electron_phonon_coupling
use Universal_Constants
use Objects
use Fermi_Bose
use Dealing_with_arrays
implicit none


 contains

subroutine Electron_phonon(Te_grid, mu_grid, Matter, Mat_DOS, ksc, Error_message, gfactor)  ! coupling parameter
    real(8), dimension(:), intent(in) :: Te_grid, mu_grid    ! [K], [eV]
    real(8), dimension(:), allocatable, intent(out) :: gfactor    ! [J/s/m^3/K]
    type(Solid), intent(in) ::  Matter
    type(Density_of_states), intent(in) :: Mat_DOS  ! materail DOS
    real(8), dimension(:), intent(in) :: ksc    ! screening [1/m]
    type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
    real(8), dimension(:), allocatable :: g_ch    ! [J/s/m^3/K] temporary
    real(8) T, mu_temp, g
    integer i, N, my_id, OMP_GET_THREAD_NUM
    N = size(Te_grid)
    if (.not. allocated(gfactor)) allocate(gfactor(N))
    gfactor = 0.0d0
    if (.not. allocated(g_ch)) allocate(g_ch(N))
    g_ch = 0.0d0
    !print*, 'Start of the parallel region'
    !$omp parallel &
    !$omp private (i, my_id, T, mu_temp, g)
    !$omp do reduction( + : g_ch)
    do i = 1, N ! all gridpoints
        T = Te_grid(i)/g_kb  ! [eV]
        mu_temp = mu_grid(i) ! [eV]
        call get_gfactor(T, Matter, Mat_DOS, mu_temp, Te_grid, ksc, g, Error_message, mu_grid(1))
        g_ch(i) = g  ! save it in the array
        my_id = 1 + OMP_GET_THREAD_NUM() ! identify which thread it is
        print*, Te_grid(i), Matter%Ta, g,  'G'
        !pause 'PAUSE'
    enddo
    !$omp end do
    !$omp end parallel
    gfactor = g_ch  ! save it
    !print*, gfactor(:)
    if (allocated(g_ch)) deallocate(g_ch)
end subroutine

subroutine get_gfactor(T, Matter, Mat_DOS, mu, T_grid, ksc, g, Error_message, Efermi)
    real(8), intent(in) :: T, mu    ! [eV] electronic temperature
    real(8), intent(out) :: g       ! [J/s/m^3/K]
    type(Solid), intent(in) ::  Matter
    type(Density_of_states), intent(in) :: Mat_DOS  ! materail DOS
    real(8), dimension(:), intent(in) :: T_grid, ksc    ! screening [1/m]
    type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
    real(8), intent(in) :: Efermi	! [eV] Fermi level (just for test calculations)
    real(8) Ta, dEdt
    integer i
    Ta = Matter%Ta    ! [K] lattice temperature
    !Calculate:
    if (abs(T*g_kb - Ta) < 1.0d-6) then	! too close temperatures, shift a little:
       call get_dEdt(T+1.0d0, Matter, Mat_DOS, mu, T_grid, ksc, dEdt, Error_message, Efermi)
       g = dEdt/(1.0d0) !*Matter%At_dens*1d6*Matter%Ne
    else	! temperatures are far enough
       call get_dEdt(T, Matter, Mat_DOS, mu, T_grid, ksc, dEdt, Error_message, Efermi)
       g = dEdt/(T*g_kb - Ta) !*Matter%At_dens*1d6*Matter%Ne
    endif
end subroutine


subroutine get_dEdt(T, Matter, Mat_DOS, mu, T_grid, ksc, dEdt, Error_message, Efermi)
    real(8), intent(in), target :: T, mu    ! [eV], [eV]
    real(8), intent(out) :: dEdt    ! [J/s]
    type(Solid), intent(in) ::  Matter
    type(Density_of_states), intent(in) :: Mat_DOS  ! materail DOS
    real(8), dimension(:), intent(in) :: T_grid, ksc    ! screening [1/m]
    type(Error_handling), intent(inout) :: Error_message  ! save data about error if any
    real(8), intent(in) :: Efermi	! [eV] Fermi level (just for test calculations)
    !---------------------------------------------
    real(8), pointer :: T_cur, mu_cur
    real(8) M2, Ta, dEdt_cur, third, third2, E, dE, f, f0, f1, fmin, fmax, dEmin, dEmax
    real(8) E2, E3, E4, Ieph1, Ieph2, Ieph3, Ieph4, g1, g2, g3, g4
    integer i
    T_cur => T
    mu_cur => mu
    third = 1.0d0/3.0d0
    third2 = 2.0d0/3.0d0
    dEdt_cur = 0.0d0
    fmin = 1.d-6
    fmax = 1.d-4
    ! Correct grid:
    dEmin = 1.0d-3
    ! Testing incorrect grid:
    dEmin = 0.01d0

    dEmax = 0.1d0
    
    !E = 2.0d0*Matter%E_debye
    E = Mat_DOS%E(1)	! starting point
    dE = 0.1d0
    f = 1.0d0
    f0 = 1.0d0
    do while (f .GT. 1d-30)
        ! Correct grid:
        f1 = Fermi(T_cur, mu_cur, E+dE)
        ! Correct values:
        call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E, g1)
        call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E+dE, g4)
        do while ((abs(f1*g4-f0*g1) .GT. fmax) .AND. (abs(f1*g4) .GT. fmin))
            dE = dE/2.0d0
            f1 = Fermi(T_cur, mu_cur, E+dE) ! fermi-function
            call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E+dE, g4)
            if (dE .LT. dEmin) dE = dEmin
            if (dE .LE. dEmin) exit
        enddo
        f1 = Fermi(T_cur, mu_cur, E+dE) ! fermi-function
        call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E+dE, g4)
        do while ((abs(f1*g4-f0*g1) .LT. fmin) .or. (abs(f1*g4) .LT. fmin))
            dE = dE*2.0d0
            f1 = Fermi(T_cur, mu_cur, E+dE) ! fermi-function
            call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E+dE, g4)
            if (dE .GT. dEmax) dE = dEmax
            if (dE .GE. dEmax) exit
        enddo
        ! Testing incorrect grid:
!         dE = dEmin
        
        E2 = E+third*dE
        E3 = E+third2*dE
        E4 = E+dE
        f0 = Fermi(T_cur, mu_cur, E4)
        f = f0
        ! Correct values:
        call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E, g1)
        call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E2, g2)
        call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E3, g3)
        call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E4, g4)
        ! Testing incorrect expression from Zhigilei:
!         call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, Efermi, g1)
!         g2 = g1
!         g3 = g1
!         g4 = g1
        
        call Collision_integral(T, T_grid, Matter, Mat_DOS, mu, ksc, E, Ieph1, Efermi)
        call Collision_integral(T, T_grid, Matter, Mat_DOS, mu, ksc, E2, Ieph2, Efermi)
        call Collision_integral(T, T_grid, Matter, Mat_DOS, mu, ksc, E3, Ieph3, Efermi)
        call Collision_integral(T, T_grid, Matter, Mat_DOS, mu, ksc, E4, Ieph4, Efermi)
        dEdt_cur = dEdt_cur + dE/8.0d0*(Ieph1*g1*E + 3.0d0*(Ieph2*g2*E2+Ieph3*g3*E3) + Ieph4*g4*E4)
        E = E + dE
        !print*, 'dEdt', E, dEdt_cur*g_e*(Matter%At_dens*1d6)
    enddo
    dEdt = -dEdt_cur*g_e*(Matter%At_dens*1d6) ! [J/s/m^3]
    nullify(T_cur)
    nullify(mu_cur)
end subroutine


subroutine Collision_integral(T, T_grid, Matter, Mat_DOS, mu, ksc, E, Ieph, Efermi)
    real(8), intent(in) :: T, mu    ! [eV], [eV]
    type(Solid), intent(in) ::  Matter
    type(Density_of_states), intent(in) :: Mat_DOS  ! materail DOS
    real(8), dimension(:), intent(in) :: T_grid, ksc    ! [K] temperature grid, [1/m] inverse screening length
    real(8), intent(in) :: E    ! energy of electron
    real(8), intent(out) :: Ieph    ! electron-phonon scattering integral
    real(8), intent(in) :: Efermi	! [eV] Fermi level (just for test calculations)
    !---------------------------------------------
    real(8) ke, Int_Ed, Eph, dEph, third, third2, thet, thet2, fmin, fmax, dEmin, dEmax
    real(8) fe, fe1(4), fe2(4), fa(4), fa2(4), f0, f1, f2, Mel2(4), q(4), DOSph(4), DOSe(4), ke1(4), Ftil(4), DOSe2(4), ke2(4), g0, temp
    integer i
    logical commented
    !commented = .false.  ! printouts are not commented out
    commented = .true.  ! printouts are commented out
    third = 1.0d0/3.0d0
    third2 = 2.0d0/3.0d0
    
    ! Correct expression:
    call Linear_approx(Mat_DOS%E, Mat_DOS%k, E, ke)
    ! Testing incorrect expression fromZhigilei:
!     call Linear_approx(Mat_DOS%E, Mat_DOS%k, Efermi, ke)
  
    fe = Fermi(T, mu, E)
    f0 = fe
    call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E, g0)
    
    Int_Ed = 0.0d0
    fmin = 1.d-6
    fmax = 1.d-4
    ! Correct:
    dEmin = Matter%E_debye*5.0d-4
    ! Testing incorrect grid:
    dEmin = Matter%E_debye*1.0d-3
    
!      print*, 'Matter%E_debye=', Matter%E_debye
!      pause 
    
    dEmax = Matter%E_debye/100.0d0
    dEph = dEmax
    Eph = dEph/2.0d0
    !if (.not. commented) print*, 'Matter%E_debye', Matter%E_debye
    do while (Eph .LT. Matter%E_debye)  ! integrate all phonon energies
        !if (.not. commented) print*, 'null', Eph, dEph, Int_Ed
        f1 = Fermi(T, mu, E+(Eph+dEph))
        f2 = Fermi(T, mu, E-(Eph+dEph))
        call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E+Eph+dEph, DOSe(4))
        call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E-(Eph+dEph), DOSe2(4))
        
        ! Correct:
        do while (((abs(f1*DOSe(4)-f0*g0) .GT. fmax) .AND. (abs(f1*DOSe(4)) .GT. fmin)) .AND. ((abs(f2*DOSe2(4)-f0*g0) .GT. fmax) .AND. (abs(f2*DOSe2(4)) .GT. fmin)))
            dEph = dEph/2.0d0
            f1 = Fermi(T, mu, E+(Eph+dEph)) ! fermi-function
            f2 = Fermi(T, mu, E-(Eph+dEph)) ! fermi-function
            call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E+Eph+dEph, DOSe(4))
            call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E-(Eph+dEph), DOSe2(4))
            if (dEph .LT. dEmin) dEph = dEmin
            if (dEph .LE. dEmin) exit
        enddo
        f1 = Fermi(T, mu, E+(Eph+dEph)) ! fermi-function
        f2 = Fermi(T, mu, E-(Eph+dEph)) ! fermi-function
        call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E+Eph+dEph, DOSe(4))
        call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E-(Eph+dEph), DOSe2(4))
        do while (((abs(f1*DOSe(4)-f0*g0) .LT. fmin) .or. (max(f1*DOSe(4),f0*g0) .LT. fmin)) .AND. ((abs(f2*DOSe2(4)-f0*g0) .LT. fmin) .or. (max(f2*DOSe2(4),f0*g0) .LT. fmin)))
            dEph = dEph*2.0d0
            f1 = Fermi(T, mu, E+(Eph+dEph)) ! fermi-function
            f2 = Fermi(T, mu, E-(Eph+dEph)) ! fermi-function
            call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E+Eph+dEph, DOSe(4))
            call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E-(Eph+dEph), DOSe2(4))
            if (dEph .GT. dEmax) dEph = dEmax
            if (dEph .GE. dEmax) exit
        enddo
        ! Testing incorrect grid:
!          dEph = dEmax
        
        !if (.not. commented) print*, 'raz', Eph, dEph, Int_Ed
        
        q(4) = (Eph+dEph)*g_e/(g_h*Matter%Vs) ! phonon wave-vector
        ! Correct expressions:
        call Linear_approx(Mat_DOS%E, Mat_DOS%k, E+Eph+dEph, ke1(4))
        call Linear_approx(Mat_DOS%E, Mat_DOS%k, E-(Eph+dEph), ke2(4))
        ! Testing incorrect expressions from Zhigilei:
!         call Linear_approx(Mat_DOS%E, Mat_DOS%k, Efermi, ke1(4))
!         ke2(4) = ke1(4)
        
        ! Angles:
        thet = Theta(ke, ke1(4), q(4))  ! theta-function for incoming
        thet2 = Theta(ke, ke2(4), q(4))  ! theta-function for outgoing
        
        q(1) = Eph*g_e/(g_h*Matter%Vs) ! phonon wave-vector
        q(2) = (Eph+third*dEph)*g_e/(g_h*Matter%Vs) ! phonon wave-vector
        q(3) = (Eph+third2*dEph)*g_e/(g_h*Matter%Vs) ! phonon wave-vector
        !q(4) is done above
        
        Mel2(1) = M2(T,q(1),T_grid,ksc,Matter)   ! matrix element squared
        Mel2(2) = M2(T,q(2),T_grid,ksc,Matter)   ! matrix element squared
        Mel2(3) = M2(T,q(3),T_grid,ksc,Matter)   ! matrix element squared
        Mel2(4) = M2(T,q(4),T_grid,ksc,Matter)   ! matrix element squared

        DOSph(1) = Debye_phonon_DOS(Matter%Vs, Eph)    ! phonon DOS
        DOSph(2) = Debye_phonon_DOS(Matter%Vs, (Eph+third*dEph))    ! phonon DOS
        DOSph(3) = Debye_phonon_DOS(Matter%Vs, (Eph+third2*dEph))    ! phonon DOS
        DOSph(4) = Debye_phonon_DOS(Matter%Vs, (Eph+dEph))    ! phonon DOS
        
        ! Correct expressions:
        call Linear_approx(Mat_DOS%E, Mat_DOS%k, E+Eph, ke1(1))
        call Linear_approx(Mat_DOS%E, Mat_DOS%k, E+Eph+third*dEph, ke1(2))
        call Linear_approx(Mat_DOS%E, Mat_DOS%k, E+Eph+third2*dEph, ke1(3))
        !ke1(4) is done above
        ! Testing incorrect expressions from Zhigilei:
!         call Linear_approx(Mat_DOS%E, Mat_DOS%k, Efermi, ke1(1))
!         ke1(2) = ke1(1)
!         ke1(3) = ke1(1)
        
        ! Correct expressions:
        call Linear_approx(Mat_DOS%E, Mat_DOS%k, E-Eph, ke2(1))
        call Linear_approx(Mat_DOS%E, Mat_DOS%k, E-(Eph+third*dEph), ke2(2))
        call Linear_approx(Mat_DOS%E, Mat_DOS%k, E-(Eph+third2*dEph), ke2(3))
        !ke2(4) is done above
        ! Testing incorrect expressions from Zhigilei:
!         call Linear_approx(Mat_DOS%E, Mat_DOS%k, Efermi, ke2(1))
!         ke2(2) = ke2(1)
!         ke2(3) = ke2(1)
        
        call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E+Eph, DOSe(1))
        call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E+Eph+third*dEph, DOSe(2))
        call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E+Eph+third2*dEph, DOSe(3))
        call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E+Eph+dEph, DOSe(4))
        
        call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E-Eph, DOSe2(1))
        call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E-(Eph+third*dEph), DOSe2(2))
        call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E-(Eph+third2*dEph), DOSe2(3))
        call Linear_approx(Mat_DOS%E, Mat_DOS%DOS, E-(Eph+dEph), DOSe2(4))
        
        fe1(1) = Fermi(T, mu, E+Eph)
        fe1(2) = Fermi(T, mu, E+(Eph+third*dEph))
        fe1(3) = Fermi(T, mu, E+(Eph+third2*dEph))
        fe1(4) = f1
        
        fe2(1) = Fermi(T, mu, E-Eph)
        fe2(2) = Fermi(T, mu, E-(Eph+third*dEph))
        fe2(3) = Fermi(T, mu, E-(Eph+third2*dEph))
        fe2(4) = f2
        
        fa(1) = Bose(Matter%Ta, Eph)
        fa(2) = Bose(Matter%Ta, Eph+third*dEph)
        fa(3) = Bose(Matter%Ta, Eph+third2*dEph)
        fa(4) = Bose(Matter%Ta, Eph+dEph)
        
        ! Correct:
        Ftil(:) = thet*DOSe(:)/ke1(:)*(fe1(:)*(1.0d0 - fe)*(fa(:) + 1.0d0) - fe*(1.0d0 - fe1(:))*fa(:)) + &
                thet2*DOSe2(:)/ke2(:)*(fe2(:)*(1.0d0 - fe)*(fa(:)) - fe*(1.0d0 - fe2(:))*(fa(:) + 1.0d0))
        ! Allen's factor (also correct!):
!         fa2(1) = Bose(T*g_kb, Eph)
!         fa2(2) = Bose(T*g_kb, Eph+third*dEph)
!         fa2(3) = Bose(T*g_kb, Eph+third2*dEph)
!         fa2(4) = Bose(T*g_kb, Eph+dEph)
!         Ftil(:) = (fa(:) - fa2(:)) * (thet*DOSe(:)/ke1(:) * (fe1(:) - fe) - thet2*DOSe2(:)/ke2(:) * (fe - fe2(:)))
        ! Zhigilei's approximation (also acceptable!):
!         Ftil(1) = Eph * (fa(1) - fa2(1)) * (Diff_Fermi(T, mu, E+Eph/2.0d0)*thet*DOSe(1)/ke1(1) - Diff_Fermi(T, mu, E-Eph/2.0d0)*thet2*DOSe2(1)/ke2(1))        
!         Ftil(2) = Eph * (fa(2) - fa2(2)) * (Diff_Fermi(T, mu, E+Eph/2.0d0+third*dEph)*thet*DOSe(2)/ke1(2) - Diff_Fermi(T, mu, E-(Eph/2.0d0+third*dEph))*thet2*DOSe2(2)/ke2(2))
!         Ftil(3) = Eph * (fa(3) - fa2(3)) * (Diff_Fermi(T, mu, E+Eph/2.0d0+third2*dEph)*thet*DOSe(3)/ke1(3) - Diff_Fermi(T, mu, E-(Eph/2.0d0+third2*dEph))*thet2*DOSe2(3)/ke2(3))
!         Ftil(4) = Eph * (fa(4) - fa2(4)) * (Diff_Fermi(T, mu, E+Eph/2.0d0+dEph)*thet*DOSe(4)/ke1(4) - Diff_Fermi(T, mu, E-(Eph/2.0d0+dEph))*thet2*DOSe2(4)/ke2(4))
        
        Int_Ed = Int_Ed + dEph/8.0d0*(Mel2(1)*DOSph(1)/q(1)*Ftil(1) + &
                 3.0d0*(Mel2(2)*DOSph(2)/q(2)*Ftil(2) + Mel2(3)*DOSph(3)/q(3)*Ftil(3)) + &
                 Mel2(4)*DOSph(4)/q(4)*Ftil(4))
        Eph = Eph + dEph
        if ((.not. commented) .AND. (Int_Ed .NE. 0.0d0) .OR. isnan(Int_Ed)) then
            write(*,'(a,es,es,es,es)') 'q', q(:)
            write(*,'(a,es,es,es,es)') 'Mel2', Mel2(:)
            write(*,'(a,es,es,es,es)') 'DOSph', DOSph(:)
            write(*,'(a,es,es,es,es)') 'ke1', ke1(:)
            write(*,'(a,es,es,es,es)') 'ke2', ke2(:)
            write(*,'(a,es,es,es,es)') 'DOSe', DOSe(:)
            write(*,'(a,es,es,es,es)') 'DOSe2', DOSe2(:)
            write(*,'(a,es,es,es,es)') 'fe1', fe1(:)
            write(*,'(a,es,es,es,es)') 'fe2', fe2(:)
            write(*,'(a,es,es,es,es)') 'fa', fa(:)
            write(*,'(a,es,es,es,es)') 'Ftil', Ftil(:)
            print*, 'dva', Eph, dEph, Int_Ed
            print*, 'theta 1 n 2:', thet, thet2
            !pause 'PAUSE DVA'
        endif
    enddo
    Ieph = g_Pi*g_Pi*g_Pi/(g_h*ke)*Int_Ed*(Matter%At_dens*1d6)
    !if (.not. commented) then
        !print*, 'Ieph = ', Ieph
        !pause 'PAUSE DVA'
    !endif
end subroutine


function Theta(k, k1, q)
    real(8) Theta, k, k1, q
    real(8) t
    t = (q*q+k*k-k1*k1)/(2.0d0*k*q)
    if (abs(t) .LE. 1.0d0) then
        Theta = 1.0d0
    else
        Theta = 0.0d0
    endif
end function


function M2(Te,q,T,ksc,Matter)
    real(8) :: M2   ! matrix element squared
    real(8), intent(in) :: Te   ! [K] current electron temperature
    real(8), intent(in) :: q    ! [1/m] phonon wave-vector
    real(8), dimension(:), intent(in) :: T, ksc    ! [K] temperature grid, [1/m] inverse screening length
    type(Solid), intent(in) :: Matter   ! material parameters (for speed of sound)
    real(8) k0  ! current screening value
    real(8) Eph ! phonon energy
    Eph = Matter%Vs*g_h*q   ! for accustic phonon [J]
    ! Correct screening parameter:
    call Linear_approx(T, ksc, Te, k0)  ! find screening k0 as approximation from ksc-array (ksc vs T)
    ! Testing Fermi screening:
!     k0 = ksc(1)

    M2 = g_e*g_e/(2.0d0*g_e0)*Eph/(q*q + k0*k0)
    !write(*,'(a,es,es)') 'M2-screening:', k0, q
end function


END MODULE Electron_phonon_coupling