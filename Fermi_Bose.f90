!1111111111111111111111111111111111111111111111111111111111111
! This module contains Fermi and Bose function definitions, and related subroutines

MODULE Fermi_Bose
use Universal_Constants

implicit none

 contains

function Fermi(Te, mu, E)    ! Fermi-function
    real(8), intent(in) :: Te, mu, E   ! [eV], temperature, chem.potential, and energy
    real(8) :: Fermi   ! Fermi-function
    if ((E - mu)/Te >= log(HUGE(mu))) then ! exp(x) -> infinity
       Fermi = 0.0d0
    else
       Fermi = 1.0d0/(1.0d0 + dexp((E - mu)/Te))
    endif
end function

function Bose(Ta, E) ! Bose-function
    real(8), intent(in) :: Ta, E   ! [K], [eV], temperature and energy
    real(8) :: Bose   ! Bose-function
    Bose = 1.0d0/(-1.0d0 + dexp(E/(Ta/g_kb)))
end function

function Debye_phonon_DOS(Vs, Eph)
    real(8) Vs, Eph ! [m/s] speed of sound, [eV] energy state
    real(8) Debye_phonon_DOS
    !Debye_phonon_DOS = g_e*g_e*Eph*Eph/(2.0d0*g_Pi*g_Pi*g_h*g_h*g_h*Vs*Vs*Vs)
    Debye_phonon_DOS = 3.0d0*g_e*g_e*Eph*Eph/(2.0d0*g_Pi*g_Pi*g_h*g_h*g_h*Vs*Vs*Vs)
end function

function Diff_Fermi(Te, mu, E)    ! Derivative of the Fermi-function over E
    real(8), intent(in) :: Te, mu, E   ! [eV], temperature, chem.potential, and energy
    real(8) :: Diff_Fermi   ! Derivative of the Fermi-function
    real(8) F, Ex
    F = Fermi(Te, mu, E)
    if (F < 1.0d-10) then
       Diff_Fermi = 0.0d0
    else
       Ex = (E - mu)/Te
       if (abs(Ex) >= log(HUGE(mu))) then ! exp(x) -> infinity ! too small or too large exp(Ex)
          Diff_Fermi = 0.0d0
       else
          Diff_Fermi = -dexp(Ex)/Te*F*F
       endif
    endif
end function

function Diff_Fermi_Te(Te, mu, dmu, E, partial)
    real(8), intent(in) :: Te, mu, dmu, E   ! [eV], temperature, chem.potential, and energy
    real(8) :: Diff_Fermi_Te   ! Derivative of the Fermi-function
    logical, intent(in), optional :: partial ! partial derivative, or a full one?
    real(8) F, dmu_part
    F = Fermi(Te, mu, E)
    if (F <= 1.0d-20) then
       Diff_Fermi_Te = 0.0d0
    else
       !F = 1.0d0/(1.0d0 + dexp((E - mu)/Te))
       if (present(partial)) then
          if (partial) then
             dmu_part = 0.0d0
          else
             dmu_part = dmu
          endif
       else
          dmu_part = dmu
       endif
       Diff_Fermi_Te = dexp((E - mu)/Te)*F*F*(E - mu + Te*dmu_part)/(Te*Te)
    endif
!     if (Diff_Fermi_Te < 0.0d0) Diff_Fermi_Te = 0.0d0 ! TEST
end function



END MODULE Fermi_Bose
