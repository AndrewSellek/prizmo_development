module prizmo_rates
  use prizmo_commons
  use prizmo_shielding
  use prizmo_utils
contains

  ! ************************
  subroutine compute_rates(x, Tgas, Tdust)
    implicit none
    real*8,intent(in)::x(nspecies), Tgas, Tdust
    real*8::Tgas32, Tgas14, invTgas, invTgas32, sqrTgas, sqrTgas32
    real*8::Tgas32_m06, Tgas32_m05, Tgas32_m04, Tgas32_m03, invTgas14_065
    real*8::sticking, nu_debye
    real*8::n_h, ione_h

    !! PREPROCESS_PROTOTYPES_DEFINE
    !! PREPROCESS_END

    ! temperature shortcuts
    Tgas32 = Tgas / 3d2
    Tgas14 = Tgas / 1d4
    invTgas = 1d0 / Tgas
    invTgas32 = 1d0 / Tgas32
    sqrTgas = sqrt(Tgas)
    sqrTgas32 = sqrt(Tgas32)
    Tgas32_m06 = Tgas32**(-0.6)
    Tgas32_m05 = Tgas32**(-0.5)
    Tgas32_m04 = Tgas32**(-0.4)
    Tgas32_m03 = Tgas32**(-0.3)
    invTgas14_065 = 1d0 / Tgas14**0.65

    ! dust parameters
    sticking = 1d0
    nu_debye = 1d12  ! 1/s

    ! secondary ionisation
    n_h = get_Hnuclei(x)
    ione_h = 13.60*ev2erg

    !! PREPROCESS_PROTOTYPES
    !! PREPROCESS_END

    !! PREPROCESS_RATES
    !! PREPROCESS_END

  end subroutine compute_rates

end module prizmo_rates
