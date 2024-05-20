module prizmo_attenuate
  use prizmo_commons
contains

  ! *******************
  subroutine attenuate(x, Tgas, jflux, ds)
    use prizmo_utils
    implicit none
    real*8,intent(in)::x(nspecies), Tgas, ds
    real*8,intent(inout)::jflux(nphoto)
    real*8::tau(nphoto), rhogas, nH

    tau = 0d0

    !! PREPROCESS_ATTENUATE
    !! PREPROCESS_END

    rhogas = get_rho(x)
    tau = tau + dust_kappa_opacity * rhogas * d2g

    nH = get_Hnuclei(x)
    tau = tau + PAH_sigma_opacity * nH * fPAH * XCPAH

    jflux = jflux * exp(-tau * ds)

  end subroutine

end module prizmo_attenuate
