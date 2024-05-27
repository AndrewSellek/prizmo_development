module prizmo_heating_photo
  use prizmo_commons
contains

  ! ***************
  function heating_photo(x, Tgas, Tdust, jflux, ntot) result(heat)
    implicit none
    real*8,intent(in)::x(nspecies), Tgas, Tdust, jflux(nphoto), ntot
    real*8::heat

    heat = 0d0

    !! PREPROCESS_PHOTOHEATING
    !! PREPROCESS_END

    secondion = fLoss_ion * heat
    heat = heat * max(1d0-fLoss_ion-fLoss_rad, 0d0)

  end function heating_photo

end module prizmo_heating_photo
