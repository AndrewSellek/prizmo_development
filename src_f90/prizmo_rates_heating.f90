module prizmo_rates_heating
  use prizmo_commons
  use prizmo_utils
contains

  ! ************************
  subroutine compute_photorates_heating(jflux)
    implicit none
    real*8,intent(in)::jflux(nphoto)
    real*8::f(nphoto), kernel(nphoto)
    real*8::Ephotoelectron(nphoto)
    real*8::loss_ion(nphoto), loss_rad(nphoto)

    kall_heat = 0d0
    kall_secondIon = 0d0

    kernel = jflux / energy / hplanck

    !!! EXAMPLE
    !! H -> H+ + E
    !Ephotoelectron = max(energy - energy_threshold(272), 0d0)
    !loss_ion = calc_loss_ion_E(Ephotoelectron)
    !loss_rad = calc_loss_rad_E(Ephotoelectron)
    !loss_ion = lossH2(1) * (sign(0.5d0,Ephotoelectron-EthH2(1)) + 0.5d0) * xH2 + lossH(1) * (sign(0.5d0,Ephotoelectron-EthH(1)) + 0.5d0) * (1d0-xH2)
    !loss_ion = max(min( loss_ion, 1d0), 0d0)
    !loss_rad = (lossH2(2) * (sign(0.5d0,Ephotoelectron-EthH2(2)) + 0.5d0)  + lossH2(3) * (sign(0.5d0,Ephotoelectron-EthH2(3)) + 0.5d0)) * xH2 + lossH(2) * (sign(0.5d0,Ephotoelectron-EthH(2)) + 0.5d0) * (1d0-xH2)
    !loss_rad = max(min( loss_rad, 1d0), 0d0)
    !f(:) = photo_xsecs(:, 272) * Ephotoelectron * kernel * max(1d0-loss_ion-loss_rad, 0d0)
    !kall_heat(272) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0
    !f(:) = photo_xsecs(:, 272) * Ephotoelectron * kernel * loss_ion
    !kall_secondIon(272) = sum((f(2:nphoto) + f(1:nphoto-1)) * delta_energy) / 2d0

    !! PREPROCESS_PHOTOHEATING_RATE
    !! PREPROCESS_END

  end subroutine compute_photorates_heating


  ! **********************
  function calc_loss_ion_E(Ephotoelectron) result(loss_ion)
    use prizmo_commons
    implicit none
    real*8,intent(in)::Ephotoelectron(nphoto)
    real*8::loss_ion(nspecies)

    loss_ion = lossH2(1) * (sign(0.5d0,Ephotoelectron-EthH2(1)) + 0.5d0) * xH2 + lossH(1) * (sign(0.5d0,Ephotoelectron-EthH(1)) + 0.5d0) * (1d0-xH2)
    loss_ion = max(min( loss_ion, 1d0), 0d0)

  end function calc_loss_ion_E


  ! **********************
  function calc_loss_rad_E(Ephotoelectron) result(loss_rad)
    use prizmo_commons
    implicit none
    real*8,intent(in)::Ephotoelectron(nphoto)
    real*8::loss_rad(nspecies)

    loss_rad = lossH2(1) * (sign(0.5d0,Ephotoelectron-EthH2(1)) + 0.5d0) * xH2 + lossH(1) * (sign(0.5d0,Ephotoelectron-EthH(1)) + 0.5d0) * (1d0-xH2)
    loss_rad = max(min( loss_ion, 1d0), 0d0)

  end function calc_loss_rad_E


  ! ************************
  subroutine compute_secondary_loss(x)
    implicit none
    real*8,intent(in)::x(nspecies)
    real*8::xe !,xH2
    real*8::f0H(2),f0H2(3),alphaH2(3),betaH(2),betaH2(3),gammaH(2)
    !real*8::ionH,ionH2,radH,radH2LW,radH2nu

    xe  = min(x(idx_E)/get_Hnuclei(x),1d0)
    xH2 = 2d0*x(idx_H2)/get_Hnuclei(x)

    f0H     = (/0.4412, 0.4766/)
    f0H2    = (/0.4221, 0.1913, 0.0755/)
    alphaH2 = (/6.72, 7.0, 23500.0/)
    betaH   = (/0.4092, 0.2735/)
    betaH2  = (/0.824, 0.8, 0.955/)
    gammaH  = (/1.7592, 1.5221/)
    !EthH    = (/13.60, 10.2/)
    !EthH2   = (/15.12, 11.37, 0.516/)

    !ionH  = 0.4412 * (1-xe**0.4092)**1.7592
    !radH  = 0.4766 * (1-xe**0.2735)**1.5221
    lossH = f0H * (1d0 - xe**betaH)**gammaH
    
    !ionH2 = 0.4211 / (1+6.72*xe**0.824)
    !radH2LW = 0.1913 / (1+7d0*xe**0.8)
    !radH2nu = 0.0755 / (1+23500d0*xe**0.955)
    lossH2 = f0H2 / (1d0 + alphaH2 * xe**betaH2)

    ! Save limiting values
    fLoss_ion = lossH2(1)*xH2 + lossH(1)*(1d0-xH2)
    fLoss_ion = max(min( fLoss_ion, 1d0), 0d0)
    fLoss_rad = (lossH2(2)+lossH2(3))*xH2 + lossH(2)*(1d0-xH2)
    fLoss_rad = max(min( fLoss_rad, 1d0), 0d0)

  end subroutine compute_secondary_loss

end module prizmo_rates_heating
