from prizmo_commons import clight, hplanck, nphoto, kboltzmann, erg2ev, ev2erg, \
    hplanck_eV, energy_min, energy_max, print_title, photo_logspacing, natom2name, name2natom, sp2idx, idx2sp, \
    fuv_energy1, fuv_energy2, radiation_type, data_dir
from prizmo_preprocess import preprocess
import numpy as np
from glob import glob
from scipy.interpolate import interp1d
from scipy import integrate
import sys

# Calculates the cross-section according to Verner for outer and inner shells
def calc_crossSection_Verner(Z, Q, dN=1, maxdN=10, inner=True):
    maxdN = min(maxdN,10)
    # Line counting
    #linesOuter = {}
    #linesOuterCount = 0
    #for Z in natom2name.keys():
    #    linesOuter[Z] = linesOuterCount
    #    linesOuterCount += Z

    # Shell identification
    dataFor = ["H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","S","Ar","Ca","Fe"]
    dataLineOld = int(Z*(Z-1)/2)+Q
    dataLine = np.sum([name2natom[el] for el in dataFor[:dataFor.index(natom2name[Z])]], dtype=int)+Q
    ne = Z-Q
    shelle = [1,1,2,2,3,3,3,3,3,3,4,4,5,5,5,5,5,5,6,6,6,6,6,6,6,6,6,6,7,7]  # Shell that nth electron occupies
    shellName = {1: '1s', 2: '2s', 3: '2p', 4: '3s', 5: '3p', 6: '3d', 7: '4s'}
    shellOuter = shelle[ne-1]

    # Load Cross-section fit parameters
    outerData = np.genfromtxt('../data/photoionization_ions/ph2.dat', names=True)
    innerData = np.genfromtxt('../data/photoionization_ions/ph1b.dat', names=True)
    AugerMeitnerLookup = np.genfromtxt('../data/photoionization_ions/auger.dat', names=['Z','stage','shell','maxdN'], usecols=(0,1,2,3))
    AugerMeitnerYield = np.genfromtxt('../data/photoionization_ions/auger.dat', usecols=(4,5,6,7,8,9,10,11,12,13))#, names=['Z','stage','shell','maxdN']+['p{}'.format(n) for n in range(1,11)])

    # Initalise Energies and sigma
    E = np.geomspace(outerData['Eth'][dataLine], energy_max*erg2ev, 1000)
    print(np.min(E), np.max(E))
    sigma = np.zeros_like(E)

    # Outer Shell (Verner et al. 1996); outer shell cannot have an Auger-Meitner effect and therefore only applies for dN=1.
    if dN==1:
        sigma0 = outerData['sigma0'][dataLine]
        x = E/outerData['E0'][dataLine] - outerData['y0'][dataLine]
        y = np.sqrt(x**2+outerData['y1'][dataLine]**2)
        F = ((x-1)**2 + outerData['yw'][dataLine]**2) * np.power(y,0.5*outerData['P'][dataLine]-5.5) * (1 + np.sqrt(y/outerData['ya'][dataLine]))**(-outerData['P'][dataLine])
        sigma += 1e-18*sigma0*F

    # Inner Shells (Verner & Yakovlev 1995)
    if inner:
        dataLines = np.nonzero((innerData['Z']==Z) * (innerData['N']==ne))[0]
        for shell in range(1,shellOuter):
            dataLine = dataLines[-shell]
            sigma0 = innerData['sigma0'][dataLine]
            y = E/innerData['E0'][dataLine]
            F = ((y-1)**2 + innerData['yw'][dataLine]**2) * np.power(y,0.5*innerData['P'][dataLine]-5.5-innerData['l'][dataLine]) * (1 + np.sqrt(y/innerData['ya'][dataLine]))**(-innerData['P'][dataLine])
            AMYields = AugerMeitnerYield[(AugerMeitnerLookup['Z']==Z)*(AugerMeitnerLookup['stage']==Q+1)*(AugerMeitnerLookup['shell']==shell),:][0]
            TotYield = np.sum(AMYields[0:maxdN])
            AMYield = AMYields[dN-1]/TotYield
            sigma += 1e-18*sigma0*F*(E>innerData['Eth'][dataLine])*AMYield

    return E, sigma

# Photoionization of H2
def calc_branching_Chung(energy):
    # Calculates fraction of H2 photoionization events resulting in dissociative ionization (i.e. H2 -> H+ + H + E)
    E = [18.076, 18.1, 18.15, 18.2, 18.3, 18.4, 18.5, 18.6, 18.8, 19.0, 19.5, 20.0, 20.5, 21.0, 21.5, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 35.0, 35.5, 36.0, 36.5, 37.0, 37.5, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 48.0, 50.0, 52.0, 54.0, 56.0, 58.0, 60.0, 62.0, 64.0, 66.0, 68.0, 70.0, 72.0, 74.0, 76.0, 78.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 115.0, 120.0, 124]
    branch = [0, 0.0022, 0.0046, 0.0063, 0.0086, 0.0102, 0.0115, 0.0128, 0.0148, 0.0164, 0.0190, 0.0200, 0.0206, 0.0212, 0.0216, 0.0220, 0.0228, 0.0235, 0.0242, 0.0252, 0.0275, 0.0333, 0.0410, 0.0508, 0.065, 0.081, 0.096, 0.110, 0.121, 0.123, 0.117, 0.105, 0.101, 0.100, 0.102, 0.112, 0.122, 0.132, 0.142, 0.149, 0.152, 0.155, 0.158, 0.169, 0.182, 0.195, 0.208, 0.222, 0.235, 0.247, 0.256, 0.264, 0.271, 0.277, 0.282, 0.284, 0.285, 0.284, 0.283, 0.282, 0.279, 0.276, 0.273, 0.271, 0.268, 0.265, 0.262, 0.260, 0.258]

    interpor = interp1d(np.log10(E), np.log10(branch), fill_value=-np.inf, bounds_error=False)
    return 1e1**interpor(np.log10(energy))

def calc_crossSection_Yan(E=None, smooth=True, noLo=True):
    # H2 photoionization cross-section (total) according to Yan et al. 1998
    # Given in units of barn
    barn = 1e-24

    # Initialise energies
    if E is None:
        E = np.geomspace(15.4+2.6*noLo, energy_max*erg2ev, 1000)
    x = E/15.4
    
    # Piecewise functions
    def lo(x):
        return 1e8   * (-37.895 + 99.723*x - 87.227*x**2 + 25.400*x**3)
    def mi(E):
        return 2e7   * (0.071 - 0.673/x + 1.977/x**2 - 0.692/x**3) * x**-0.252
    hi = 45.57 * (1 - 2.003/x**0.5 - 4.806/x + 50.577/x**1.5 - 171.044/x**2 + 231.608/x**2.5 -81.885/x**3)/(E/1000)**3.5

    # Construct cross-section
    sigma = hi*(E>85)*barn
    if smooth:
        sigma += np.minimum(mi(x),hi)*(18<E)*(E<85)*barn
    else:
        sigma += mi(x)*(18<E)*(E<85)*barn
    if not noLo:
        sigma += lo(x)*(15.4<E)*(E<18)*barn        

    return E, sigma

# returns the energy bins in eV
def prepare(photo_limits, species):

    print_title("photochemistry")

    data_all = dict()
    # Leiden format
    for fname in glob("../data/xsecs/*.dat"):
        reactant_name = fname.split("/")[-1].split("__")[0]
        if sp2idx(reactant_name) not in species:
            continue
        header = [x for x in open(fname).read().split("\n") if x.strip().startswith("#")][-1]
        header = header.replace("#", "").split()
        data = np.loadtxt(fname).T
        data = {h: data[i][::-1] for i, h in enumerate(header)}
        data["energy"] = clight * hplanck / (data["wavelength"] * 1e-7)  # nm -> erg
        data["file"] = fname

        name_file = fname.split("/")[-1].replace(".dat", "")
        data_all[name_file] = data

    def get_plus(val):
        if val > 200:
            plus = "+" + str(val)
        else:
            plus = "+"*val
        return plus

    for idx in species:
        spec=idx2sp(idx)
        if spec=='H2' and 'H2__H2+_E' in data_all.keys():
            Eev = data_all['H2__H2+_E']['energy']*erg2ev
            _, data_all['H2__H2+_E']['photoionisation'][(Eev>18)] = calc_crossSection_Yan(Eev[(Eev>18)])
            Emax = np.max(Eev)
            Eext, PIext = calc_crossSection_Yan()
            Etot = np.append(Eext[(Eext>Emax)][::-1]*ev2erg,data_all['H2__H2+_E']['energy'])
            PItot = np.append(PIext[(Eext>Emax)][::-1],data_all['H2__H2+_E']['photoionisation'])
            PDtot = np.append(np.zeros_like(PIext[(Eext>Emax)][::-1]),data_all['H2__H2+_E']['photodissociation'])
            data_all['H2__H2+_E']['energy'] = Etot
            data_all['H2__H2+_E']['photoionisation'] = PItot
            data_all['H2__H2+_E']['photodissociation'] = PDtot
            continue
        if not spec.replace('+','') in natom2name.values():
            continue

        Z = name2natom[spec.replace('+','')]
        Q = spec.count('+')
        if Q>=Z:
            continue
        assert Q<Z, "Can't have charge >= atomic number, or ionization when charge == atomic number"
        N = Z-Q

        atom = natom2name[Z] + get_plus(Z - N)
        assert atom==spec, "Identification of atoms in photoionisation not working"
        maxdN=N
        for dN in range(N,0,-1):
            atom_ionized = natom2name[Z] + get_plus(Z - N + dN)
            if not sp2idx(atom_ionized) in species:
                maxdN = dN-1
                continue
            name_file = "%s__%s" % (atom, atom_ionized) + "_E"*dN
            print("Calculate Verner cross sections for ", name_file)

            erange_eV, xsecs = calc_crossSection_Verner(Z, Q, dN=dN, maxdN=maxdN)
            erange = erange_eV * ev2erg
            data_all[name_file] = {"energy": erange,
                                   "photoionisation": xsecs,
                                   "file": name_file + ".dat"}

    energy = find_energy(data_all, photo_limits)
    save_xsecs(energy, data_all)
    if "BB@" in radiation_type:
        print("field: {} K blackbody".format(tbb))
        tbb = float(radiation_type.replace("BB@", ""))
        prepare_bb_radiation(energy, tbb=tbb)
    elif radiation_type.lower() == "draine":
        print("field: draine")
        prepare_draine_radiation(energy)
    elif radiation_type.lower() == "xdr":
        print("field: {} xdr")
        prepare_xdr(energy)
    elif ".dat" in radiation_type:
        print("field: {}".format(radiation_type))
        spectrum = radiation_type.split("@")[0]
        luminosity = radiation_type.split("@")[1]
        band_lo, band_hi = radiation_type.split("@")[2].split("-")
        addBB = not "noBB" in radiation_type
        if "@Lacc" in radiation_type:
            Lacc = float(radiation_type.split("@")[3][4:])
        else:
            Lacc = None
        if "_photon" in spectrum:
            #spectrum=spectrum.split('_')[0]
            fluxUnits='photon'
        else:
            fluxUnits='energy'
        prepare_external_spec(energy, spectrum, L_X=float(luminosity), X_lo=float(band_lo), X_hi=float(band_hi), add_BB=addBB, Lacc=Lacc, fluxUnits=fluxUnits)
    else:
        sys.exit("ERROR: unknown radiation type %s" % radiation_type)
    compute_habing_flux(energy)

    return energy


def find_energy(data_all, photo_limits):
    #emin = 1e99
    #emax = -1e99
    #for file_name, data in data_all.items():
    #    rr, pp = file_name.split("__")
    #    if rr.count("+") != pp.count("+"):
    #        what = "photoionisation"
    #    else:
    #        what = "photodissociation"
    #    xsecs = data[what]
    #    emin = min(emin, data["energy"][xsecs > 1e-40].min())
    #    emax = max(emax, data["energy"][xsecs > 1e-40].max())
    #    print(data["file"], emin * erg2ev, emax * erg2ev)

    emin = energy_min  # min(emin, energy_min)
    emax = energy_max  # min(emax, energy_max)
    print("radiation field range, eV", emin * erg2ev, emax * erg2ev)

    denergy = 1e-3 * ev2erg
    photo_limits = np.concatenate([photo_limits, photo_limits+denergy, photo_limits-denergy])
    photo_limits = np.array([x for x in photo_limits if emin < x < emax])
    
    nphoto_left = nphoto - len(photo_limits) - (fuv_energy1>emin) - (fuv_energy2>emin)
    if photo_logspacing:
        energy = np.logspace(np.log10(emin), np.log10(emax), nphoto_left)
    else:
        energy = np.linspace(emin, emax, nphoto_left)

    # Add photo_limits into energy list
    energy = np.sort(np.concatenate([energy, photo_limits]))
    # Add fuv energies into energy list
    if fuv_energy1>emin:
        energy = np.sort(np.concatenate([energy, [fuv_energy1]]))
    if fuv_energy2>emin:
        energy = np.sort(np.concatenate([energy, [fuv_energy2]]))

    if len(energy) != nphoto:
        sys.exit("ERROR: number of created energy points is different from expected!")

    np.savetxt(data_dir+"energy.dat", energy)
    return energy


def save_xsecs(energy, data_all):
    for file_name, data in data_all.items():
        rr, pp = file_name.split("__")
        if rr.count("+") != pp.count("+"):
            what = "photoionisation"
        else:
            what = "photodissociation"
        f = interp1d(np.log10(data["energy"]), np.log10(data[what] + 1e-80), bounds_error=False, fill_value=-80.)
        xsecs = 1e1**f(np.log10(energy))
        np.savetxt(data_dir+"photo_xsecs_%s.dat" % file_name, xsecs)


def fplanck(energy, tbb):
    nu = energy / hplanck  # erg -> 1/s
    exp = np.minimum(hplanck * nu / kboltzmann / tbb, 1e2)
    return 2e0 * hplanck * nu ** 3 / clight ** 2 / (np.exp(exp) - 1e0)


def prepare_bb_radiation(energy, tbb=7775.):
    bfield = fplanck(energy, tbb)
    np.savetxt(data_dir+"radiation_field.dat", bfield)


def get_draine_field(energy):
    energy_ev = energy * erg2ev
    bfield = (1.658e6 * energy_ev - 2.152e5 * energy_ev**2 + 6.919e3 * energy_ev**3) * hplanck_eV * energy_ev / erg2ev
    return bfield * ((energy_ev >= 6e0) & (energy_ev <= 13.6)).astype(float)


def prepare_draine_radiation(energy):
    bfield = get_draine_field(energy)
    np.savetxt(data_dir+"radiation_field.dat", bfield)


def prepare_xdr(energy):
    # read: eV, eV/cm2/Hz/s
    energy_xdr, flux_xdr = np.loadtxt("../data/spectra/xdr_spectrum.dat").T
    fxdr = interp1d(np.log10(energy_xdr * ev2erg), np.log10(flux_xdr * ev2erg), fill_value=-99., bounds_error=False)
    np.savetxt(data_dir+"radiation_field.dat", 1e1**fxdr(np.log10(energy)))


def prepare_external_spec(energy, spectrum, L_X=1e30, X_lo=1e2, X_hi=1e4, rstar=6.957e10, add_BB=False, Lacc=None, fluxUnits='energy'):
    if fluxUnits=='energy':
        ## Assume MOCASSIN-style spectrum where columns are lambda and F_lambda
        # Import
        lambdaA, F_lambda = np.genfromtxt("../data/spectra/"+spectrum, usecols=(0,1), unpack=True, comments='!')
        # Convert lambda -> nu
        F_nu = F_lambda * lambdaA**2
        F_nu = F_nu[::-1]
        nuHz = clight/(lambdaA*1e-8)    # Hz
        nuHz = nuHz[::-1]
        # Do all the normalisation
        if 'WG17' in spectrum:
            bfield = np.zeros_like(energy)
            for eV_nonzero, F_nonzero in zip((nuHz * hplanck_eV)[F_nu>0], F_nu[F_nu>0]):
                i_nearest = np.argmin(np.abs(np.log10(energy * erg2ev) - np.log10(eV_nonzero)))
                bfield[i_nearest] = F_nonzero
            #F_interp = interp1d(np.log10(nuHz * hplanck_eV), F_nu, bounds_error=False, fill_value=0.0)
            #bfield = F_interp(np.log10(energy * erg2ev))#/(4*np.pi**2*rstar**2)
            print(F_nu[F_nu>0])
            print(bfield[bfield>0])
            print((energy * erg2ev)[bfield>0])
        else:
            F_interp = interp1d(np.log10(nuHz * hplanck_eV), np.log10(F_nu), bounds_error=False, fill_value=-np.inf)
            bfield = 1e1**F_interp(np.log10(energy * erg2ev))
        X_band = (energy * erg2ev > X_lo) * (energy * erg2ev < X_hi)        
        L_band = integrate.trapz(bfield[X_band], energy[X_band]/hplanck)
        Multiplier = L_X/L_band
        bfield *= Multiplier
        bfield /= (4*np.pi**2*rstar**2)
    elif fluxUnits=='photon':
        ## Assume Nakatani-style spectrum where columns are photon energy, energy bin width, and photon flux
        # Import
        EkeV, F_photon = np.genfromtxt("../data/spectra/"+spectrum, usecols=(0,2), unpack=True, comments='!')
        # Convert to energy flux
        F_nu = EkeV*F_photon
        # Do all the normalisation
        F_interp = interp1d(np.log10(EkeV*1000), np.log10(F_nu), bounds_error=False, fill_value="extrapolate")
        bfield = 1e1**F_interp(np.log10(energy * erg2ev))
        X_band = (energy * erg2ev > X_lo) * (energy * erg2ev < X_hi)        
        L_band = integrate.trapz(bfield[X_band], energy[X_band]/hplanck)
        Multiplier = L_X/L_band
        bfield *= Multiplier
        bfield /= (4*np.pi**2*rstar**2)        
    else:
        raise NotImplementedError("Spectrum units not recognised: should either be in erg/s/A(/cm^2/sr) or photon/s/keV(/cm^2)")

    # Add stellar blackbody
    if add_BB:
        bfield += fplanck(energy, 5000.)
    # Add accretion blackbody
    if Lacc:
        TFUV = 12000
        ffill = Lacc * (TFUV/5780)**-4 * (rstar/6.957e10)**-2
        bfield += ffill * fplanck(energy, TFUV)

    np.savetxt(data_dir+"radiation_field.dat", bfield)


def compute_habing_flux(energy):
    field = get_draine_field(energy) * 4. * np.pi
    cond = (energy >= fuv_energy1) & (energy <= fuv_energy2)
    hf = np.trapz(field[cond] / energy[cond], energy[cond])

    idx1 = np.argmin(np.abs(energy - fuv_energy1))+1
    idx2 = np.argmin(np.abs(energy - fuv_energy2))+1

    rconst = "! FUV interval photobins indexes\n"
    rconst += "integer,parameter::fuv_idx1=%d\n" % idx1
    rconst += "integer,parameter::fuv_idx2=%d\n" % idx2
    rconst += "real*8,parameter::habing_flux=%s\n" % ("%.18e" % hf).replace("e", "d")

    preprocess("../prizmo_commons.f90", {"RADIATION_CONSTANTS": rconst})
    print("Habing flux:", hf, idx1, idx2)
