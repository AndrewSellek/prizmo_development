import os.path
import numpy as np
import matplotlib.pyplot as plt
from prizmo_commons import clight, hplanck, erg2ev, print_title, fuv_energy1, fuv_energy2, plotOn, data_dir, NCPAH
from prizmo_preprocess import preprocess

def Drude_profile(wav, wav_j, gamma_j, sigma_int_j):
    # wav_j is in micron -> need wav in micron
    return 2/np.pi * gamma_j*wav_j*sigma_int_j / (np.pi * (wav/wav_j-wav_j/wav)**2 + gamma_j**2) * 1e-24

def HperC(NC):
    if NC<=25:
        return 0.5
    elif NC>=100:
        return 0.25
    else:
        return 0.5*(25/NC)**0.5

def cutoff(wav, NC, charged=False):
    # Number of benzoid rings
    if NC>=40:
        M = 0.4*NC
    else:
        M = 0.3*NC
    # Cut-off wavelength
    if charged:
        wav_c = 1/(2.282/np.sqrt(M)+0.889)
    else:
        wav_c = 1/(3.804/np.sqrt(M)+1.052)
    # Desert+1990 cutoff function
    y = wav_c/wav
    return 1/np.pi*np.arctan(1e3*(y-1)**3/y)+0.5

def prepare(user_energy, charged=False, NC=100):
    
    print_title("PAH opacity")
    fname = data_dir+"sigma_PAH.dat"

    if charged:
        resonanceData = np.genfromtxt('../data/PAH_opacity/drude_parameters_DL07.dat', usecols=(0,1,2,4), names=True, dtype=(int,float,float,'U9'))
        resonanceData.dtype.names = resonanceData.dtype.names[:-1]+('sigma_int_j',)
    else:
        resonanceData = np.genfromtxt('../data/PAH_opacity/drude_parameters_DL07.dat', usecols=(0,1,2,3), names=True, dtype=(int,float,float,'U9'))

    wav_um = hplanck*clight/user_energy * 1e4
    sigma_PAH = np.zeros_like(user_energy)

    # wav_um<1/17.25 -> EUV/X-ray: behaves like graphite
    sigma_PAH += np.zeros_like(wav_um) * (1/wav_um>17.25)

    # 1/17.25<wav_um<1/15 -> EUV: continuum
    sigma_PAH += (126.0-6.4943/wav_um)*1e-18 * (1/wav_um>15) * (1/wav_um<17.25)
    
    # 1/15<wav_um<1/10 -> FUV/EUV: continuum + j=1 feature
    sigma_PAH += (-3.0+1.35/wav_um)*1e-18 * (1/wav_um>10) * (1/wav_um<15)
    sigma_PAH += Drude_profile(wav_um, resonanceData['lambda_j'][0], resonanceData['gamma_j'][0], float(resonanceData['sigma_int_j'][0])) * (1/wav_um>10) * (1/wav_um<15)

    # 1/10<wav_um<1/7.7 -> FUV: continuum
    sigma_PAH += (66.302-24.367/wav_um+2.950/wav_um**2-0.1057/wav_um**3)*1e-18 * (1/wav_um>7.7) * (1/wav_um<10)
    
    # 1/7.7<wav_um<1/5.9 -> FUV: continuum + j=2 feature
    sigma_PAH += (1.8687+0.1905/wav_um+0.4175*(1/wav_um-5.9)**2+0.04370*(1/wav_um-5.9)**3)*1e-18 * (1/wav_um>5.9) * (1/wav_um<7.7)
    sigma_PAH += Drude_profile(wav_um, resonanceData['lambda_j'][1], resonanceData['gamma_j'][1], float(resonanceData['sigma_int_j'][1])) * (1/wav_um>5.9) * (1/wav_um<7.7)
        
    # 1/5.9<wav_um<1/3.3 -> FUV: continuum + j=2 feature
    sigma_PAH += (1.8687+0.1905/wav_um)*1e-18 * (1/wav_um>3.3) * (1/wav_um<5.9)
    sigma_PAH += Drude_profile(wav_um, resonanceData['lambda_j'][1], resonanceData['gamma_j'][1], float(resonanceData['sigma_int_j'][1])) * (1/wav_um>3.3) * (1/wav_um<5.9)
        
    # wav_um>1/3.3 -> NUV/optical/IR features: continuum + j>=3 features
    sigma_PAH += 34.58*np.power(10,-18-3.431*wav_um) * cutoff(wav_um, NC, charged) * (1/wav_um<3.3)
    for j in range(3,len(resonanceData['j'])):
        Hband = "(H/C)" in resonanceData['sigma_int_j'][j]
        resonanceData['sigma_int_j'][j] = resonanceData['sigma_int_j'][j].replace("(H/C)","")
        if Hband:
            sigma_PAH += Drude_profile(wav_um, resonanceData['lambda_j'][j], resonanceData['gamma_j'][j], float(resonanceData['sigma_int_j'][j])) * (1/wav_um<3.3) * HperC(NC)
        else:
            sigma_PAH += Drude_profile(wav_um, resonanceData['lambda_j'][j], resonanceData['gamma_j'][j], float(resonanceData['sigma_int_j'][j])) * (1/wav_um<3.3)

    if plotOn:
        plt.loglog(wav_um, sigma_PAH)
        plt.xlabel('$\lambda\,/\,\mathrm{\mu m}$')
        plt.ylabel('$\sigma\,/\,\mathrm{cm^{2}/C}$')
        plt.show()
        
    # Scale to number of C atoms
    #sigma_PAH*=NC

    # Save
    np.savetxt(fname, sigma_PAH)

