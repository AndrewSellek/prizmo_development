import numpy as np

natom2name = {1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 10: "Ne", 11: "Na", 12: "Mg",
              13: "Al", 14: "Si", 15: "P", 16: "S", 17: "Cl", 18: "Ar", 19: "K", 20: "Ca", 21: "Sc", 22: "Ti", 23: "V",
              24: "Cr", 25: "Mn", 26: "Fe", 27: "Co", 28: "Ni", 29: "Cu", 30: "Zn", 31: "Ga", 32: "Ge", 33: "As",
              34: "Se", 35: "Br", 36: "Kr", 37: "Rb", 38: "Sr", 39: "Y", 40: "Zr", 41: "Nb", 42: "Mo", 43: "Tc",
              44: "Ru", 45: "Rh", 46: "Pd", 47: "Ag", 48: "Cd", 49: "In", 50: "Sn", 51: "Sb", 52: "Te", 53: "I",
              54: "Xe", 55: "Cs", 56: "Ba", 57: "La", 58: "Ce", 59: "Pr", 60: "Nd", 61: "Pm", 62: "Sm", 63: "Eu",
              64: "Gd", 65: "Tb", 66: "Dy", 67: "Ho", 68: "Er", 69: "Tm", 70: "Yb", 71: "Lu", 72: "Hf", 73: "Ta",
              74: "W", 75: "Re", 76: "Os", 77: "Ir", 78: "Pt", 79: "Au", 80: "Hg", 81: "Tl", 82: "Pb", 83: "Bi",
              84: "Po", 85: "At", 86: "Rn", 87: "Fr", 88: "Ra", 89: "Ac", 90: "Th", 91: "Pa", 92: "U", 93: "Np"}

def get_plus(val):
  if val > 200:
      plus = "+" + str(val)
  else:
      plus = "+"*val
  return plus



for row in open("cfit.dat"):
  srow = row.strip()

  # Tmin/eV, Tmax/keV
  Z, N, dE, P, A, X, K, Tmin, Tmax = [float(x) for x in srow[6:].split()]
  Z = int(Z)
  N = int(N)

  atom = natom2name[Z] + get_plus(Z - N)
  atom_ionized = natom2name[Z] + get_plus(Z - N + 1)

  verbatim = "%s + E -> %s + E + E" % (atom,atom_ionized)

  line = verbatim.ljust(40) + " [] " + "%e * (1e0 + %e*invsqrTe) / (%e + %e * invTe) * %e * invTe**(%.2f) * exp(-%.2f * invTe)" % (A, P * np.sqrt(dE), X, dE, dE**K, K, dE)

  line = line.replace(" + 0.000000e+00*invsqrTe", "").replace("(1e0)", "1e0")

  print(line)
