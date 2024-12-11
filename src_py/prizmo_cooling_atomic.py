import os.path
import numpy as np
import matplotlib.pyplot as plt
# import ChiantiPy.core as ch
from tqdm import tqdm
from prizmo_commons import idx2sp, kboltzmann, clight, hplanck, print_title, plotOn, sp2spj, data_dir
from prizmo_preprocess import preprocess
from scipy.interpolate import interp1d
import sys
import shutil

def prepare_line_functions(atom, nlevels):
    lineFunc = lineFunc_c = "! *****************\n"

    lineFunc += "function prizmo_get_atomic_lines_%s(x, Tgas) result(lines)\n" % sp2spj(atom)
    lineFunc += "  use prizmo_lines_atomic\n"
    lineFunc += "  implicit none\n"
    lineFunc += "  real*8,intent(in)::x(nspecies), Tgas\n"
    lineFunc += "  real*8::lines(%d)\n\n" % int(nlevels*(nlevels-1)/2)
    lineFunc += "  lines = atomic_lines_%s(x, log10(Tgas))\n\n" % sp2spj(atom)
    lineFunc += "end function prizmo_get_atomic_lines_%s\n\n" % sp2spj(atom)

    lineFunc_c += "subroutine prizmo_get_atomic_lines_%s_c(x, Tgas, lines) bind(C)\n" % sp2spj(atom)
    lineFunc_c += "  use prizmo_commons\n"
    lineFunc_c += "  real(C_DOUBLE),intent(in)::x(nspecies), Tgas\n"
    lineFunc_c += "  real(C_DOUBLE),intent(out)::lines(%d)\n\n" % int(nlevels*(nlevels-1)/2)
    lineFunc_c += "  lines = prizmo_get_atomic_lines_%s(x, Tgas)\n\n" % sp2spj(atom)
    lineFunc_c += "end subroutine prizmo_get_atomic_lines_%s_c\n\n" % sp2spj(atom)

    return lineFunc, lineFunc_c


def prepare_atomic_cooling(species_indexes, H2_inc, fname="../data/atomic_cooling/krome_data.dat"):
    print_title("atomic cooling")
    prepare_atomic_cooling_levels(H2_inc, fname=fname)


def prepare_atomic_cooling_levels(H2_inc, fname="../data/atomic_cooling/krome_data.dat"):
    atoms = ["C", "O", "C+", "O+", "Ne+", "Ne++", "Ar+", "Ar++"]

    funcs = lines = loaders = commons = cool_tot = lineFuncs = lineFuncs_c = ""
    cool_arr = "cools(1) = atomic_cooling_H(x, 1d1**log_Tgas)\n"
    icount = 0
    for atom in atoms:
        print(atom)
        data, _ = krome_cooling(atom, fname=fname)
        if data["nlevels"] in [2, 3, 5]:
            fs, As, ls, cs, ct = prepare_xlevel(data, atom, data["nlevels"], H2_inc)
            cool_arr += "cools(%d) = atomic_cooling_%s(x, log_Tgas)\n" % (icount + 2, sp2spj(atom))
            lFs, lFcs = prepare_line_functions(atom, data["nlevels"])
            icount += 1
        else:
            continue
        funcs += fs
        lines += As
        loaders += ls
        commons += cs
        cool_tot += ct
        lineFuncs += lFs
        lineFuncs_c += lFcs

    commons += "real*8,parameter::natomic_cools=%d\n" % (icount + 1)


    preprocess("prizmo_loaders.f90", {"LOAD_ATOMIC_COOLING": loaders})
    preprocess("prizmo_commons.f90", {"ATOMIC_COOLING_COMMONS": commons})
    preprocess("prizmo_cooling_atomic.f90", {"ATOMIC_COOLING_FUNCTIONS": funcs,
                                                "ATOMIC_COOLING": cool_tot,
                                                "ATOMIC_COOLING_ARRAY": cool_arr})
    preprocess("prizmo_lines_atomic.f90", {"ATOMIC_LINES_FUNCTIONS": lines})
    preprocess("prizmo.f90", {"GET_ATOMIC_LINES": lineFuncs})
    preprocess("prizmo_c.f90", {"GET_ATOMIC_LINES_C": lineFuncs_c})


def prepare_atomic_cooling_tables(species_indexes):

    nt = 100  # number of temperature grid points
    ne = 30  # electrons grid points
    ni = 30  # proton grid points
    nn = 30  # H grid points

    # avail_coolants = ["C", "O", "C+", "O+", "H", "O++", "O+++", "O++++",
    #                  "C++", "C+++", "C++++"]
    # method = ["KROME", "KROME", "KROME", "KROME", "CHIANTI", "CHIANTI", "CHIANTI", "CHIANTI",
    #          "CHIANTI", "CHIANTI", "CHIANTI"]
    # colliders = [["H+", "H", "e"], ["H+", "H", "e"], ["H", "e"], ["e"], ["e"], ["e"], ["e"], ["e"],
    #             ["e"], ["e"], ["e"]]

    avail_coolants = ["H", "O++", "O+++", "O++++", "C++", "C+++", "C++++"]
    method = ["CHIANTI", "CHIANTI", "CHIANTI", "CHIANTI", "CHIANTI", "CHIANTI", "CHIANTI"]
    colliders = [["e"], ["e"], ["e"], ["e"], ["e"], ["e"], ["e"], ["e"]]

    t = np.logspace(0, 6, nt)

    atoms = []
    atomic_cooling2 = atomic_cooling3 = atomic_cooling4 = ""
    atomic_cooling_loader2 = []
    atomic_cooling_loader3 = []
    atomic_cooling_loader4 = []
    icount2 = icount3 = icount4 = 0
    for idx in species_indexes:
        a = idx2sp(idx)
        if a not in avail_coolants:
            continue
        atoms.append(a)
        ii = avail_coolants.index(a)
        ncolliders = len(colliders[ii])
        if ncolliders == 1:
            icount2 += 1
            atomic_cooling2 += "cool = cool + coola2d(%d) * x(%s)\n" % (icount2, idx)
            atomic_cooling_loader2.append("\"%s\"" % ("runtime_data/cool_%s.dat" % a).ljust(30))
        elif ncolliders == 2:
            icount3 += 1
            atomic_cooling3 += "cool = cool + coola3d(%d) * x(%s)\n" % (icount3, idx)
            atomic_cooling_loader3.append("\"%s\"" % ("runtime_data/cool_%s.dat" % a).ljust(30))
        elif ncolliders == 3:
            icount4 += 1
            atomic_cooling4 += "cool = cool + coola4d(%d) * x(%s)\n" % (icount4, idx)
            atomic_cooling_loader4.append("\"%s\"" % ("runtime_data/cool_%s.dat" % a).ljust(30))
        else:
            sys.exit("ERRORS colliders number %d is unknown, it should be 1-3" % ncolliders)

    atomic_cooling = atomic_cooling2 + atomic_cooling3 + atomic_cooling4

    atomic_cooling_loader_str = "fnames2d = (/" + ", &\n".join(atomic_cooling_loader2) + "/)\n\n"
    atomic_cooling_loader_str += "fnames3d = (/" + ", &\n".join(atomic_cooling_loader3) + "/)\n\n"
    atomic_cooling_loader_str += "fnames4d = (/" + ", &\n".join(atomic_cooling_loader4) + "/)\n\n"

    number_of_atomic_coolants = ""
    for d in range(2, 5):
        nv = [len(x) for x in colliders].count(d-1)
        number_of_atomic_coolants += "integer,parameter::atomic_cooling_nvec%dd=%d\n" % (d, nv)

    preprocess("prizmo_commons.f90", {"ATOMIC_COOLING_NVEC": number_of_atomic_coolants})
    preprocess("prizmo_cooling_atomic.f90", {"ATOMIC_COOLING": atomic_cooling})
    preprocess("prizmo_loaders.f90", {"LOAD_ATOMIC_COOLING": atomic_cooling_loader_str})

    colors = ["r", "g", "b", "tab:orange"]
    lss = ["-", ":", "--", "-."]
    for a in atoms:
        ii = avail_coolants.index(a)
        plt.clf()
        fname = data_dir+"cool_%s.dat" % a
        fname_database = "../data/atomic_cooling/cool_%s.dat" % a
        if os.path.isfile(fname):
            print("skipping, cooling file found", fname)
            continue
        if os.path.isfile(fname_database):
            print("skipping, cooling file found in data", fname_database)
            shutil.copyfile(fname_database, fname)
            continue
        print("Atom " + a)
        icount = 0
        if colliders[ii] == ["H+", "H", "e"]:
            tab = np.zeros((5, nt * ne * ni * nn))
            for k, dn in enumerate(tqdm(np.logspace(-6, 20, nn))):
                for j, di in enumerate(np.logspace(-6, 20, ni)):
                    for i, de in enumerate(np.logspace(-6, 20, ne)):
                        em = emission(method, avail_coolants, a, dn, di, de, t)
                        tab[0, icount:icount+nt] = dn
                        tab[1, icount:icount+nt] = di
                        tab[2, icount:icount+nt] = de
                        tab[3, icount:icount+nt] = t
                        tab[4, icount:icount+nt] = em
                        # plt.loglog(t, tab[3, icount:icount+nt], color=colors[ii], label="H+=%.1e e=%.1e" % (di, d), ls=lss[ii])
                        icount += nt
                        # plt.loglog(t, 7.3e-19*d*np.exp(-118400/t), color="k", ls="--")
        elif colliders[ii] == ["H", "e"]:
            tab = np.zeros((4, nt * ne * nn))
            for k, dn in enumerate(tqdm(np.logspace(-6, 20, nn))):
                for i, de in enumerate(np.logspace(-6, 20, ne)):
                    di = 0e0
                    em = emission(method, avail_coolants, a, dn, di, de, t)
                    tab[0, icount:icount+nt] = dn
                    tab[1, icount:icount+nt] = de
                    tab[2, icount:icount+nt] = t
                    tab[3, icount:icount+nt] = em
                    icount += nt
        elif colliders[ii] == ["e"]:
            tab = np.zeros((3, nt * ne))
            for i, de in enumerate(tqdm(np.logspace(-6, 20, ne))):
                di = dn = 0e0
                em = emission(method, avail_coolants, a, dn, di, de, t)
                tab[0, icount:icount+nt] = de
                tab[1, icount:icount+nt] = t
                tab[2, icount:icount+nt] = em
                icount += nt
        else:
            sys.exit("ERROR: unknown colliders combination " + ", ".join(colliders[ii]))

        np.savetxt(fname, tab.T)

        if plotOn:
            plt.ylim(bottom=1e-30)
            plt.show()
        else:
            plt.close()


def emission(method, avail_coolants, a, dn, di, d, t):
    if method[avail_coolants.index(a)] == "CHIANTI":
        chianti_atom = (a.replace("+", "") + "_" + str(a.count("+") + 1)).lower()
        h1 = ch.ion(chianti_atom, temperature=t, eDensity=d, abundance=1e0, pDensity=np.ones_like(t) * di)
        h1.emiss()
        em = h1.Emiss["emiss"] * 4e0 * np.pi  # CHIANTI is in erg/s/cm3/sr
        em = np.clip(np.sum(em, axis=0), 1e-40, 1e99)
    elif method[avail_coolants.index(a)] == "KROME":
        data = krome_cooling(a)
        em = np.array([cool(data, ["H", "H+", "e"], [dn, di, d], x) for x in t])
    elif method[avail_coolants.index(a)] == "LAMDA":
        chianti_atom = (a.replace("+", "") + "_" + str(a.count("+") + 1)).lower()
        data = lamda_cooling("../data/lamda/%s.dat" % chianti_atom)
        em = np.array([cool(data, ["H", "H+", "e"], [dn, di, d], x) for x in t])
    else:
        sys.exit("ERROR: unknown method " + method[avail_coolants.index(a)])

    return em


def rate2fit(expr, species, collider, gu, strength=False):
    expr = expr.lower()
    expr = expr.replace("d", "e")
    expr = expr.replace("exp(", "np.exp(")
    expr = expr.replace("expinvt2", "np.exp(1e2/tgas)")
    expr = expr.replace("sqrt(", "np.sqrt(")
    expr = expr.replace("log10(", "np.log10(")
    expr = expr.replace("log(", "np.log(")
    expr = expr.replace("lnt4", "(np.log(tgas*1e-4))")
    expr = expr.replace("lnt", "(np.log(tgas))")
    expr = expr.replace("t2", "(tgas*1e-2)")
    expr = expr.replace("t4", "(tgas*1e-4)")
    expr = expr.replace("invtgas", "(1e0/tgas)")
    expr = expr.replace("; ", "\n")
    expr = "f[i] = " + expr
    expr = expr.replace("):", "): f[i] = ")

    def w(x, x0, s):
        return (np.tanh(s * (x - x0)) + 1.) / 2.

    trange = np.logspace(0, 6, 10000)
    f = np.zeros_like(trange)
    for i, tgas in enumerate(trange):
        exec(expr)
    if strength:
        f*=8.629e-8/gu/np.sqrt(trange/1.0e4)

    if plotOn:
        plt.title("%s + %s" % (species, collider))
        plt.loglog(trange, f)
        plt.show()
    else:
        plt.close()

    return trange, f


# ************************
def prepare_xlevel(data, atom, nlevels, H2_inc, nt=10000, multiplet_hierachy=False):

    from prizmo_commons import sp2spj, sp2idx, py2f90
    multipletE, uniqueE = np.unique(data['multipletE'], return_index=True)
    multiplets = np.array(data['multiplets'])[uniqueE]
    nmultiplets = len(multiplets)

    aa = [[[] for _ in range(nlevels)] for _ in range(nlevels)]

    for i in range(nlevels):
        for j in range(nlevels):
            if j > i:
                aa[i][j] = [py2f90(data["Aul"][j, i])]
            if j == i:
                aa[i][j] = [py2f90(data["Aul"][j, k]) for k in range(j)]

    defs = []
    fits = loader = commons = cool_tot = ""

    has_ortho_para = False

    for collider in data["rates"]:
        if (collider == 'H2or' or collider =='H2pa') and not H2_inc:
            continue
        nlev_coll = data["collider_max"][collider]+1
        assert nlev_coll<=nlevels
        print(atom, collider, nlev_coll)

        if "or" in collider or "pa" in collider:
            has_ortho_para = True

        spj = sp2spj(collider)
        idx = sp2idx(collider)

        fhk = [[None for _ in range(nlev_coll)] for _ in range(nlev_coll)]
        vnames = []
        for i in range(1, nlev_coll):
            for j in range(nlev_coll):
                if j != i:
                    kname = "k%d%d" % (j, i)
                    # fhk[j][i] = open(data_dir+"cool_%s_%s_%s.dat" % (atom, collider, kname), "w")
                else:
                    kname = "".join(["k"+str(i)+str(k) for k in range(nlev_coll) if k != i])
                vnames.append(kname)

                ffname = data_dir+"cool_%s_%s_%s.dat" % (atom, collider, kname)
                if os.path.isfile(ffname):
                    print("skipping, cooling file found", ffname)
                    continue
                fhk[j][i] = open(ffname, "w")
                
        if collider in data["rates_m"] and multiplet_hierachy:
            fhk_m = [[None for _ in range(nmultiplets)] for _ in range(nmultiplets)]
            vnames_m = []
            for i, mi in enumerate(multiplets[1:]):
                for j, mj in enumerate(multiplets):
                    if mj != mi:
                        kname = "k%s%s" % (mj, mi)
                        # fhk[j][i] = open(data_dir+"cool_%s_%s_%s.dat" % (atom, collider, kname), "w")
                    else:
                        kname = "".join(["k"+mi+mk for mk in multiplets if mk != mi])
                    vnames_m.append(kname)

                    ffname = data_dir+"cool_%s_%s_%s.dat" % (atom, collider, kname)
                    if os.path.isfile(ffname):
                        print("skipping, cooling file found", ffname)
                        continue
                    print("Opening file {},{} for edits:".format(j,i+1), ffname)
                    fhk_m[j][i+1] = open(ffname, "w")

        fnames = ["\"runtime_data/cool_%s_%s_%s.dat\"" % (atom, collider, x) for x in vnames]

        loader += "fnames_%dlev = (/%s/)\n\n" % (nlev_coll, ", &\n".join(fnames))

        loader += "call load_1d_fit_vec(fnames_%dlev, atomic_cooling_%dlev_nvec, atomic_cooling_n1, &\n" \
                  "atomic_cooling_table_%s_%s, do_log=.false.)\n\n" % (nlev_coll, nlev_coll, sp2spj(atom), spj)

        if collider in data["rates_m"] and multiplet_hierachy:
            fnames_m = ["\"runtime_data/cool_%s_%s_%s.dat\"" % (atom, collider, x) for x in vnames_m]

            loader += "fnames_%dlev = (/%s/)\n\n" % (nmultiplets, ", &\n".join(fnames_m))

            loader += "call load_1d_fit_vec(fnames_%dlev, atomic_cooling_%dlev_nvec, atomic_cooling_n1, &\n" \
                      "atomic_cooling_table_%s_%s, do_log=.false.)\n\n" % (nmultiplets, nmultiplets, sp2spj(atom), spj)

        commons += "type(fit1d_data_vec(nv=atomic_cooling_%dlev_nvec, n1=atomic_cooling_n1))::atomic_cooling_table_%s_%s\n" \
                   % (nlev_coll, sp2spj(atom), spj)

        rates = data["rates"][collider]
        if collider in data["rates_m"] and multiplet_hierachy:
            rates_m = data["rates_m"][collider]
        for tgas in np.linspace(0, 6, nt):
            for i in range(1, nlev_coll):
                for j in range(nlev_coll):
                    if fhk[j][i] is None:
                        continue
                    if j != i:
                        kk = rates[j, i](tgas)
                    else:
                        kk = np.log10(sum([1e1**rates[i, k](tgas) for k in range(nlev_coll) if i != k]))
                    fhk[j][i].write("%17.8e %17.8e\n" % (tgas, kk))
            if collider in data["rates_m"] and multiplet_hierachy:
                for i, mi in enumerate(multiplets[1:]):
                    for j, mj in enumerate(multiplets):
                        if fhk_m[j][i+1] is None:
                            continue
                        if mj != mi:
                            kk = rates_m[j, i+1](tgas)
                        else:
                            kk = np.log10(sum([1e1**rates_m[i+1, k](tgas) for k in range(nmultiplets) if i+1 != k]))
                        fhk_m[j][i+1].write("%17.8e %17.8e\n" % (tgas, kk))

        for i in range(1, nlev_coll):
            for j in range(nlev_coll):
                if fhk[j][i] is None:
                    continue
                fhk[j][i].close()
                
        if collider in data["rates_m"] and multiplet_hierachy:
            for i in range(1,nmultiplets):
                for j in range(nmultiplets):
                    if fhk_m[j][i] is None:
                        continue
                    fhk_m[j][i].close()

        defs.append("kfit_%s(atomic_cooling_%dlev_nvec)" % (spj, nlev_coll))

        fits += "  kfit_%s = 1d1**interp_1dfit_vec(log_Tgas, atomic_cooling_table_%s_%s, & \n" % (spj, sp2spj(atom), spj)
        fits += "    atomic_cooling_%dlev_nvec, atomic_cooling_n1)\n\n" % nlev_coll

        count = 0
        for i in range(1, nlev_coll):
            for j in range(nlev_coll):
                kfit_name = "kfit_%s(%d) * x(%s)" % (spj, count + 1, idx)
                aa[i][j].append(kfit_name)
                count += 1

    fun = "! *****************\n"
    fun += "function atomic_cooling_%s(x, log_Tgas) result(cool)\n" % sp2spj(atom)
    lin = "! *****************\n"
    lin += "function atomic_lines_%s(x, log_Tgas) result(cool)\n" % sp2spj(atom)
    head = "  use prizmo_commons\n"
    head += "  use prizmo_linear_solver\n"
    head += "  implicit none\n"
    head += "  real*8,intent(in)::x(nspecies), log_Tgas\n"
    head += "  real*8::" + ", ".join(np.unique(defs)) + "\n"
    head += "  real*8::b(%d), A(%d, %d), n(%d), H2or, H2pa\n" % (nlevels, nlevels, nlevels, nlevels)
    fun += head
    lin += head
    fun += "  real*8::cool\n\n"
    lin += "  real*8::lines(atomic_cooling_%dlev_nvec/2)\n\n" % nlevels

    if has_ortho_para:
        fun += "  H2or = x(idx_H2) * ortho_to_para / (ortho_to_para + 1d0)\n"
        fun += "  H2pa = x(idx_H2) / (ortho_to_para + 1d0)\n\n"
        lin += "  H2or = x(idx_H2) * ortho_to_para / (ortho_to_para + 1d0)\n"
        lin += "  H2pa = x(idx_H2) / (ortho_to_para + 1d0)\n\n"

    #fun += "#ifdef NO%sCOOL\n\n" % sp2spj(atom)
    #fun += "  cool = 0d0\n\n"
    #fun += "#else\n\n"

    fun += fits
    lin += fits

    A = ""
    for i in range(1, nlevels):
        for j in range(nlevels):
            if i == j:
                sgn = " - "
            else:
                sgn = " + "
            A += "  A(%d, %d) =%s" % (i+1, j+1, sgn) + sgn.join(aa[i][j]) + "\n"
    A += "\n"

    A = A.replace("x(idx_H2pa)", "H2pa")
    A = A.replace("x(idx_H2or)", "H2or")
    fun += A
    lin += A

    fun += "  b = (/1d0, %s/)\n" % ", ".join(["0d0" for _ in range(nlevels - 1)])
    fun += "  n = linear_solver_n%d(A, b)\n\n" % nlevels
    lin += "  b = (/1d0, %s/)\n" % ", ".join(["0d0" for _ in range(nlevels - 1)])
    lin += "  n = linear_solver_n%d(A, b)\n\n" % nlevels

    fun += "  cool = 0d0\n"
    for i in range(1, nlevels):
        esum = []
        for j in range(i):
            de = data["deltaE"][i] - data["deltaE"][j]
            esum.append("%s * %s" % (py2f90(data["Aul"][i, j]), py2f90(de)))
        fun += "  cool = cool + n(%d) * (%s)\n" % (i+1, " + ".join(esum))
    fun += "  cool = cool * x(%s)\n" % sp2idx(atom)
    fun += "  cool = max(cool, 0d0)\n\n"

    nline=0
    for i in range(1, nlevels):
        for j in range(i):
            nline+=1
            de = data["deltaE"][i] - data["deltaE"][j]
            lin += "  lines(%d) = n(%d) * %s * %s\n" % (nline, i+1, py2f90(data["Aul"][i, j]), py2f90(de))

    fun = fun.replace("x(x(idx_H2pa))", "x(idx_H2) / (ortho_to_para + 1d0)")
    fun = fun.replace("x(x(idx_H2or))", "x(idx_H2) * ortho_to_para / (ortho_to_para + 1d0)")
    lin = lin.replace("x(x(idx_H2pa))", "x(idx_H2) / (ortho_to_para + 1d0)")
    lin = lin.replace("x(x(idx_H2or))", "x(idx_H2) * ortho_to_para / (ortho_to_para + 1d0)")

    #fun += "#endif\n\n"

    fun += "end function atomic_cooling_%s\n\n" % sp2spj(atom)
    lin += "end function atomic_lines_%s\n\n" % sp2spj(atom)

    # cool_tot += "if(x(%s) > xlimit) then\n" % sp2idx(atom)
    cool_tot += " cool = cool + atomic_cooling_%s(x, log_Tgas)\n" % sp2spj(atom)
    # cool_tot += "end if\n"

    return fun, lin, loader, commons, cool_tot


def krome_cooling(species, fname="../data/atomic_cooling/krome_data.dat"):
    in_metal = False
    data = {"nlevels": 0,
            "weights": [],
            "deltaE": [],
            "terms": [],
            "multiplets": [],
            "multipletE": [],
            "rates": dict(),
            "rates_m": dict(),
            "collider_max": dict()}

    def fzero(arg):
        return -99.

    with open(fname) as file:
        rows = file.read()

    rows = rows.replace("\t", " ")
    while "  " in rows:
        rows = rows.replace("  ", " ")
    rows = rows.replace("\n if(", "; if(").split("\n")

    for row in rows:
        srow = row.strip()
        if srow.startswith("#") or srow == "":
            continue
        if srow.replace(" ", "") == "metal:" + species:
            in_metal = True
            continue
        if srow == "endmetal":
            in_metal = False
        if not in_metal:
            continue

        if srow.startswith("level:"):

            if len(srow.split(",")) == 4:
                _, deltaE, g, term = srow.split(",")
                term = term.replace(' ','')
            else:
                raise Exception("Atomic energy level not specified in expected format (n, deltaE, g, term):", srow)
                #_, deltaE, g = srow.split(",")
            data["nlevels"] += 1
            data["weights"].append(float(g))
            data["deltaE"].append(float(deltaE) * kboltzmann)  # erg
            data["terms"].append(term)
            data["multiplets"].append(''.join([term[c] if (term[c].isalpha() or term[:c+1].isnumeric()) else '' for c in range(len(term))]))
        elif "-> " in srow:
            if "Aul" not in data:
                nlevels = data["nlevels"]
                data["Aul"] = np.zeros((nlevels, nlevels))
            up, low, Aul = srow.replace("->", ",").split(",")
            data["Aul"][int(up), int(low)] = float(Aul.replace("d", "e"))

        else:
            nlevels = data["nlevels"]
            collider, up, low, rate = srow.split(",", 3)
            strength = False
            if collider == "E":
                collider = "e"
            elif collider == "Es":
                collider = "e"
                strength = True
            up = int(up)
            low = int(low)

            if collider not in data["rates"]:
                data["rates"][collider] = np.full((nlevels, nlevels), fzero, dtype=object)
                data["collider_max"][collider] = -np.inf
            data["collider_max"][collider] = max(data["collider_max"][collider],up)
            trange, kul = rate2fit(rate, species, collider, data["weights"][up], strength)
            delta = data["deltaE"][up] - data["deltaE"][low]
            klu = kul * data["weights"][up] / data["weights"][low] \
                * np.exp(-delta / kboltzmann / trange)
            klu += 1e-99
            kul += 1e-99
            if kul.min() < 0e0 or klu.min() < 0e0:
                print("ERROR: negative rate coefficient in %s cooling, collider %s!" % (species, collider))
                print(kul.min(), klu.min())
                plt.clf()
                plt.loglog(trange, kul)
                plt.loglog(trange, klu)
                plt.show()
            data["rates"][collider][up, low] = interp1d(np.log10(trange), np.log10(kul))
            data["rates"][collider][low, up] = interp1d(np.log10(trange), np.log10(klu))
            
    multiplets = np.unique(data['multiplets'])
    nmultiplets = len(multiplets)
    if nmultiplets==1:
        return data, nmultiplets
    trange = np.logspace(0, 6, 10000)
    Eeff = {}
    totalweights = {}
    for multiplet in multiplets:
        terms = np.nonzero(np.array(data["multiplets"])==multiplet)[0]
        #print(multiplet, terms)
        deltaE  = np.array(data["deltaE"])[terms]
        weights = np.array(data["weights"])[terms]
        totalweights[multiplet] = np.sum(weights)
        Elim = np.sum(weights*deltaE)/np.sum(weights)
        if len(terms)==1:
            Eeff[multiplet] = np.full_like(trange,Elim)
        else:    
            totalweight = np.sum(weights)
            Eeff[multiplet] = Elim - kboltzmann * trange * np.log(np.sum([g*np.exp(-(E-Elim)/kboltzmann/trange) for g, E in zip(weights, deltaE)], axis=0) / np.sum(weights))
            #plt.loglog(trange, Eeff/kboltzmann)
            #for E in deltaE:
            #    plt.axhline(E/kboltzmann)
            #plt.show()
    for multiplet in data['multiplets']:
        data['multipletE'].append(Eeff[multiplet][-1])
    multiplets = multiplets[np.argsort([E[-1] for E in Eeff.values()])]
    for collider in data["rates"].keys():
        #data["rates_m"][collider] = np.full((nmultiplets, nmultiplets), fzero, dtype=object)
        for l, ml in enumerate(multiplets):
            terms_l = np.nonzero(np.array(data["multiplets"])==ml)[0]
            for u, mu in enumerate(multiplets):
                terms_u = np.nonzero(np.array(data["multiplets"])==mu)[0]
                delta = Eeff[mu] - Eeff[ml]
                if delta[-1]<=0:
                    continue
                kul = np.zeros_like(trange)
                for low in terms_l:
                    for up in terms_u:
                        kpart = data["rates"][collider][up, low](np.log10(trange))
                        if np.max(kpart) > -99:
                            #print("add term {},{} to {},{} for collider {}".format(up,low,mu,ml,collider))
                            kul += 10**kpart
                if np.max(kul)==0:
                    #print("No significant rates for collider", collider)
                    continue
                if collider not in data["rates_m"]:
                    data["rates_m"][collider] = np.full((nmultiplets, nmultiplets), fzero, dtype=object)
                klu = kul * totalweights[mu] / totalweights[ml] \
                    * np.exp(-delta / kboltzmann / trange)
                klu += 1e-99
                kul += 1e-99
                if kul.min() < 0e0 or klu.min() < 0e0:
                    print("ERROR: negative rate coefficient in %s cooling, collider %s!" % (species, collider))
                    print(kul.min(), klu.min())
                    plt.clf()
                    plt.loglog(trange, kul)
                    plt.loglog(trange, klu)
                    plt.show()
                data["rates_m"][collider][u, l] = interp1d(np.log10(trange), np.log10(kul))
                data["rates_m"][collider][l, u] = interp1d(np.log10(trange), np.log10(klu))

    return data, nmultiplets


def lamda_cooling(fname):
    rows = [x.strip().replace("\t", " ") for x in open(fname)]
    rows = [x if x.startswith("!") else x.split("!")[0].strip() for x in rows]

    flags = ["BETWEEN", "TEMPS", "COLLRATES", "NUMBER OF COLL TRANS", "WEIGHT +",
             "EINSTEINA", "NUMBER OF ENERGY LEVELS", "NUMBER OF RADIATIVE TRANSITIONS"]

    flag = ""
    data = dict()
    collider = temps = nlevels = None
    for row in rows:
        while "  " in row:
            row = row.replace("  ", " ")
        for fl in flags:
            if fl in row and row.startswith("!"):
                flag = fl
        if row.startswith("!") and flag not in row:
            flag = ""
        if not row.startswith("!") and flag != "":
            if flag == "NUMBER OF RADIATIVE TRANSITIONS":
                data["nlevels"] = nlevels = int(row)
            if flag == "NUMBER OF ENERGY LEVELS":
                data["nlevels"] = nlevels = int(row)
            if flag == "WEIGHT +":
                if "weights" not in data:
                    data["weights"] = []
                    data["deltaE"] = []
                _, deltaE, g, J = row.split()
                data["weights"].append(float(g))
                data["deltaE"].append(hplanck * clight * float(deltaE))  # cm-1 -> erg

            if flag == "EINSTEINA":
                if "Aul" not in data:
                    data["Aul"] = np.zeros((nlevels, nlevels))
                _, up, low, Aul, _, _ = row.split()
                up = int(up) - 1
                low = int(low) - 1
                data["Aul"][up, low] = float(Aul)

            if flag == "BETWEEN":
                temps = None
                collider = row.split()[-1]
                data[collider] = dict()
            if flag == "TEMPS":
                temps = np.array([float(x) for x in row.split()])
            if flag == "COLLRATES":
                if "rates" not in data[collider]:
                    data[collider]["rates"] = np.empty((nlevels, nlevels), dtype=object)
                _, up, low, k = row.split(" ", maxsplit=3)
                up = int(up) - 1
                low = int(low) - 1
                kul = np.array([float(x) for x in k.split()]) + 1e-40
                delta = data["deltaE"][up] - data["deltaE"][low]
                klu = kul * data["weights"][up] / data["weights"][low] \
                    * np.exp(-delta / kboltzmann / temps) + 1e-40
                t = np.log10(temps)
                data[collider]["rates"][up, low] = interp1d(t, np.log10(kul), fill_value="extrapolate")
                data[collider]["rates"][low, up] = interp1d(t, np.log10(klu), fill_value="extrapolate")

    return data


def cool(data, colliders, xcolliders, tgas):
    nlevels = len(data["deltaE"])
    M = np.zeros((nlevels, nlevels))

    log_tgas = np.log10(tgas)
    for i in range(nlevels):
        for j in range(nlevels):
            M[i, j] += data["Aul"][j, i]
            M[i, i] -= data["Aul"][i, j]
            for c, y in zip(colliders, xcolliders):
                if c not in data:
                    continue

                if data["rates"][c][i, j] is None:
                    continue
                k = 1e1**data["rates"][c][j, i](log_tgas)
                M[i, j] += y * k

                k = 1e1**data["rates"][c]["rates"][i, j](log_tgas)
                M[i, i] -= y * k

    M[-1, :] = 1e0
    b = np.zeros(nlevels)
    b[-1] = 1e0
    x = np.linalg.solve(M, b)

    cool_tot = 0e0
    for i in range(nlevels):
        for j in range(nlevels):
            deltaE = data["deltaE"][i] - data["deltaE"][j]
            cool_tot += data["Aul"][i, j] * x[i] * deltaE

    return cool_tot


