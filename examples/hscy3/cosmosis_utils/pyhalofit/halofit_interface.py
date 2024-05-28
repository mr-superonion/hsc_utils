import pyhmcode
from cosmosis.datablock import option_section, names

iDv_Mead = 4
flag_ucold = 3
iDv_HMcode2016 = 3

hmod_pars = names.halo_model_parameters
cosmo_pars = names.cosmological_parameters


def setup(options):
    config = {}
    config["verbose"] = options.get_int(option_section, "verbose", 1)
    ver = options.get_string(option_section, "version", default="HMcode2020")
    config["version"] = ver
    # if not config["verbose"]:
    # logging.getLogger("INIT_HALOMOD").setLevel(logging.ERROR)
    return config


def setup_hmcode2016(hmod):
    hmod.ihm = 51
    hmod.iDv = iDv_HMcode2016
    hmod.ip2h = 1
    hmod.i1hdamp = 1
    hmod.idc = 3
    hmod.iconc = 1
    hmod.ieta = 1
    hmod.iAs = 1
    hmod.i2hdamp = 2
    hmod.itrans = 1
    hmod.iDolag = 2
    hmod.zD = 10.0
    hmod.DMONLY_neutrino_halo_mass_correction = False
    hmod.Dv0 = 418.0
    hmod.Dv1 = -0.352
    hmod.dc0 = 1.59
    hmod.dc1 = 0.0314
    hmod.As = 3.13
    hmod.eta0 = 0.603
    hmod.eta1 = 0.300
    hmod.f0 = 0.0095
    hmod.f1 = 1.37
    hmod.ks = 0.584
    hmod.alp0 = 3.24
    hmod.alp1 = 1.85
    hmod.Dvnu = 0.916
    hmod.dcnu = 0.262
    # Cold un-normalised produces better massive-neutrino results
    hmod.flag_sigma = flag_ucold
    hmod.mmin = 1e2  # lower-mass limit
    hmod.mmax = 1e17  # upper-mass limit
    hmod.n = 129  # number of points in mass
    return


def setup_hmcode2020(hmod, feedback=False):
    hmod.ip2h = 3  # 3 - Linear two-halo term with damped wiggles
    hmod.i1hdamp = 3  # 3 - k^4 at large scales for one-halo term
    hmod.itrans = 1  # 1 - HMcode alpha-neff smoothing
    hmod.i2hdamp = 3  # 3 - fdamp for perturbation theory
    hmod.idc = 4  # 4 - delta_c from Mead (2017) fit
    hmod.iDv = iDv_Mead
    hmod.iconc = 1
    hmod.iDolag = 3  # 3 - Dolag c(M) correction with sensible z evolution
    hmod.iAs = 2  # 2 - Vary c(M) relation prefactor with sigma8 dependence
    hmod.ieta = 3  # 3 - HMcode 2020 eta bloating
    hmod.zD = 10.0  # z=10 vs 100 does make a difference for EDE-type cosmologies
    hmod.flag_sigma = (
        flag_ucold  # Cold un-normalised produces better massive-neutrino results
    )
    hmod.DMONLY_neutrino_halo_mass_correction = (
        True  # Correct haloes for missing neutrino mass
    )
    hmod.mmin = 1e2  # lower-mass limit
    hmod.mmax = 1e17  # upper-mass limit
    hmod.n = 129  # number of points in mass
    if not feedback:
        hmod.f0 = 0.2695822
        hmod.f1 = 0.9403087
        hmod.ks = 0.0561778
        hmod.kp = -1.0131066
        hmod.kd = 0.0569871
        hmod.kdp = -1.0890162
        hmod.nd = 2.8534197
        hmod.eta0 = 0.1281210
        hmod.eta1 = -0.3643963
        hmod.alp0 = 1.8751465
        hmod.alp1 = 1.6029913
        hmod.As = 5.1958429
        hmod.ihm = 123
    else:
        hmod.DMONLY_baryon_recipe = True
        hmod.response_baseline = 123
        hmod.response_denominator = 125
        hmod.ihm = 124
    return


def execute(block, config):

    pk_lin = block[names.matter_power_lin, "p_k"]
    k_h = block[names.matter_power_lin, "k_h"]
    z = block[names.matter_power_lin, "z"]

    # setup cosmology
    c = pyhmcode.Cosmology()
    c.om_m = block[cosmo_pars, "omega_m"]
    c.om_b = block[cosmo_pars, "omega_b"]
    c.om_v = block[cosmo_pars, "omega_lambda"]
    c.h = block[cosmo_pars, "h0"]
    c.ns = block[cosmo_pars, "n_s"]
    c.sig8 = block[cosmo_pars, "sigma_8"]
    c.m_nu = block[cosmo_pars, "mnu"]
    c.set_linear_power_spectrum(k_h, z, pk_lin)

    # setup hmcode
    ver = config["version"]
    if "HMcode" in ver:
        if ver in ["HMcode2016", "HMcode2016_1par"]:
            hmod = pyhmcode.Halomodel(pyhmcode.HMcode2016, verbose=False)
            setup_hmcode2016(hmod)
            if "1par" in ver:
                eta0 = 0.98 - 0.12 * block[hmod_pars, "A_bary"]
            else:
                eta0 = block[hmod_pars, "eta"]
            hmod.eta0 = eta0
            hmod.As = block[hmod_pars, "A_bary"]
        elif ver == "HMcode2020":
            hmod = pyhmcode.Halomodel(pyhmcode.HMcode2020, verbose=False)
            setup_hmcode2020(hmod, feedback=True)
        elif ver == "HMcode2020_feedback":
            hmod = pyhmcode.Halomodel(pyhmcode.HMcode2020_feedback, verbose=False)
            setup_hmcode2020(hmod, feedback=True)
            hmod.theat = 10.0 ** block[hmod_pars, "logT_AGN"]
        else:
            raise ValueError("version cannot be %s" % config["mode"])
    else:
        raise ValueError("version cannot be %s" % config["mode"])

    pk_nl = pyhmcode.calculate_nonlinear_power_spectrum(
        c,
        hmod,
        verbose=config["verbose"],
    )
    block.put_grid(
        names.matter_power_nl,
        "z",
        z,
        "k_h",
        k_h,
        "p_k",
        pk_nl,
    )
    return 0


def cleanup(config):
    pass
