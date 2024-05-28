import warnings
import emu
import numpy as np
from cosmosis.datablock import names, option_section as opt

warnings.filterwarnings("ignore", category=UserWarning)
As_try = 2.1841e-9  # temporary A_s


cosmo_pars = names.cosmological_parameters
distance_pars = names.distances
growth_pars = names.growth_parameters

def setup(options):
    zmax = options.get_double(opt, "zmax", default=4.00)
    nz = options.get_int(opt, "nz", default=150)
    # power spectrum
    nk = options.get_int(opt, "nk", default=500)
    # these are k [h/Mpc]
    logkhmin = np.log10(options.get_double(opt, "kmin", default=1e-4))
    logkhmax = np.log10(options.get_double(opt, "kmax", default=50.0))
    cc = emu.PsConfig(
            nz=nz, zmax=zmax,
            logkh_min=logkhmin,
            logkh_max=logkhmax,
            nk=nk,
            )
    emulator = emu.BaccoEmuLin(cc)
    sample_param = options.get_string(
            opt, "sample_param", default="A_s"
            )
    return emulator, sample_param


def save_derived_parameters(cdict, block):
    for cosmosis_name, cname, scaling in [
        ("h0", "h_0", 1),
        ("hubble", "h_0", 100),
        ("omnuh2", "omnuh2", 1),
        ("ommh2", "ommh2", 1),
        ("omega_nu", "Omega_nu", 1),
        ("omega_b", "Omega_b", 1),
        ("omega_c", "Omega_c", 1),
        ("omega_m", "Omega_m", 1),
        ("omega_lambda", "Omega_v", 1),
    ]:
        cvale = cdict[cname] * scaling
        if block.has_value(cosmo_pars, cosmosis_name):
            input_value = block[cosmo_pars, cosmosis_name]
            if not np.isclose(input_value, cvale, rtol=2e-3):
                warnings.warn(
                    f"Parameter {cosmosis_name} inconsistent: \
                        input was {input_value} but value is now {cvale}."
                )
        block[cosmo_pars, cosmosis_name] = cvale
    return


def save_matter_power_lin(emu, pklin_table, block):
    # Linear version of the spectrum
    section_name = "matter_power_lin"
    sigma8 =  emu.get_sigma8(pklin_table[0, :])
    if block.has_value(cosmo_pars, "sigma_8"):
        assert np.isclose(sigma8, block[cosmo_pars, "sigma_8"])
    block[cosmo_pars, "sigma_8"] = sigma8
    block[cosmo_pars, "S_8"] = sigma8 * np.sqrt(block[cosmo_pars, "omega_m"] / 0.3)
    block.put_grid(
            section_name, "z", emu.z, "k_h", emu.ks,
            "p_k", pklin_table,
            )
    return

def get_input_pars(block):
    omv = 1.0 - block[cosmo_pars, "Omega_m"] - block[cosmo_pars, "Omega_k"]
    params = {
        'A_s': block[cosmo_pars, "A_s"],
        'n_s': block[cosmo_pars, "n_s"],
        'Omega_b': block[cosmo_pars, "Omega_b"],
        'Omega_m': block[cosmo_pars, "Omega_m"],
        'Omega_v': omv,
        'h_0': block[cosmo_pars, "h0"],
        "mnu": block.get_double(names.cosmological_parameters, "mnu", 0.06),
        "w": block.get_double(names.cosmological_parameters, "w", -1.0),
        "wa": block.get_double(names.cosmological_parameters, "wa", 0.0),
    }
    return params

def execute(block, config):
    # preparation
    emulator, sample_param = config
    if sample_param == "sigma_8":
        block[cosmo_pars, "A_s"] = As_try
    params = get_input_pars(block)

    # some derived parameters
    cpars = emu.make_pars_dict(params)
    save_derived_parameters(cpars, block)

    # linear power spectrum
    epars = emulator.get_emu_params(params)
    pklin_table = emulator.compute_pklin_table(epars)
    if sample_param == "sigma_8":
        sigma8_try =  emulator.get_sigma8(pklin_table[0, :])
        rs = (block[cosmo_pars, "sigma_8"] / sigma8_try) ** 2
        pklin_table = pklin_table * rs
        block[cosmo_pars, "A_s"] = As_try * rs

    save_matter_power_lin(emulator, pklin_table, block)
    return 0


def cleanup(config):
    pass
