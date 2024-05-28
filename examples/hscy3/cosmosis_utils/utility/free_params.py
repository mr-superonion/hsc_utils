from cosmosis.datablock import names, option_section as opt
import numpy as np

cosmo = names.cosmological_parameters
hmod_pars = names.halo_model_parameters


def setup(options):
    sample_choice = options.get_string(opt, "sample_choice", default="")
    return sample_choice

def execute(block, config):
    sample_choice = config
    if "lnAs" in sample_choice:
        block[cosmo, "A_s"] = np.exp(block[cosmo, "ln_as1e10"]) * 1e-10
    elif "logAs" in sample_choice:
        block[cosmo, "A_s"] = pow(10, block[cosmo, "log_as1e9"]) * 1e-9
    elif "S_8" in sample_choice:
        block[cosmo, "sigma_8"] =  block[cosmo, "S_8"] / (
                        block[cosmo, "omega_m"] / 0.3
                        )**0.5
    if "camb" in sample_choice:
        block[hmod_pars, "eta"] = 0.98 - 0.12*block[hmod_pars, "A_bary"]
        block[hmod_pars, "A"] = block[hmod_pars, "A_bary"]
    if sample_choice == "kidsdes":
        block[cosmo, "sigma_8_input"] =  block[cosmo, "S_8"] / (
                        block[cosmo, "omega_m"] / 0.3
                        )**0.5
    return 0


def cleanup(config):
    # nothing to do here!  We just include this
    # for completeness
    return 0
