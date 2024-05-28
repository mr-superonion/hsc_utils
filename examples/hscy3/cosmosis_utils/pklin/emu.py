import baccoemu
import numpy as np
from typing import Union
from dataclasses import dataclass

inv_neutrino_mass_fac = 0.010752306648969645
def make_pars_dict(cdict):
    omnuh2 = cdict["mnu"] * inv_neutrino_mass_fac
    h2 = cdict["h_0"] ** 2.0
    omnu = omnuh2 / h2
    omc = (cdict["Omega_m"] - cdict["Omega_b"] - omnu)
    omch2 =  omc * h2
    omk = 1.0 - cdict["Omega_m"] - cdict["Omega_v"]
    ln10p10As = float(np.log(cdict["A_s"] * 1.0e10))
    odict = cdict.copy()
    odict["ombh2"] = cdict["Omega_b"] * h2
    odict["omch2"] = omch2
    odict["omnuh2"] = omnuh2
    odict["ommh2"] = cdict["Omega_m"]* h2
    odict["Omega_c"] = omc
    odict["Omega_k"] = omk
    odict["Omega_nu"] = omnu
    odict["ln10p10As"] = ln10p10As
    assert np.isclose(
            odict["ommh2"], odict["omch2"] + odict["ombh2"] + odict["omnuh2"],
            )
    assert np.isclose(
            odict["Omega_m"],
            odict["Omega_c"] + odict["Omega_b"] + odict["Omega_nu"],
            )
    return odict

@dataclass(frozen=True)
class GeoConfig:
    zmax: Union[float, int] = 4.0
    nz: int = 401

@dataclass(frozen=True)
class PsConfig(GeoConfig):
    logkh_min: Union[float, int] = -4
    logkh_max: Union[float, int] = 1.5
    nk: int = 700

def window_tophat(kR):
    return np.float64(3.0 * (np.sin(kR) - kR * np.cos(kR)) / kR**3.0)

class BaccoEmuLin:
    def __init__(self, config):
        if not isinstance(config, PsConfig):
            raise TypeError("input config is not a power spectrum configuration")
        self.ks = np.logspace(
            config.logkh_min,
            config.logkh_max,
            config.nk,
        )  # [h/Mpc]
        self.z = np.linspace(
            0.0,
            config.zmax,
            config.nz,
        )
        self.a = 1.0 / (1 + self.z)
        # dark emulator based linear class
        self.pkL = baccoemu.Matter_powerspectrum(
            nonlinear_boost=False,
            baryonic_boost=False,
            verbose=False,
        )
        return

    def get_emu_params(self, lss_pars):
        """Gets the parameters

        Args:
          lss_pars (dict) : cosmology dictionary
        """
        if not np.isclose(1.0 - lss_pars["Omega_v"] - lss_pars["Omega_m"],
                0.0, atol = 1.e-5):
            raise ValueError("BaccoEmu does not support Omega_k != 0")
        # init dark emulator
        params = {
            "omega_matter": lss_pars["Omega_m"],
            "A_s": lss_pars["A_s"],
            "omega_baryon": lss_pars["Omega_b"],
            "ns": lss_pars["n_s"],
            "hubble": lss_pars["h_0"],
            "neutrino_mass": lss_pars["mnu"],
            "w0": lss_pars["w"],
            "wa": lss_pars["wa"],
            "expfactor": self.a,
        }
        return params

    def compute_pklin_table(self, params):
        """Computes the linear power spectrum table.
        rows of the table: z; columns of the table: k
        """
        _, pklin_table = self.pkL.get_linear_pk(k=self.ks, cold=False, **params)
        return pklin_table

    def get_sigmaR(self, R, pklin0):
        """Calculates `sigma_R`, the RMS linear matter fluctuation in
        spheres of radius R

        Args:
            R (float):              radius of the spheres [Mpc/h]
            pklin0 (ndarray):       linear power spetrum at z=0
            nint (int)              number of points for integral [default: 100]
        """
        assert pklin0.shape == self.ks.shape
        logks = np.log(self.ks)
        kR = self.ks * R
        integrant = self.ks**3.0 * pklin0 * window_tophat(kR) ** 2
        sigma = np.sqrt(np.trapz(integrant, logks) / (2.0 * np.pi**2))
        return sigma

    def get_sigma8(self, pklin0):
        """Calculates the RMS linear matter fluctuation in spheres of radius 8
        Mpc/h
        """
        return self.get_sigmaR(8., pklin0)
