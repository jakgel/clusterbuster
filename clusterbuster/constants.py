#!/usr/bin/env python

class ConstantsCGS :
   def __init__(self):
      self.gravity                = 6.672e-8       # cm3 g-1 s-2
      self.speed_of_light         = 2.997925e10    # cm s-1
      self.planck                 = 6.626176e-27   # erg s-1
      self.electron_charge        = 4.80325e-10    # esu
      self.boltzmann              = 1.38066e-16    # erg K-1
      self.stefan_boltzmann       = 5.6703e-5      # erg cm-2 K-4 s-1
      self.radiation_dens         = 7.5646e-15     # erg cm-3 K-4
      self.gas_constant           = 8.31425e7      # g cm3 s-2 K-1
      self.rydberg_energy         = 2.17992e-11    # erg
      self.rydberg_frequency      = 3.289842e15    # s-1
      self.rydberg_wave_vector    = 1.09737e5      # cm-1
      self.bohr_radius            = 5.29177e-9     # cm
      self.thomson_cross_section  = 6.6524616e-25  # cm2
      self.proton_mass            = 1.67261e-24    # g
      self.neutron_mass           = 1.67482e-24    # g
      self.hydrogen_mass          = 1.6747e-24     # g
      self.electron_mass          = 9.1096e-28     # g
      self.atomic_mass_unit       = 1.66055e-24    # g
      self.parsec                 = 3.085678e18    # cm
      self.light_year             = 9.463e17       # cm
      self.astronomical_unit      = 1.496e13       # cm
      self.year                   = 3.156e7        # s
      self.electron_volt          = 1.6022e-12     # erg
      self.barn                   = 1.0e-24        # cm2
      self.coulomb                = 2.99796e9      # esu
      self.kilometer_per_sec      = 1.0e5          # cm s-1
      self.solar_mass             = 1.989e33       # g
      self.solar_lum              = 3.826e33       # erg s-1
      self.solar_radius           = 9.96e10        # cm
      self.solar_temperature      = 5.780e3        # K
      self.solar_metallizity      = 0.02           #  estimate !!!
      self.hubble                 = 3.2407789e-18  # h s-1
      self.Y_helium               = 0.24           # relative mass dens
      self.T_cmb_0                = 2.73           # T of CMB in Kelvin
      self.sigma_effect_HI        = 6.3e-18        # cm2, note from Jan
      self.sigma_effect_HeI       = 7.4e-18        # cm2, note from Jan
      self.sigma_effect_HeII      = 1.6e-18        # cm2, note from Jan
      # critical density = 3 H_0 / 8 pi G,  values from Sparke, Gallagher 
      # multiply the values by h^2    
      self.critical_density       = 1.9e-29        # g cm-3
      self.critical_density_cosm  = 2.8e11         # M_sun Mpc-3const 



#  Units, conversion factors, and input variables
#  Many are here for historic reasons and should eventually become swapped out
radiounit_A  =  1                     # W/Hz        --- Unit of particle luminousity in .radio snaps, formerly 1e23 for Sebastians cubes
Jy           =  1e-23                 # erg/s/Hz/cm^2
Jy_SI        =  1e-26                 # W/Hz/m^2
H0           =  70                    # 0.7  ... Hubble constant used in MUSIC-2
FWHM2sigma   =  1/2.354               # convert FWHM to sigma for a gaussian
