import astropy.constants as cst
import astropy.units as u
"""
radiation_unit_conversion
This module is a python implementation of the radiation unit
conversion functions found at https://www.stsci.edu/~strolger/docs/UNITS.txt.

The fnu2flambda and flambda2fnu functions take arbitrary astropy radiation units
and will convert it to the specified output. Any functions with a wavelength
dependence require wavelength units rather than flux units, though conversion
functions are also provided.
"""


def fnu2flambda(input_flux, input_wavelength, output_units):
    """
    Converts flux from the f_nu (flux per unit frequency) form to the f_lambda 
    (flux per unit wavelength) form.

    Parameters:
    -----------
    input_flux : `astropy.units.Quantity`
        Flux in units of W/m^2/Hz or similar frequency-based units (e.g., erg/s/cm^2/Hz).

    input_wavelength : `astropy.units.Quantity`
        Wavelength in units of meters, or similar wavelength-based units (e.g., microns, nanometers).

    output_units : `astropy.units.Unit`
        The desired output units for the flux (e.g., W/m^2/micron, erg/s/cm^2/angstrom).

    Returns:
    --------
    `astropy.units.Quantity`
        The converted flux in the desired output units (e.g., W/m^2/micron, erg/s/cm^2/angstrom).

    Notes:
    ------
    The conversion from f_nu to f_lambda uses the relationship:
        f_lambda = f_nu * lambda^2 / c
    where `lambda` is the wavelength and `c` is the speed of light.
    """
    # Convert input flux to W/m^2/Hz
    input_converted = input_flux.to(u.W/u.m**2/u.Hz)

    # Convert wavelength to microns
    wlen_converted = input_wavelength.to(u.micron)

    # Define the speed of light constant (in microns per second)
    constant = cst.c.to(u.micron/u.second)

    # Convert the flux from f_nu to f_lambda using the formula
    output = input_converted.value * wlen_converted.value**2 / constant.value

    # Multiply by the appropriate units for flux per wavelength
    output *= u.W / u.m**2 / u.micron

    # Convert to the specified output units
    output = output.to(output_units)

    return output


def flambda2fnu(input_flux, input_wavelength, output_units):
    """
    Converts flux from the f_lambda (flux per unit wavelength) form to the f_nu 
    (flux per unit frequency) form.

    Parameters:
    -----------
    input_flux : `astropy.units.Quantity`
        Flux in units of W/m^2/micron or similar wavelength-based units (e.g., erg/s/cm^2/angstrom).

    input_wavelength : `astropy.units.Quantity`
        Wavelength in units of meters, or similar wavelength-based units (e.g., microns, nanometers).

    output_units : `astropy.units.Unit`
        The desired output units for the flux (e.g., W/m^2/Hz, erg/s/cm^2/Hz).
    Returns:
    --------
    `astropy.units.Quantity`
        The converted flux in the desired output units (e.g., W/m^2/Hz, erg/s/cm^2/Hz).

    Notes:
    ------
    The conversion from f_lambda to f_nu uses the relationship:
        f_nu = c * f_lambda / lambda^2
    where `lambda` is the wavelength and `c` is the speed of light.
    """
    # Convert input flux to W/m^2/micron
    input_converted = input_flux.to(u.W / u.m**2 / u.micron)

    # Convert wavelength to microns
    wlen_converted = input_wavelength.to(u.micron)

    # Define the speed of light constant (in microns per second)
    constant = cst.c.to(u.micron / u.second)

    # Convert the flux from f_lambda to f_nu using the formula
    output = constant.value * input_converted.value / wlen_converted.value**2

    # Multiply by the appropriate units for flux per frequency
    output *= u.W / u.m**2 / u.Hz

    # Convert to the specified output units
    output = output.to(output_units)

    return output


def wavelength2frequency(input_wavelength, output_units):
    input_wavelength.to(u.m)
    return (cst.c / input_wavelength).to(output_units)


def frequency2wavelength(input_frequency, output_units):
    input_frequency.to(u.m)
    return (cst.c / input_frequency).to(output_units)


## STScI direct conversion functions
def watt_metersquared2erg_cmsquared_second(input_flux):
    # [Y erg/cm^2/s] = 1000 * [X W/m^2]
    return 1000 * input_flux


def erg_cmsquared_second2watt_metersquared(input_flux):
    return input_flux / 1000


def watt_metersquared_hertz2erg_cmsquared_second_hertz(input_flux):
    # [Y erg/cm^2/s/Hz] = 1000 * [X W/m^2/Hz]
    return 1000 * input_flux


def erg_cmsquared_secondhertz2watt_metersquaredhertz(input_flux):
    return input_flux / 1000


def watt_metersquared_hertz2erg_cmsquared_second_angstrom(input_flux, input_wavelength):
    # [Y erg/cm^2/s/A] = 2.99792458E+21 * [X1 W/m^2/Hz] / [X2 A]^2
    constant = 2.99792458e21
    return constant * input_flux / input_wavelength**2


def erg_cmsquared_second_angstrom2watt_metersquared_hertz(input_flux, input_wavelength):
    constant = 2.99792458e21
    return  input_flux * input_wavelength**2 / constant


def watt_metersquared_micron2watt_metersquared_hertz(input_flux, input_wavelength):
    # [Y W/m^2/um] = 2.99792458E+14 * [X1 W/m^2/Hz] / [X2 um]^2
    constant = 2.99792458e14
    return constant * input_flux / input_wavelength**2


def watt_metersquared_hertz2watt_metersquared_micron(input_flux, input_wavelength):
    constant = 2.99792458e14
    return input_flux*input_wavelength**2 / constant


def erg_cmsquared_hertz2photon_cmsquared_second_micron(input_flux, input_wavelength):
    # [Y photon/cm^2/s/um] = 1.50918896E+22 * [X1 erg/cm^2/Hz] / [X2 um]
    constant = 1.50918896e22
    return constant * input_flux/input_wavelength


def photon_cmsquared_second_micron2erg_cmsquared_hertz(input_flux, input_wavelength):
    constant = 1.50918896e22
    return input_flux * input_wavelength / constant 


def watt_metersquared_micron2photon_cmsquared_second_micron(input_flux,input_wavelength):
    # [Y photon/cm^2/s/um] = 5.03411250E+14 * [X1 W/m^2/um] * [X2 um]
    constant = 5.03411250e14
    return constant * input_flux * input_wavelength


def photon_cmsquared_second_micron2watt_metersquared_micron(input_flux, input_wavelength):
    constant = 5.03411250e14
    return input_flux / input_wavelength / constant


def erg_cmsquared_second_angstrom2photon_cmsquared_second_angstrom(input_flux, input_wavelength):
    # [Y photon/cm^2/s/A] = 5.03411250E+07 * [X1 erg/cm^2/s/A] * [X2 A]
    constant = 5.03411250E+07 
    return constant * input_flux * input_wavelength


def photon_cmsquared_second_angstrom2erg_cmsquared_second_angstrom(input_flux, input_wavelength):
    constant = 5.03411250E+07 
    return input_flux / input_wavelength / constant


def watt_metersquared_hertz2jansky(input_flux):
    # [Y Jy] = 1.0E+26 * [X W/m^2/Hz]
    constant = 1e26
    return constant * input_flux


def jansky2watt_metersquared_hertz(input_flux):
    constant = 1e26
    return input_flux / constant


def erg_cmsquared_second_hertz2jansky(input_flux):
    # [Y Jy] = 1.0E+23 * [X erg/cm^2/s/Hz]
    constant = 1e23
    return input_flux * constant


def jansky2erg_cmsquared_second_hertz(input_flux):
    constant = 1e23
    return input_flux / constant


def erg_cmsquared_second_angstrom2jansky(input_flux, input_wavelength):
    # [Y Jy] = 3.33564095E+04 * [X1 erg/cm^2/s/A] * [X2 A]^2
    constant = 3.33564095e04 
    return constant * input_flux * input_wavelength**2


def jansky2erg_cmsquared_second_angstrom(input_flux, input_wavelength):
    constant = 3.33564095e04 
    return input_flux / input_wavelength**2 / constant


def watt_metersquared_micron2jansky(input_flux, input_wavelength):
    # [Y erg/cm^2/s/A] = 2.99792458E-05 * [X1 Jy] / [X2 A]^2
    constant = 3.33564095e04 
    return constant * input_flux * input_wavelength**2


def jansky2watt_metersquared_micron(input_flux, input_wavelength):
    fnu = watt_metersquared_micron2watt_metersquared_hertz(input_flux, input_wavelength)
    return watt_metersquared_hertz2jansky(fnu)


def jansky2photon_cmsquared_second_angstrom(input_flux, input_wavelength):
    # [Y photon/cm^2/s/A] = 1.50918896E+03 * [X1 Jy] / [X2 A]
    constant = 1.50918896e03
    return constant * input_flux / input_wavelength


def photon_cmsquared_second_angstrom2jansky(input_flux, input_wavelength):
    constant = 1.50918896e03
    return input_flux * input_wavelength  / constant


def rayleigh2photon_cmsquared_second_angstrom_steradian(input_flux):
    # [Y photon/cm^2/s/A/sr] = 7.9577539E+04 [X Rayleigh]
    constant = 7.9577539e04 
    return constant * input_flux


def photon_cmsquared_second_angstrom_steradian2rayleigh(input_flux):
    constant = 7.9577539e04 
    return  input_flux / constant


def rayleigh2photon_cmsquared_angstrom_degreesquared(input_flux):
    # [Y photon/cm^2/s/A/deg^2] = 2.4240705E+01 [X Rayleigh]
    constant = 2.4240705e01
    return constant * input_flux


def photon_cmsquared_angstrom_degreesquared2rayleigh(input_flux):
    constant = 2.4240705e01
    return  input_flux / constant


def rayleigh2photon_cmsquared_second_angstrom_arcsecondsquared(input_flux):
    # [Y photon/cm^2/s/A/arcsec^2] = 1.8704247E-06 [X Rayleigh]
    constant = 1.8704247e06 
    return constant * input_flux


def photon_cmsquared_second_angstrom_arcsecondsquared2rayleigh(input_flux):
    constant = 1.8704247e06 
    return constant / input_flux

