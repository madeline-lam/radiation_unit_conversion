import unittest
import numpy as np
import astropy.units as u
import radiation_unit_conversion.units as units

test_flux = np.linspace(1, 10, 10)
test_wavelength = np.linspace(1, 5, 10)


class TestUnitsConversions(unittest.TestCase):


    def setUp(self):
        """Set up test data"""
        self.test_flux = np.linspace(1, 10, 10)  # Array of flux values from 1 to 10 W/m^2
        self.test_wavelength = np.linspace(1, 5, 10)  # Array of wavelength values from 1 to 5 micron

    def test_watt_metersquared2erg_cmsquared_second(self):
        """Test conversion from W/m^2 to erg/cm^2/s"""
        expected_output = self.test_flux * 1000  # erg/cm^2/s
        result = units.watt_metersquared2erg_cmsquared_second(self.test_flux)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_erg_cmsquared_second2watt_metersquared(self):
        """Test conversion from erg/cm^2/s to W/m^2"""
        test_flux_erg = self.test_flux * 1000  # erg/cm^2/s
        expected_output = self.test_flux  # W/m^2
        result = units.erg_cmsquared_second2watt_metersquared(test_flux_erg)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_round_trip_conversion(self):
        """Test that a round-trip conversion returns the original values."""
        intermediate_flux = units.watt_metersquared2erg_cmsquared_second(self.test_flux)
        final_flux = units.erg_cmsquared_second2watt_metersquared(intermediate_flux)
        np.testing.assert_allclose(self.test_flux, final_flux, rtol=1e-6)

    def test_watt_metersquared_hertz2erg_cmsquared_second_hertz(self):
        """Test conversion from W/m^2/Hz to erg/cm^2/s/Hz"""
        expected_output = self.test_flux * 1000
        result = units.watt_metersquared_hertz2erg_cmsquared_second_hertz(self.test_flux)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_erg_cmsquared_secondhertz2watt_metersquaredhertz(self):
        """Test conversion from erg/cm^2/s/Hz to W/m^2/Hz"""
        test_flux_erg = self.test_flux * 1000
        expected_output = self.test_flux
        result = units.erg_cmsquared_secondhertz2watt_metersquaredhertz(test_flux_erg)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_watt_metersquared_hertz2erg_cmsquared_second_angstrom(self):
        """Test conversion from W/m^2/Hz to erg/cm^2/s/\u00c5ngström"""
        constant = 2.99792458e21
        expected_output = constant * self.test_flux / self.test_wavelength**2
        result = units.watt_metersquared_hertz2erg_cmsquared_second_angstrom(self.test_flux, self.test_wavelength)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_erg_cmsquared_second_angstrom2watt_metersquared_hertz(self):
        """Test conversion from erg/cm^2/s/\u00c5ngström to W/m^2/Hz"""
        constant = 2.99792458e21
        test_flux_erg = constant * self.test_flux / self.test_wavelength**2
        expected_output = self.test_flux
        result = units.erg_cmsquared_second_angstrom2watt_metersquared_hertz(test_flux_erg, self.test_wavelength)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_watt_metersquared_micron2watt_metersquared_hertz(self):
        """Test conversion from W/m^2/\u03bcm to W/m^2/Hz"""
        constant = 2.99792458e14
        expected_output = constant * self.test_flux / self.test_wavelength**2
        result = units.watt_metersquared_micron2watt_metersquared_hertz(self.test_flux, self.test_wavelength)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_watt_metersquared_hertz2watt_metersquared_micron(self):
        """Test conversion from W/m^2/Hz to W/m^2/\u03bcm"""
        constant = 2.99792458e14
        test_flux_hz = constant * self.test_flux / self.test_wavelength**2
        expected_output = self.test_flux
        result = units.watt_metersquared_hertz2watt_metersquared_micron(test_flux_hz, self.test_wavelength)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_photon_cmsquared_second_micron2watt_metersquared_micron(self):
        """Test conversion from photon/cm^2/s/μm to W/m^2/μm"""
        constant = 5.03411250e14
        expected_output = self.test_flux / self.test_wavelength / constant
        result = units.photon_cmsquared_second_micron2watt_metersquared_micron(self.test_flux, self.test_wavelength)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_erg_cmsquared_second_angstrom2photon_cmsquared_second_angstrom(self):
        """Test conversion from erg/cm^2/s/Ångström to photon/cm^2/s/Ångström"""
        constant = 5.03411250E+07
        expected_output = constant * self.test_flux * self.test_wavelength
        result = units.erg_cmsquared_second_angstrom2photon_cmsquared_second_angstrom(self.test_flux, self.test_wavelength)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_photon_cmsquared_second_angstrom2erg_cmsquared_second_angstrom(self):
        """Test conversion from photon/cm^2/s/Ångström to erg/cm^2/s/Ångström"""
        constant = 5.03411250E+07
        expected_output = self.test_flux / self.test_wavelength / constant
        result = units.photon_cmsquared_second_angstrom2erg_cmsquared_second_angstrom(self.test_flux, self.test_wavelength)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_watt_metersquared_hertz2jansky(self):
        """Test conversion from W/m^2/Hz to Jansky"""
        constant = 1e26
        expected_output = self.test_flux * constant
        result = units.watt_metersquared_hertz2jansky(self.test_flux)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_jansky2watt_metersquared_hertz(self):
        """Test conversion from Jansky to W/m^2/Hz"""
        constant = 1e26
        test_flux_jansky = self.test_flux * constant
        expected_output = self.test_flux
        result = units.jansky2watt_metersquared_hertz(test_flux_jansky)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_erg_cmsquared_second_hertz2jansky(self):
        """Test conversion from erg/cm^2/s/Hz to Jansky"""
        constant = 1e23
        expected_output = self.test_flux * constant
        result = units.erg_cmsquared_second_hertz2jansky(self.test_flux)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_jansky2erg_cmsquared_second_hertz(self):
        """Test conversion from Jansky to erg/cm^2/s/Hz"""
        constant = 1e23
        test_flux_jansky = self.test_flux * constant
        expected_output = self.test_flux
        result = units.jansky2erg_cmsquared_second_hertz(test_flux_jansky)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_erg_cmsquared_second_angstrom2jansky(self):
        """Test conversion from erg/cm^2/s/Ångström to Jansky"""
        constant = 3.33564095e04
        expected_output = constant * self.test_flux * self.test_wavelength**2
        result = units.erg_cmsquared_second_angstrom2jansky(self.test_flux, self.test_wavelength)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_jansky2erg_cmsquared_second_angstrom(self):
        """Test conversion from Jansky to erg/cm^2/s/Ångström"""
        constant = 3.33564095e04
        test_flux_jansky = constant * self.test_flux * self.test_wavelength**2
        expected_output = self.test_flux
        result = units.jansky2erg_cmsquared_second_angstrom(test_flux_jansky, self.test_wavelength)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_watt_metersquared_micron2jansky(self):
        """Test conversion from W/m^2/μm to Jansky"""
        constant = 3.33564095e04
        expected_output = constant * self.test_flux * self.test_wavelength**2
        result = units.watt_metersquared_micron2jansky(self.test_flux, self.test_wavelength)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_jansky2watt_metersquared_micron(self):
        """Test conversion from Jansky to W/m^2/μm"""
        fnu = units.watt_metersquared_micron2watt_metersquared_hertz(self.test_flux, self.test_wavelength)
        expected_output = units.watt_metersquared_hertz2jansky(fnu)
        result = units.jansky2watt_metersquared_micron(self.test_flux, self.test_wavelength)
        np.testing.assert_allclose(result, expected_output, rtol=1e-6)

    def test_flambda2fnu(self):
        """Test conversion from Jansky to W/m^2/μm"""
        expected_output = np.array([
            2.99792458e+14,
            2.87375019e+14,
            2.52074627e+14,
            2.20255683e+14,
            1.94265513e+14,
            1.73245107e+14,
            1.56090288e+14,
            1.41903223e+14,
            1.30011125e+14,
            1.19916983e+14]
        )
        result = units.flambda2fnu(self.test_flux * u.W /u.m**2 / u.micron, self.test_wavelength * u.micron,  u.W /u.m**2 / u.Hz)
        np.testing.assert_allclose(result.value, expected_output, rtol=1e-6)


if __name__ == "__main__":
    unittest.main()