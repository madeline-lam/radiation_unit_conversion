=========================
radiation_unit_conversion
=========================

**An easy-to-use Python package for converting between different radiation units**

Implementing the `STScI units page <https://www.stsci.edu/~strolger/docs/UNITS.txt>`_, since astropy doesn't like to easily convert between radiation units.

Example
=======
Using astropy units: ::
  import numpy as np
  import astropy.units as u
  import radiation_unit_conversion.units as units

  flambda_spectrum = np.linspace(1,10,10) * u.W/u.m**2/u.micron
  wavelength = np.linspace(1,5,10) * u.micron
  fnu_spectrum = units.flambda2fnu(flambda_spectrum, wavelength, u.Jy)

Without using astropy units: ::
  import numpy as np
  import radiation_unit_conversion.units as units

  flambda_spectrum = np.linspace(1,10,10) # W/m^2/micron
  wavelength = np.linspace(1,5,10) # 
  fnu_spectrum = units.watt_metersquared_micron2jansky(flambda_spectrum, wavelength)

License
=======
Copyright 2024 Evert Nasedkin

radiation_unit_conversion is available under the MIT license.
See the LICENSE file for more information.
