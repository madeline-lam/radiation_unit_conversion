# For debugging: set to debug to 'true' and optimization to 'g'

project('radiation_unit_conversion',
    'c',
    # version must be updated here as well
    version : '0.1',
    license: 'MIT',
    meson_version: '>=1.1.0',
    default_options : [
        'debug=false',
        'optimization=3',
        'warning_level=1',
        'c_std=c99'
    ],
)

py_mod = import('python')
py = py_mod.find_installation(pure: false)
py_dep = py.dependency()

# Install python sources
subdir('radiation_unit_conversion')
