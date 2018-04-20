from distutils.core import setup, Extension
from distutils import sysconfig
from glob import glob
from numpy import get_include
import sys

# require version 3.0+
if sys.version_info.major < 3:
	print('''LNA++: Please build module using Python version 3.0 or greater.
	Aborting build.''')
	sys.exit(1)

numpyPath = get_include()
lnaRootDir = 'LNA_ROOT_DIR'

# Remove the "-Wstrict-prototypes" compiler option, which isn't valid for C++ to fix warnings.
cfg_vars = sysconfig.get_config_vars()
for key, value in cfg_vars.items():
    if type(value) == str:
        cfg_vars[key] = value.replace("-Wstrict-prototypes", "")
 
module1 = Extension('myModuleLNA',
                    define_macros = [('MAJOR_VERSION', '1'),
                                     ('MINOR_VERSION', '0')],
                    include_dirs = INCLUDE_DIRS + ['/usr/local/include', 
												lnaRootDir + '/myModule/C', 
												lnaRootDir + '/include', 
												lnaRootDir + '/src', 
												lnaRootDir + '/libraries/blitz-1.0.1/install/include', 
												lnaRootDir + '/libraries/cvodes-2.7.0/install/include', 
												numpyPath],
                    libraries = ['stdc++', 'sundials_cvodes', 'blitz', 'sundials_nvecserial'], #,'python3.4'],
                    library_dirs = LIB_DIRS + ['/usr/local/lib', 
											lnaRootDir + '/libraries/blitz-1.0.1/install/lib',
											lnaRootDir + '/libraries/cvodes-2.7.0/install/lib/'],
                    sources = [lnaRootDir + '/src/computeLinearNoise.cpp', lnaRootDir + '/myModule/myModule_LNA.cpp']\
						+ glob(lnaRootDir + '/myModule/C/*.c'), 
					runtime_library_dirs=['/usr/local/lib'],
					extra_compile_args=['-Wno-parentheses', '-Wno-unused-variable', '-Wno-unused-function'])
					

setup(name = 'myModuleLNA',
       version = '1.0',
       description = '',
       author = 'Justin Feigelman',
       author_email = 'justin.feigelman@helmholtz-muenchen.de',
       url = '',
       long_description='',
       ext_modules = [module1])

