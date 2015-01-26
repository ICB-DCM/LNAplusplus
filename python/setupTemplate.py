from distutils.core import setup, Extension
from glob import glob
from numpy import get_include
import sys

# require version 3.0+
if sys.version_info.major < 3:
	print('''LNA++: Please build module using Python version 3.0 or greater.
	Aborting build.''')
	sys.exit(1)

numpyPath = get_include()

module1 = Extension('myModuleLNA',
                    define_macros = [('MAJOR_VERSION', '1'),
                                     ('MINOR_VERSION', '0')],
                    include_dirs = ['/usr/local/include', 'myModule/C', '../include', '../src', numpyPath]+INCLUDE_DIRS,
                    libraries = ['stdc++', 'sundials_cvodes', 'blitz', 'sundials_nvecserial'], #,'python3.4'],
                    library_dirs = ['/usr/local/lib']+LIB_DIRS,
                    sources = ['../src/computeLinearNoise.cpp', 'myModule/myModule_LNA.cpp']\
						+ glob('myModule/C/*.c'), 
					runtime_library_dirs=['/usr/local/lib'],
					extra_compile_args=['-Wno-header-guard','-Wno-parentheses','-Wno-incompatible-pointer-types', '-Wno-unused-variable', 
					'-Wno-c++11-compat-deprecated-writable-strings', '-Wno-unused-function', '-Wno-reorder'])
					

setup(name = 'myModuleLNA',
       version = '1.0',
       description = '',
       author = 'Justin Feigelman',
       author_email = 'justin.feigelman@helmholtz-muenchen.de',
       url = '',
       long_description='',
       ext_modules = [module1])

