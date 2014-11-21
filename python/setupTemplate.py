from distutils.core import setup, Extension
from glob import glob

module1 = Extension('myModuleLNA',
                    define_macros = [('MAJOR_VERSION', '1'),
                                     ('MINOR_VERSION', '0')],
                    include_dirs = ['/usr/local/include', '/opt/local/Library/Frameworks/Python.framework/Versions/3.4/include/python3.4m/', 'models/myModule/C', \
					'/opt/local/Library/Frameworks/Python.framework/Versions/3.4/lib/python3.4/site-packages/numpy/core/include/', 'C++'],
                    libraries = ['stdc++.6', 'sundials_cvodes', 'blitz', 'sundials_nvecserial','python3.4'],
                    library_dirs = ['/usr/local/lib','/opt/local/Library/Frameworks/Python.framework/Versions/3.4/lib/'],
                    sources = ['C++/computeLinearNoise.cpp', 'models/myModule/myModule_LNA.cpp']\
						+ glob('models/myModule/C/*.c'), 
					runtime_library_dirs=['/usr/local/lib'])
					

setup(name = 'myModuleLNA',
       version = '1.0',
       description = '',
       author = 'Justin Feigelman',
       author_email = 'justin.feigelman@helmholtz-muenchen.de',
       url = '',
       long_description='',
       ext_modules = [module1])

