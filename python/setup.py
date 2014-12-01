from distutils.core import setup, Extension
from glob import glob

module1 = Extension('birthDeathLNA',
                    define_macros = [('MAJOR_VERSION', '1'),
                                     ('MINOR_VERSION', '0')],
                    include_dirs = ['/usr/local/include', '/opt/local/Library/Frameworks/Python.framework/Versions/3.4/include/python3.4m/', 'C', \
					'/opt/local/Library/Frameworks/Python.framework/Versions/3.4/lib/python3.4/site-packages/numpy/core/include/'], # '/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.10.sdk/System/Library/Frameworks/Python.framework/Versions/2.7/include/python2.7/'],
                    libraries = ['stdc++.6', 'sundials_cvodes', 'blitz', 'sundials_nvecserial','python3.4'],
                    library_dirs = ['/usr/local/lib','/opt/local/Library/Frameworks/Python.framework/Versions/3.4/lib/'],
                    #sources = ['../src/LNAModule.cpp', '../src/computeLinearNoise.cpp','../src/main.cpp']+glob('C/*.c'),
                    sources = ['../src/computeLinearNoise.cpp','../src/pyBirthDeath.cpp']+glob('C/*.c'), #'../src/LNAModule.cpp', 
					runtime_library_dirs=['/usr/local/lib'])


setup(name = 'birthDeathLNA',
       version = '1.0',
       description = '',
       author = 'Justin Feigelman',
       author_email = 'justin.feigelman@helmholtz-muenchen.de',
       url = '',
       long_description='',
       ext_modules = [module1])

