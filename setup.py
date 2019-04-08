#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pkg_resources import parse_version
from setuptools    import setup, find_packages
from os            import path

needed_modules = {'astropy': 1.3,
                  'numpy'  : 1.8,
                  'scipy'  : 0.16,
                  'pandas' : 0.22,
                  'opencv' : 3.4,
                  'NFW'    : 0.2,
                  'aplpy'  : 1.1,
                  'ephem'  : 3.7,
                  'PyPDF2' : 1.26}


            
#These are just estimated because I do not know all side-effects of these packages

def get_misc_status(neededmodule):
    """
    Returns a dictionary containing a boolean specifying whether the needed software
    is up-to-date, along with the version string (empty string if
    not installed).
    """

    module_status = {}
    try:
        test = __import__(neededmodule)
        '''
        import numpy
        numpy_version = numpy.__version__
        '''
        module_version = test.__version__
        module_status['up_to_date'] = parse_version(module_version) >= parse_version(needed_modules[neededmodule])
        module_status['version']= module_version
    except ImportError:
        module_status['up_to_date'] = False
        module_status['version'] = ""
    return module_status

def setup_clusterbuster():

    here = path.abspath(path.dirname(__file__))

    # Get the long description from the README file
    with open(path.join(here, 'README.md')) as f:  #, encoding='utf-8'
        long_description = f.read()
        
#    for module in needed_modules.keys():
#        module_status = get_misc_status(module)
#        module_req_str = "clusterbuster requires {0} >= {1}.\n".format(module,needed_modules[module])
#        
#        
#        if module_status['up_to_date'] is False:
#            if module_status['version']:
#                raise ImportError("Your installation of {0}""{1} is out-of-date.\n{2}".format(module,module_status['version'],module_req_str))
#            else:
#                raise ImportError("{0} is not installed.\n{0}".format(module,module_req_str))
                                  
    setup(
        name='clusterbuster',
        version='0.9.0',

        # Author details
        author="Jakob Gelszinnis",
        author_email="jgelszinnis@gmail.com ",

        # The project's main homepage.
        url="https://github.com/jakgel/clusterbuster",

        description='A Python implementation of a joint analysis framework for radio relic surveys and simulations',
        long_description=long_description,
        license='MIT',

        # See https://pypi.python.org/pypi?%3Aaction=list_classifiers
        classifiers=[
            # How mature is this project? Common values are
            #   3 - Alpha
            #   4 - Beta
            #   5 - Production/Stable
            'Development Status :: 4 - Beta',


            'Environment :: Console',
            'Operating System :: Linux; Windows unknown',

            # Indicate who your project is intended for
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Astronomy',

            # Pick your license as you wish (should match "license" above)
            'License :: MIT License',

            # Specify the Python versions you support here. In particular, ensure
            # that you indicate whether you support Python 2, Python 3 or both.
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
        ],

        # What does your project relate to?
        keywords='astrophysics astronomy galaxy-clusters relics nfw navarro-frenk-white',

        # List run-time dependencies here.  These will be installed by pip when
        # your project is installed. For an analysis of "install_requires" vs
        # pip's requirements files see:
        # https://packaging.python.org/en/latest/requirements.html
        install_requires=['NumPy (>=1.8)',
                          'pandas (>=0.22)',
                          'scipy (>=0.16)',
                          'opencv (>=3.4)',
                          'aplpy (>=1.1)',
                          'NFW (>=0.2)',
                          'astropy (>=1.3)',
                          'ephem (>=3.7)',
                          'PyPDF2 (>=1.26)'],

        # You can just specify the packages manually here if your project is
        # simple. Or you can use find_packages().
        #packages=find_packages(exclude=['contrib', 'docs', 'tests']),
        packages=["clusterbuster", "surveyexample", "inference", "surveysim"],
        )

        
if __name__ == '__main__':
    setup_clusterbuster()

