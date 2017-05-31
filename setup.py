import sys

if not sys.version_info[0] == 2:
    print("Invalid version of python, use python 2.")
    sys.exit(1)

#import ez_setup
#ez_setup.use_setuptools()

from setuptools import setup, Extension

setup(name="tram", version="0.1",
        install_requires=['intervaltree','pysam'],
        scripts = ["tram.py"],
        ext_modules=[],
        entry_points = {
        'console_scripts': [
            'tram = tram:main'
            ]
         }
        )

