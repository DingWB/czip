# -*- coding: utf-8 -*-
"""
@author: DingWB
"""
try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup, find_packages
from pathlib import Path
# from Cython.Build import cythonize
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()
# __version__ = "0.4.1"

setup(
    name="czip",
	setup_requires=['setuptools_scm'],
	use_scm_version=True,  # version=__version__,
    description="czip: Chunk based ZIP",
    long_description=long_description,
    long_description_content_type='text/markdown',
    author="Wubin Ding",
    author_email="ding.wu.bin.gm@gmail.com",
    url="https://github.com/DingWB/czip",
    packages=find_packages(exclude=('docs',)),
    install_requires=['pandas', 'fire', 'numpy', 'cython', 'fast-fisher'],
    include_package_data=True,
    package_data={
        '': ['*.txt', '*.tsv', '*.csv', '*.fa', '*Snakefile', '*ipynb']
    },
    entry_points={
        'console_scripts':
            [
                'czip=czip:main',
            ],
    },

    # ext_modules=cythonize("czip/cz.pyx",language_level = "3"), #python setup.py build_ext --inplace
)