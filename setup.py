# -*- coding: utf-8 -*-
"""
@author: DingWB
"""
try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup, find_packages
# from Cython.Build import cythonize
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()
__version__ = "0.2"

setup(
    name="czip",
    version=__version__,
    description="czip: Chunk based ZIP",
    long_description=long_description,
    long_description_content_type='text/markdown',
    author="Wubin Ding",
    author_email="ding.wu.bin.gm@gmail.com",
    url="https://github.com/DingWB/czip",
    packages=find_packages(exclude=('docs',)),
    # install_requires=['pandas', 'fire', 'numpy'],
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

    # ext_modules=cythonize("czip/utils.py"), #python setup.py build_ext --inplace
)