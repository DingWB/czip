# -*- coding: utf-8 -*-
"""
@author: DingWB
"""
try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup, find_packages
from Cython.Build import cythonize
from pathlib import Path
__version__=0.1
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
	name="bmzip",
	version=__version__,
	description="bmzip: A Python package for single cell DNA methylation storage",
	long_description=long_description,
	long_description_content_type='text/markdown',
	author="Wubin Ding",
	author_email="ding.wu.bin.gm@gmail.com",
	url="https://github.com/DingWB/bmzip",
	packages=find_packages(exclude=('doc',)),
	install_requires=['pandas','fire','numpy'],
	include_package_data=True,
	package_data={
		'': ['*.txt', '*.tsv', '*.csv', '*.fa', '*Snakefile', '*ipynb','*yaml']
	},
	entry_points={
		'console_scripts':
			[
				'bmzip=bmzip:main',
			],
		},

	# ext_modules=cythonize("bmzip/utils.py"), #python setup.py build_ext --inplace
)