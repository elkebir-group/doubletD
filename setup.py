#! /usr/bin/env python
# -*- coding: utf-8 -*-
import setuptools

with open('README.md') as f:
    long_description = f.read()

setuptools.setup(
    name='doubletD',
    packages=["doubletD"],
    description="Doublet detection in scDNA-seq data",
    long_description=long_description,
    long_description_content_type='text/markdown',
    version='0.1.0',
    url='http://github.com/elkebir-group/doubletD',
    author='Leah Weber, Palash Sashittal',
    author_email='leahlw2@illinois.edu',
    python_requires='>=3.6',
    scripts=[
        'scripts/doubletD',
    ],
    install_requires=[
        "numpy",
        "pandas",
    ],
)

