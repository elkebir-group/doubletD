#! /usr/bin/env python
# -*- coding: utf-8 -*-
import setuptools

with open('README.md') as f:
    long_description = f.read()

setuptools.setup(
    name='doubletd',
    packages=["doubletd"],
    description="detecting doublets in single-cell DNA sequencing data",
    long_description=long_description,
    long_description_content_type='text/markdown',
    version='0.1.0',
    url='http://github.com/elkebir-group/doubletD',
    author='Palash Sashittal and Leah Weber',
    author_email='sashitt2@illinois.edu',
    python_requires='>=3.6',
    scripts=[
        'doubletd/scripts/doubletd',
    ],
    install_requires=[
        "numpy",
        "pandas",
    ],
)
