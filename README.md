About
======

graphslope is a Python module dedicated to solve the Graph-Slope
optimization problem as introduced in [1]. It is ditributed under
the 3-Clause BSD license.

Important links
===================

- [Source code repo](https://github.com/svaiter/gslope\_oracle\_inequality/)
- [Issue tracker](https://github.com/svaiter/gslope\_oracle\_inequality/issues)

Dependencies
===============

The required dependencies to build the software are Python >= 3.5,
setuptools, Numpy >= 1.17 and scikit-learn >= 0.20.

For running the notebooks, you will also needs osmnx >= 0.5.1.

Install
=========

This package uses distutils, which is the default way of installing
python modules. To install it locally, use::

  pip install -e .

To install for all users on Unix/Linux::

  python setup.py build

  python setup.py install

[1]: A sharp oracle inequality for Graph-Slope. Pierre C. Bellec, Joseph Salmon, and Samuel Vaiter. Electron. J. Statist., 11(2):4851–4870, 2017.
