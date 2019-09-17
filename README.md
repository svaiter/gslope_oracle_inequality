About
======

graphslope is a Python module dedicated to solve the Graph-Slope
optimization problem as introduced in [1]. It is ditributed under
the 3-Clause BSD license.

WARNING: This is a python2 project and is deprecated as python2 itself.

Important links
===================

- [Source code repo](https://github.com/svaiter/gslope\_oracle\_inequality/)
- [Issue tracker](https://github.com/svaiter/gslope\_oracle\_inequality/issues)

Dependencies
===============

The required dependencies to build the software are Python >= 2.6,
setuptools, Numpy >= 1.3 and scikit-learn >= 0.18.

For running the notebooks, you will also needs osmnx >= 0.5.1.

Install
=========

This package uses distutils, which is the default way of installing
python modules. To install it locally, use::

  pip install -e .

To install for all users on Unix/Linux::

  python setup.py build
  sudo python setup.py install
