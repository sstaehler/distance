#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Python tool to predict quake distance, given a Taup file and a travel time
Uses TauP as implemented for ObsPy.
:copyright:
    Simon Stähler (mail@simonstaehler.com), 2020
:license:
    None
'''
from setuptools import setup, find_packages

setup(
      name='taup_distance',
      version='0.1',
      description='Python tool to predict taup_distance given travel times.',
      url='github.com/sstaehler/taup_distance',
      author='Simon C. Staehler',
      author_email='staehler@erdw.ethz.ch',
      license='None',
      packages=find_packages())
