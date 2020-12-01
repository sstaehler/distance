#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Python tool to predict distance, given a Taup file and a travel time

:copyright:
    Simon Stähler (mail@simonstaehler.com), 2020
:license:
    None
'''
from setuptools import setup, find_packages

setup(
      name='distance',
      version='0.1',
      description='Python tool to predict distance given travel times.',
      url='github.com/sstaehler/distance',
      author='Simon Stähler',
      author_email='staehler@erdw.ethz.ch',
      license='None',
      packages=find_packages())
