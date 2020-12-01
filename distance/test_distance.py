#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

:copyright:
    Simon Stähler (mail@simonstaehler.com), 2020
:license:
    None
"""
from distance.distance import get_distance
from obspy.taup import TauPyModel
model = TauPyModel('iasp91')

dist = get_distance(tSmP=400,
                    depth=50.)
print(dist)
