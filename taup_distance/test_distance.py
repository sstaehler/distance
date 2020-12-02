#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

:copyright:
    Simon Stähler (mail@simonstaehler.com), 2020
:license:
    None
"""
from taup_distance.taup_distance import get_dist
from obspy.taup import TauPyModel
model = TauPyModel('iasp91')

dist = get_dist(model=model,
                tSmP=400,
                depth=50.)
print(dist)
