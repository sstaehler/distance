#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

:copyright:
    Simon Stähler (mail@simonstaehler.com), 2020
:license:
    None
"""
import numpy as np
from matplotlib import pyplot as plt
from obspy.taup.helper_classes import TauModelError
from obspy.taup import TauPyModel
from collections import OrderedDict
from argparse import ArgumentParser
from obspy.taup.taup_create import build_taup_model
from os.path import split as psplit


def define_arguments():
    helptext = 'Compute taup_distance from TauP model, given a S-P traveltime difference'
    parser = ArgumentParser(description=helptext)

    helptext = "Input TauP file"
    parser.add_argument('fnam_nd', help=helptext)

    helptext = "S-P arrival time differences"
    parser.add_argument('times', nargs='+', type=float, help=helptext)

    helptext = "Plot convergence in T-X plot (default: False)"
    parser.add_argument('--plot', action="store_true", default=False,
                        help=helptext)
    return parser.parse_args()


def get_dist(model: TauPyModel,
             tSmP: float,
             depth: float,
             phase_list=('P', 'S'),
             plot=False):
    """
    Get taup_distance of an event, given difference between two phase arrival times (default: P and S)
    :param model: TauPyModel object for given velocity model
    :param tSmP: time difference between arrivals in seconds
    :param depth: depth of event in km
    :param phase_list: list-like object with two phase names
    :param plot: plot convergence (default: False)
    :return: taup_distance of event in degree
    """
    from scipy.optimize import newton

    # Reasonable guess 1
    dist0 = tSmP / 6.5
    try:
        dist = newton(func=_get_TSmP, fprime=_get_SSmP,
                      x0=dist0, args=(model, tSmP, phase_list, plot, depth),
                      maxiter=10)
    except RuntimeError:
        dist = None
    if dist is None:
        # Reasonable guess 2
        dist0 = tSmP / 8.
        try:
            dist = newton(func=_get_TSmP, fprime=_get_SSmP,
                          x0=dist0,
                          args=(model, tSmP, phase_list, plot, depth),
                          maxiter=10)
        except RuntimeError:
            dist = None

    if plot:
        plt.axhline(tSmP)
        plt.show()
    return dist


def _get_TSmP(distance: float,
              model: TauPyModel,
              tmeas: float,
              phase_list: tuple,
              plot: bool,
              depth: float):
    """
    Compute travel time difference between two phases (minus tmeas)
    :param distance: taup_distance in degree
    :param model: TauPy model
    :param tmeas: measured travel time difference
    :param phase_list: phase names
    :param plot: plot convergence
    :param depth: depth of event in km
    :return: travel time difference between phases minus tmeas
    """
    if len(phase_list) != 2:
        raise ValueError('Only two phases allowed')
    tP = None
    tS = None
    try:
        arrivals = model.get_travel_times(source_depth_in_km=depth,
                                          distance_in_degree=distance,
                                          phase_list=phase_list)
    except (ValueError, TauModelError):
        pass
    else:
        for arr in arrivals:
            if arr.name == phase_list[0] and tP is None:
                tP = arr.time
            elif arr.name == phase_list[1] and tS is None:
                tS = arr.time
    if tP is None or tS is None:
        if plot:
            plt.plot(distance, -1000, 'o')
        return -1000.
    else:
        if plot:
            plt.plot(distance, tS - tP, 'o')
        return (tS - tP) - tmeas


def _get_SSmP(distance: float,
              model: TauPyModel,
              tmeas: float,
              phase_list: tuple,
              plot: bool,
              depth: float):
    """
    Compute derivative of travel time difference between two phases (minus tmeas)
    Uses slowness of the two phases.
    Used to compute derivative for fitting. Some unused arguments, but argument list needs to be identical to _get_TSmP.
    :param distance: taup_distance in degree
    :param model: TauPy model
    :param tmeas: measured travel time difference
    :param phase_list: phase names
    :param plot: plot convergence
    :param depth: depth of event in km
    :return: travel time difference between phases minus tmeas
    """

    if len(phase_list) != 2:
        raise ValueError('Only two phases allowed')
    sP = None
    sS = None
    try:
        arrivals = model.get_travel_times(source_depth_in_km=depth,
                                          distance_in_degree=distance,
                                          phase_list=phase_list)
    except (ValueError, TauModelError):
        pass
    else:
        for arr in arrivals:
            if arr.name == phase_list[0] and sP is None:
                sP = np.deg2rad(arr.ray_param)
            elif arr.name == phase_list[1] and sS is None:
                sS = np.deg2rad(arr.ray_param)

    if sP is None or sS is None:
        return -10000.
    else:
        return sS - sP


def main(fnam_nd, times, phase_list=('P', 'S'), depth=40., plot=False):
    fnam_npz = './taup_tmp/' \
               + psplit(fnam_nd)[-1][:-3] + '.npz'
    build_taup_model(fnam_nd,
                     output_folder='./taup_tmp'
                     )
    cache = OrderedDict()
    model = TauPyModel(model=fnam_npz, cache=cache)

    for tSmP in times:
        dist = get_dist(model, tSmP=tSmP, depth=depth, phase_list=phase_list, plot=plot)
        if dist is None:
            print(f'{fnam_nd}, S-P time: {tSmP:5.1f}: NO SOLUTION FOUND!')
        else:
            print(f'{fnam_nd}, S-P time: {tSmP:5.1f}, taup_distance: {dist:5.1f}')


if __name__ == '__main__':
    args = define_arguments()
    main(fnam_nd=args.fnam_nd, times=args.times, plot=args.plot)
