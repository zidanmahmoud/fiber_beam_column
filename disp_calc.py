# from math import abs
import numpy as np

def d(res):
    return res * 10000 / 3


def res(d):
    return d * 3 / 10000


def calculate_loadsteps(step_size):
    l = list()
    l.append(d(0.0012) / step_size)
    l.append(d(0.0012 + 0.0010) / step_size + l[-1])
    l.append(d(0.0010 + 0.0020) / step_size + l[-1])
    l.append(d(0.0020 + 0.0016) / step_size + l[-1])
    l.append(d(0.0016 + 0.000175) / step_size + l[-1])
    l.append(d(0.000175 + 0.001) / step_size + l[-1])
    l.append(d(0.001 - 0.0005) / step_size + l[-1])
    return [round(number) for number in l]

def calculate_loadsteps2(step_size):
    disps = [0, 1, -1.4, 1.5, -1.7, 2.5, -2.1, 2.4, -2.2, 2.1, -2.6, 2, -2.4, 2.3]
    l = np.cumsum(np.abs(np.diff(disps)/step_size))
    return [int(round(number)+1) for number in l]
