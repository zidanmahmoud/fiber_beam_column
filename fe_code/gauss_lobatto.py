"""
Gaus-Lobatto rule
"""

import numpy as np
from scipy.special import legendre


def gauss_lobatto(num):
    """
    note: stable up to num=37
    http://mathworld.wolfram.com/LobattoQuadrature.html
    """
    if not isinstance(num, int):
        raise ValueError(f"Value must be of type int. given type: {type(num)}")
    if num < 2 or num > 37:
        raise ValueError(f"num can only be between 2 and 37")

    x = list(np.sort(legendre(num - 1).deriv().roots))
    x.insert(0, -1)
    x.append(1)
    x = np.array(x)
    w = 2 / (num * (num - 1) * (legendre(num - 1)(x)) ** 2)
    return x, w
