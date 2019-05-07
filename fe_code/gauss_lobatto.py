import numpy as np
from scipy.special import legendre


def GaussLobatto(Num):
    """
    note: stable up to Num=37
    http://mathworld.wolfram.com/LobattoQuadrature.html
    """
    if not isinstance(Num, int):
        raise ValueError(f"Value must be of type int. given type: {type(Num)}")
    if Num < 2 or Num > 37:
        raise ValueError(f"Num can only be between 2 and 37")

    x = list(np.sort(legendre(Num - 1).deriv().roots))
    x.insert(0, -1)
    x.append(1)
    x = np.array(x)
    w = 2 / (Num * (Num - 1) * (legendre(Num - 1)(x)) ** 2)
    return x, w
