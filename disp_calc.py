def d(res):
    return res * 10000 / 3


def res(d):
    return d * 3 / 10000


def calculate_loadsteps(step_size):
    l = list()
    l.append(d(0.0012) / step_size)
    l.append(d(0.0012 + 0.0010) / step_size + l[-1])
    l.append(d(0.0010 + 0.00185) / step_size + l[-1])
    l.append(d(0.00185 + 0.0016) / step_size + l[-1])
    l.append(d(0.0016 + 0.000175) / step_size + l[-1])
    l.append(d(0.000175 + 0.0011) / step_size + l[-1])
    l.append(d(0.0011 - 0.0005) / step_size + l[-1])
    return [round(number) for number in l]
