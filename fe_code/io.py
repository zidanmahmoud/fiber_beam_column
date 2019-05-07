"""Printing functions"""
def debug(message, *args, **kwargs):
    """ normal print function """
    print(message, *args, **kwargs)


def warning(message, *args, **kwargs):
    """ print warning in yellow upper case """
    print("\33[93m" + "WARNING: " + message.upper() + "\33[0m", *args, **kwargs)
