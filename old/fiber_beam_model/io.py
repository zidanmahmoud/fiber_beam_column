import inspect
import re

def debug(variable):
    """ prints variable name and value """
    frame = inspect.currentframe()
    try:
        context = inspect.getframeinfo(frame.f_back).code_context
        caller_lines = ''.join([line.strip() for line in context])
        search = re.search(r'debug\s*\((.+?)\)$', caller_lines)
        if search:
            caller_lines = search.group(1)
        print(caller_lines, ":\n", variable, end="")
        input()
    finally:
        del frame
