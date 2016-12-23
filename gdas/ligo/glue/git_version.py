id = "0.9.8"
date = ""
branch = ""
tag = ""
if tag == "None":
    tag = None
author = ""
builder = ""
committer = ""
status = ""
version = id
verbose_msg = """Branch: 
Tag: 
Id: 0.9.8

Builder: 
Build date: 
Repository status: """

import warnings

class VersionMismatchError(ValueError):
    pass

def check_match(foreign_id, onmismatch="raise"):
    """
    If foreign_id != id, perform an action specified by the onmismatch
    kwarg. This can be useful for validating input files.

    onmismatch actions:
      "raise": raise a VersionMismatchError, stating both versions involved
      "warn": emit a warning, stating both versions involved
    """
    if onmismatch not in ("raise", "warn"):
        raise ValueError, onmismatch + " is an unrecognized value of onmismatch"
    if foreign_id == "0.9.8":
        return
    msg = "Program id (0.9.8) does not match given id (%s)." % foreign_id
    if onmismatch == "raise":
        raise VersionMismatchError, msg

    # in the backtrace, show calling code
    warnings.warn(msg, UserWarning)

