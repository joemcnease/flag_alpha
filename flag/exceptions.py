"""
Exceptions for flag functions.
"""


class UnknownModelError(ValueError):
    pass


class UnsupportedPropertyError(ValueError):
    pass


def not_implemented():
    raise NotImplementedError()
