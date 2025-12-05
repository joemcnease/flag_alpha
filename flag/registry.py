"""
Register and dispatch fluid models.
"""
from collections.abc import Callable

from .types import TFVec
from .exceptions import UnknownModelError, UnsupportedPropertyError


class ModelRegistry:
    """ Register and dispatch models.
    """
    def __init__(self, models: dict[str, dict[str, Callable]]):
        self._models = models

    @property
    def names(self) -> list[str]:
        return list(self._models.keys())

    def dispatch(self, model: str, prop: str) -> Callable:
        if model not in self._models:
            raise UnknownModelError(f"Unknown model {model}. "
                                    f"Available models: {', '.join(self._models)}")
        if prop not in self._models[model]:
            raise UnsupportedPropertyError(f"Model {model} does not support property {prop}.")

        return self._models[model][prop]

    def call(self, model: str, prop: str, *args: TFVec) -> TFVec:
        return self.dispatch(model, prop)(*args)
