import pathlib


def detect(predicated, collection, if_not_none=lambda element: element, default=None):
    return next((if_not_none(element) for element in collection if predicated(element)), default)


def root_path():
    return pathlib.Path(__file__).parent.parent.absolute()
