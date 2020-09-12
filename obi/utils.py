import pathlib


def detect(predicated, collection, if_not_none=lambda element: element, default=None, if_none=lambda: None):
    result = next((if_not_none(element) for element in collection if predicated(element)), default)
    if not result:
        result = if_none()

    return result


def root_path():
    return pathlib.Path(__file__).parent.parent.absolute()
