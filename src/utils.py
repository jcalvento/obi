def detect(predicated, collection, if_not_none=lambda element: element, default=None):
    return next((if_not_none(element) for element in collection if predicated(element)), default)
