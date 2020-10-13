import os
import pathlib
import shutil


def detect(predicated, collection, if_not_none=lambda element: element, default=None, if_none=lambda: None):
    result = next((if_not_none(element) for element in collection if predicated(element)), default)
    return result or if_none()


def flat_map(function, collection):
    return [y for ys in collection for y in function(ys)]


def root_path():
    return pathlib.Path(__file__).parent.parent.absolute()


def create_results_dir(file, root=root_path()):
    filename = file.split('/')[-1].split('.')[0]
    dir_path = f'{root}/results/{filename}'
    if os.path.isdir(dir_path):
        shutil.rmtree(dir_path)
    if not os.path.isdir(dir_path):
        os.mkdir(dir_path)
    return dir_path


def get_element(elements, index, default=None):
    try:
        return elements[index]
    except IndexError:
        return default
