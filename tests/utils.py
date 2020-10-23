import pathlib


def get_resource(file_path):
    return f"{str(pathlib.Path(__file__).parent.absolute())}/resources/{file_path}"


def results_path(file_path):
    return f"{results_dir()}/{file_path}"


def results_dir():
    return f"{str(pathlib.Path(__file__).parent.absolute())}/tests_result_files"
