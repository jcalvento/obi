import logging


def logger(base_path='.'):
    obi_logger = logging.getLogger('obi')
    if not obi_logger.handlers:
        obi_logger.setLevel(logging.INFO)
        fh = logging.FileHandler(f'{base_path}/logs.log')
        fh.setLevel(logging.INFO)
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        fh.setFormatter(formatter)
        ch.setFormatter(formatter)
        obi_logger.addHandler(fh)
        obi_logger.addHandler(ch)

    return obi_logger


def info(message, base_path='.'):
    logger(base_path).info(message)
