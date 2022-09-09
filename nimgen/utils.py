# Authors: Federico Raimondo <f.raimondo@fz-juelich.de>
#          Sami Hamdan <s.hamdan@fz-juelich.de>
#          Vera Komeyer <v.komeyer@fz-juelich.de>
# License: AGPL
import logging
import json

logger = logging.getLogger(name='nimgen')
logger.setLevel(logging.INFO)
# console = logging.StreamHandler()
# logger.addHandler(console)
logger.propagate = False


def save_as_json(data, file):
    with open('data.txt', 'w') as outfile:
        json.dump(data, file, sort_keys=True, indent=4,
                  ensure_ascii=False)


class dotdict(dict):
    """dot.notation access to dictionary attributes"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__
