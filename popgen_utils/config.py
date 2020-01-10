"""Configurable defaults for flow."""
from __future__ import print_function
from builtins import input
import copy
import json
import os
import shutil

CONFIG_FILE = 'popgen_utils.cfg'
DEFAULT_FILE = 'popgen_utils.cfg.default'
CONFIG_PATHS = [
    os.path.expanduser('~/.config/popgen_utils'),
    os.path.join(os.path.dirname(__file__)),
    os.environ.get("POPGEN_UTILS_CONF"),
]

_params = None


def params(reload_=False):
    """Return a copy of the parameters dictionary.
    This is the primary function that should be used to access user-specific
    parameters.
    For 'defaults' and 'colors', they are initialized with the default values
    in this file, but overwritten by any settings in the user's config file.
    Parameters
    ----------
    reload_ : bool
        If True, reload the config file from disk.
    """
    global _params
    if reload_ or _params is None:
        _params = _load_config()
    return copy.deepcopy(_params)


def session_parameter(key, val):
    """
    Overwrite a parameter only for the current session.
    Parameters
    ----------
    key : str
        Name of a parameter
    val : variable
        New value of the parameter
    """

    global _params
    if _params is None:
        _params = _load_config()

    _params[key] = val


def reconfigure():
    """Re-set user-configurable parameters."""
    config_path = _find_config()

    print("Reconfiguring popgen_utils: {}".format(config_path))
    with open(config_path, 'r') as f:
        config = json.load(f)

    print("PATHS")
    if 'paths' not in config:
        config['paths'] = {}
    data_path = input(
        'Enter path to data: [{}] '.format(config['paths'].get('data', '')))
    if len(data_path):
        config['paths']['data'] = os.path.normpath(data_path)

    output_path = input(
        'Enter path to analyzed output files: [{}] '.format(
            config['paths'].get('output', '')))
    if len(output_path):
        config['paths']['output'] = os.path.normpath(output_path)

    with open(config_path, 'w') as f:
        json.dump(config, f, sort_keys=True, indent=4, separators=(',', ': '))

    params(reload_=True)


def default():
    """Return default parameters."""
    p = params()
    return p['defaults']


def _load_config():
    config_path = _find_config()
    with open(config_path, 'r') as f:
        loaded_config = json.load(f)
        config = {}
        for key in loaded_config:
            # Just add keys in the config file other than 'defaults'
            if key not in config:
                config[key] = loaded_config[key]
            else:
                # This assumes that these keys will also contain dicts,
                # they should.
                config[key].update(loaded_config[key])
    return config


def _find_config():
    for path in CONFIG_PATHS:
        if path is None:
            continue
        if os.path.isfile(os.path.join(path, CONFIG_FILE)):
            return os.path.join(path, CONFIG_FILE)
    config_path = _initialize_config()
    return config_path


def _initialize_config():
    for path in CONFIG_PATHS:
        if path is None:
            continue
        if not os.path.isdir(path):
            try:
                os.makedirs(path)
            except OSError:
                continue
        f = os.path.join(os.path.dirname(__file__), DEFAULT_FILE)
        try:
            shutil.copy(f, os.path.join(path, CONFIG_FILE))
        except IOError:
            continue
        else:
            config_path = os.path.join(path, CONFIG_FILE)
            print("Configuration initialized to: {}".format(config_path))
            print("Run `import flow.config as cfg; cfg.reconfigure()` " +
                  "to update.")
            return config_path
    print("Unable to find writable location.")
    return DEFAULT_FILE


if __name__ == '__main__':
    params()
# from pudb import set_trace; set_trace()