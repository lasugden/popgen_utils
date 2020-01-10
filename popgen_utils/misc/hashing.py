import hashlib


def hash(obj, use_new=True):
    """Hash any Python object."""
    out = hashlib.md5(repr(obj).encode())
    return out.hexdigest()
