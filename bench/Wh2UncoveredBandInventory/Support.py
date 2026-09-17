"""Small Linux-only diagnostic helpers; no historical receipt rebinding."""
import ctypes as T
import hashlib
import json
import os
from pathlib import Path
import stat

V, U32, U64, INT = T.c_void_p, T.c_uint32, T.c_uint64, T.c_int


def require(ok, reason):
    if not ok:
        raise ValueError(reason)


def exact(actual, expected, reason):
    require(type(actual) is type(expected) and actual == expected, reason)


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def canonical(value):
    return (json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode()


def decode(raw):
    def unique(pairs):
        result = {}
        for key, value in pairs:
            require(key not in result, 'duplicate JSON key')
            result[key] = value
        return result
    return json.loads(raw, object_pairs_hook=unique,
                      parse_constant=lambda _: require(False, 'nonfinite JSON'))


def read_regular(path, cap):
    with Path(path).open('rb') as stream:
        info = os.fstat(stream.fileno())
        require(stat.S_ISREG(info.st_mode) and info.st_size <= cap, 'bounded regular file')
        raw = stream.read(cap+1)
        require(len(raw) == info.st_size, 'complete bounded read')
        return raw


def pin(path):
    path = Path(path).resolve(strict=True)
    raw = read_regular(path, 256*1024**2)
    return dict(path=str(path), bytes=len(raw), sha256=sha(raw))


def publish(path, value):
    path = Path(path)
    with path.open('xb') as stream:
        stream.write(canonical(value))
    path.chmod(0o400)


class DlInfo(T.Structure):
    _fields_ = [('filename', T.c_char_p), ('base', V), ('symbol', T.c_char_p), ('address', V)]


class Library:
    """Check exact DSO bytes and every public function's actual loader owner.

    This does not claim full transitive linker/runtime input qualification.
    """
    def __init__(self, path, expected):
        self.path = Path(path).resolve(strict=True)
        exact(sha(read_regular(self.path, 16*1024**2)), expected, 'exact diagnostic baseline DSO')
        self.expected = expected
        self.loader = T.CDLL(None)
        self.dladdr = self.loader.dladdr
        self.dladdr.argtypes, self.dladdr.restype = (V, T.POINTER(DlInfo)), INT
        self.library = T.CDLL(str(self.path), mode=os.RTLD_NOW | os.RTLD_LOCAL)
        self.base = None
        self.functions = {}
        exact(self.call('wirehair_init_', INT, INT)(2), 0, 'public codec initialization')

    def location(self, function):
        address = T.cast(function, V).value
        info = DlInfo()
        require(address and self.dladdr(address, T.byref(info)) and info.filename and info.base,
                'resolved public function')
        exact(Path(os.fsdecode(info.filename)).resolve(strict=True), self.path, 'public function DSO owner')
        if self.base is None:
            self.base = info.base
        exact(info.base, self.base, 'single actual DSO base')
        require(address > self.base, 'public function offset')
        return dict(address=address, offset=address-self.base)

    def call(self, name, result, *arguments):
        require(name not in self.functions, 'unique public API binding')
        function = getattr(self.library, name)
        function.restype, function.argtypes = result, arguments
        self.functions[name] = (function, self.location(function))
        return function

    @property
    def report(self):
        return dict(path=str(self.path), sha256=self.expected, base=self.base,
                    exports={name: record for name, (_, record) in self.functions.items()})

    def check_bindings(self):
        for function, record in self.functions.values():
            exact(self.location(function), record, 'unchanged public binding')
        exact(sha(read_regular(self.path, 16*1024**2)), self.expected, 'unchanged loaded DSO')
