"""Read-only exact parity with installed K8 evidence; no codec/controller imports."""
from pathlib import Path

BUNDLE = Path('/var/tmp/wh2-k8-public-recovery-r0')
QUALIFIED = Path('/tmp/wh2-k8-public-recovery-qualified.i7v7bPnq')
PINS = {
    BUNDLE/'COMPLETE.json': '699415eef17fcc486aeb3e5487feef0a9dcd4812110497df1370014c47b7af02',
    BUNDLE/'native.raw.jsonl': 'bab96ca60d697f9f5a751fdbc74dce306c6c0d12a98e651f7e9735beb6af50e1',
    QUALIFIED/'native/fixtures.jsonl': '236469d3aa985ba69e153945c53066ff24449e2691494e317683c3f040e8fccb',
    QUALIFIED/'scalar/fixtures.jsonl': 'cc35b5553809efa34a3f1b8414f5e55427ae35931517fd9622a41755f7e9e54d',
    QUALIFIED/'asan/fixtures.jsonl': 'bb1c2ceec3a9d0ece78201e4e40cfe81b6e51599a14ab5240dfef623a02d90e8',
}
PROTOCOL = 'wirehair.wh2.k8-public-recovery-r0'


def read(path, io):
    raw = io.read_regular(path, 96*1024**2)
    io.exact(io.sha(raw), PINS[path], 'installed evidence bytes')
    return raw


def records(raw, count, mode, scope, io):
    io.require(raw.endswith(b'\n'), 'complete installed stream')
    rows = [io.decode(line) for line in raw.splitlines()]
    io.exact(len(rows), count+2, 'whole installed record stream')
    io.exact((rows[0]['type'], rows[0]['protocol'], rows[0]['backend'], rows[0]['scope']),
             ('header', PROTOCOL, mode, scope), 'installed stream identity')
    io.exact(rows[-1], dict(type='footer', records=count, checked=True), 'installed complete footer')
    return rows[1:-1]


def verify(actual, io, pin):
    """No descriptor exceptions: ordinary and explicit WHV2 must agree exactly."""
    complete = io.decode(read(BUNDLE/'COMPLETE.json', io))
    io.exact((complete['protocol'], complete['outcome']), (PROTOCOL, 'PASS'), 'installed recovery identity')
    path = BUNDLE/'native.raw.jsonl'
    raw = read(path, io)
    io.exact([p for p in complete['files'] if p['path'] == str(path)], [pin(path)],
             'installed raw COMPLETE binding')
    expected = records(raw, 6260, 'native', 'retained', io)
    io.exact(actual, expected, 'exact installed descriptor/payload/rank/status/ledger parity')


def neutral(actual, mode, io):
    io.require(mode in ('native', 'scalar', 'asan'), 'installed neutral backend')
    expected = records(read(QUALIFIED/mode/'fixtures.jsonl', io), 48, mode, 'neutral', io)
    io.exact(actual, expected, 'exact installed neutral packet/rank/ledger parity')
