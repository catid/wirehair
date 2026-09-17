"""Bounded, shell-free capture for the baseline-only recovery inventory."""
import os
import resource
import selectors
import signal
import subprocess
import time

WALL_SECONDS = 150
STDOUT_BYTES = 128*1024*1024
STDERR_BYTES = 64*1024
ENVIRONMENT = dict(PATH='/usr/bin:/bin', LANG='C', LC_ALL='C', TZ='UTC')


def limits():
    return dict(cpu_seconds=120, address_space_bytes=512*1024*1024,
                core_bytes=0, file_bytes=128*1024*1024,
                stdout_bytes=STDOUT_BYTES, stderr_bytes=STDERR_BYTES,
                wall_seconds=WALL_SECONDS)


def set_limits():
    cap = limits()
    for kind, value in ((resource.RLIMIT_CPU,cap['cpu_seconds']),
                        (resource.RLIMIT_AS,cap['address_space_bytes']),
                        (resource.RLIMIT_CORE,cap['core_bytes']),
                        (resource.RLIMIT_FSIZE,cap['file_bytes'])):
        resource.setrlimit(kind,(value,value))


def kill_group(process):
    # Kill descendants too, including a child holding pipes after parent exit.
    try:
        os.killpg(process.pid,signal.SIGKILL)
    except ProcessLookupError:
        pass
    process.wait(timeout=5)


def capture(command, stdout_path, stderr_path):
    """Preserve bounded prefixes on failure; never retry a worker.

    The controller is single-threaded. preexec_fn applies resource limits only
    in the child; start_new_session gives it a dedicated killable process group.
    The explicit environment excludes inherited loader/allocator overrides.
    """
    start = time.monotonic()
    counts = [0,0]
    failure = None
    caps = (STDOUT_BYTES,STDERR_BYTES)
    with stdout_path.open('xb') as out, stderr_path.open('xb') as err:
        with selectors.DefaultSelector() as selector:
            process = subprocess.Popen(command, stdin=subprocess.DEVNULL,
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                       env=ENVIRONMENT.copy(), close_fds=True,
                                       start_new_session=True, preexec_fn=set_limits)
            try:
                selector.register(process.stdout,selectors.EVENT_READ,(0,out))
                selector.register(process.stderr,selectors.EVENT_READ,(1,err))
                while selector.get_map():
                    remaining = WALL_SECONDS-(time.monotonic()-start)
                    if remaining <= 0:
                        failure = 'TIMEOUT'
                        break
                    for key,_ in selector.select(min(remaining,0.1)):
                        index,destination = key.data
                        block = os.read(key.fileobj.fileno(),65536)
                        if not block:
                            selector.unregister(key.fileobj)
                            continue
                        allowed = min(len(block),caps[index]-counts[index])
                        destination.write(block[:allowed])
                        counts[index] += allowed
                        if allowed != len(block):
                            failure = ('STDOUT_LIMIT','STDERR_LIMIT')[index]
                            break
                    if failure:
                        break
                if failure:
                    kill_group(process)
                else:
                    try:
                        process.wait(timeout=max(0,WALL_SECONDS-(time.monotonic()-start)))
                    except subprocess.TimeoutExpired:
                        failure = 'TIMEOUT'
                        kill_group(process)
                return dict(exit=process.returncode,failure=failure,
                            stdout_bytes=counts[0],stderr_bytes=counts[1],
                            wall_seconds=time.monotonic()-start)
            except BaseException:
                kill_group(process)
                raise
            finally:
                process.stdout.close()
                process.stderr.close()
