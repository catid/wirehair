"""Synthetic subprocess limits, environment, capture and lifetime tests."""
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock
import Launch as L


class Test(unittest.TestCase):
    def launch(self, source):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            result = L.capture([sys.executable,'-I','-c',source],root/'out',root/'err')
            out,err = (root/p for p in ('out','err'))
            self.assertEqual(result['stdout_bytes'],out.stat().st_size)
            self.assertEqual(result['stderr_bytes'],err.stat().st_size)
            return result,out.read_bytes(),err.read_bytes()

    def test_limits_and_environment(self):
        result,out,err = self.launch('import os,resource,json; '
            'print(json.dumps([dict(os.environ),[resource.getrlimit(k) for k in '
            '(resource.RLIMIT_CPU,resource.RLIMIT_AS,resource.RLIMIT_CORE,resource.RLIMIT_FSIZE)]]))')
        self.assertEqual((result['exit'],result['failure'],err),(0,None,b''))
        environment,limits = json.loads(out)
        environment.pop('LC_CTYPE',None)
        self.assertEqual(environment,L.ENVIRONMENT)
        self.assertEqual(limits,[[120,120],[512*1024*1024]*2,[0,0],[128*1024*1024]*2])
        self.assertLess(result['wall_seconds'],L.WALL_SECONDS)

    def test_env_not_inherited_and_exit(self):
        with mock.patch.dict(os.environ,{'LD_PRELOAD':'not-a-library','MALLOC_PERTURB_':'1'}):
            result,out,err = self.launch('import os,sys; print(os.getenv("LD_PRELOAD")); sys.exit(7)')
        self.assertEqual((result['exit'],result['failure'],out,err),(7,None,b'None\n',b''))

    def test_exact_output_bound(self):
        with mock.patch.object(L,'STDOUT_BYTES',32):
            result,out,err = self.launch('import os; os.write(1,b"x"*32)')
        self.assertEqual((result['exit'],result['failure'],out,err),(0,None,b'x'*32,b''))

    def test_stdout_limit(self):
        with mock.patch.object(L,'STDOUT_BYTES',32):
            result,out,_ = self.launch('import os,time; os.write(1,b"x"*33); time.sleep(10)')
        self.assertEqual((result['failure'],out),('STDOUT_LIMIT',b'x'*32))
        self.assertNotEqual(result['exit'],0)

    def test_stderr_limit(self):
        with mock.patch.object(L,'STDERR_BYTES',32):
            result,_,err = self.launch('import os,time; os.write(2,b"e"*33); time.sleep(10)')
        self.assertEqual((result['failure'],err),('STDERR_LIMIT',b'e'*32))
        self.assertNotEqual(result['exit'],0)

    def test_timeout_with_open_and_closed_pipes(self):
        for prelude in ('','os.close(1); os.close(2); '):
            with mock.patch.object(L,'WALL_SECONDS',0.3):
                result,_,_ = self.launch('import os,time; '+prelude+'time.sleep(10)')
            self.assertEqual(result['failure'],'TIMEOUT')
            self.assertNotEqual(result['exit'],0)
            self.assertLess(result['wall_seconds'],3)

    def test_descendant_holding_pipes(self):
        with mock.patch.object(L,'WALL_SECONDS',0.3):
            result,_,_ = self.launch('import os,time\nif os.fork(): os._exit(0)\ntime.sleep(10)')
        self.assertEqual((result['exit'],result['failure']),(0,'TIMEOUT'))
        self.assertLess(result['wall_seconds'],3)

    def test_no_overwrite_or_second_launch(self):
        with tempfile.TemporaryDirectory() as temporary:
            out,err = (Path(temporary)/p for p in ('out','err'))
            out.write_bytes(b'old')
            with mock.patch.object(L.subprocess,'Popen') as popen:
                with self.assertRaises(FileExistsError): L.capture(['ignored'],out,err)
                popen.assert_not_called()
            self.assertEqual(out.read_bytes(),b'old')


if __name__=='__main__': unittest.main()
