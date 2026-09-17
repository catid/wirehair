"""Pure helper and pre-loader rejection tests; no native codec work."""
from pathlib import Path
import tempfile
import unittest
from unittest import mock
import Support as A


class Test(unittest.TestCase):
    def test_recursive_types(self):
        for actual,expected in ((True,1),(1.0,1),([True],[1]),
                                ({'a':[True]}, {'a':[1]}),({'a':(1.0,)},{'a':(1,)}),
                                ({True},{1}),({(1.0,):0},{(1,):0})):
            with self.subTest(actual=actual,expected=expected):
                with self.assertRaises(ValueError): A.exact(actual,expected,'typed')
        value={'a':[1,True,1.0], 'b':(None,{'c','d'})}
        A.exact(value,value,'same')
        A.exact({'a':1,'b':2},{'b':2,'a':1},'unordered')

    def test_json(self):
        value={'a':[1,True,None],'b':'text'}
        self.assertEqual(A.decode(A.canonical(value)),value)
        for raw in (b'{"a":1,"a":2}',b'{"x":{"a":1,"a":2}}',b'NaN',b'Infinity'):
            with self.assertRaises(ValueError): A.decode(raw)
        with self.assertRaises(ValueError): A.canonical(float('nan'))

    def test_files_and_publication(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'proof.json'
            A.publish(path,{'x':1})
            raw=A.canonical({'x':1})
            self.assertEqual(A.read_regular(path,len(raw)),raw)
            with self.assertRaises(ValueError): A.read_regular(path,len(raw)-1)
            with self.assertRaises(FileExistsError): A.publish(path,{'x':2})
            self.assertEqual(A.pin(path),dict(path=str(path),bytes=len(raw),sha256=A.sha(raw)))
            self.assertEqual(path.stat().st_mode & 0o777,0o400)

    def test_wrong_hash_fails_before_loader(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'not-a-library.so'; path.write_bytes(b'not a DSO')
            with mock.patch.object(A.T,'CDLL') as loader:
                with self.assertRaisesRegex(ValueError,'exact diagnostic baseline DSO'):
                    A.Library(path,'0'*64)
                loader.assert_not_called()


if __name__=='__main__': unittest.main()
