from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import Screen as S


class Test(unittest.TestCase):
    def fixture(self, root):
        (root/'candidate').mkdir()
        (root/'CMakeCache.txt').write_text('CMAKE_HOME_DIRECTORY:INTERNAL='+str(S.HERE)+'\n')
        (root/'candidate'/'Screen.cpp').write_text('Fixture matched = fs[0];\nWork(apis[arm], matched,\nfs[a].shape,matched.source.data(),p\n')

    def test_fixed_new_namespace_and_protocol(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            self.fixture(root)
            with patch.object(S.S, 'run') as runner:
                S.run(root)
                runner.assert_called_once_with(root, output=S.OUTPUT, here=S.HERE,
                    protocol='wh2-small-payload2-outlined-screen-r0')
            self.assertNotEqual(S.OUTPUT, S.S.OUTPUT)

    def test_wrong_or_prefix_build_rejected(self):
        for cache in ('other', 'CMAKE_HOME_DIRECTORY:INTERNAL='+str(S.HERE)+'-wrong\n'):
            with tempfile.TemporaryDirectory() as name:
                root = Path(name)
                self.fixture(root)
                (root/'CMakeCache.txt').write_text(cache)
                with patch.object(S.S, 'run') as runner:
                    with self.assertRaisesRegex(ValueError, 'outlined experiment build'):
                        S.run(root)
                    runner.assert_not_called()

    def test_old_buffer_mapping_rejected(self):
        with tempfile.TemporaryDirectory() as name:
            root = Path(name)
            self.fixture(root)
            (root/'candidate'/'Screen.cpp').write_text('output[side]')
            with patch.object(S.S, 'run') as runner:
                with self.assertRaisesRegex(ValueError, 'matched workspace worker'):
                    S.run(root)
                runner.assert_not_called()


if __name__ == '__main__':
    unittest.main()
