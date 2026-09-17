"""The checkpoint must not execute codecs or consume an experiment namespace."""
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import Control as C
import Inventory as I


class TestDraft(unittest.TestCase):
    def test_all_controller_entries_fail_before_side_effects(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            with mock.patch.object(C, 'OUTPUT', root/'science'), \
                 mock.patch.object(C, 'QUALIFICATION', root/'neutral'), \
                 mock.patch.object(C, 'codec_pins') as pins, \
                 mock.patch.object(C.A, 'read_regular') as read, \
                 mock.patch.object(C.Launch, 'capture') as capture, \
                 mock.patch.object(I, 'run_inventory') as inventory:
                calls = [(C.neutral, ('native',)), (C.neutral, ('portable',)),
                         (C.inputs, ()), (C.worker, ('a'*64,)),
                         (C.run, ()), (C.replay, ())]
                for function, arguments in calls:
                    with self.subTest(entry=function.__name__, arguments=arguments):
                        with self.assertRaisesRegex(ValueError, 'unqualified draft: launch disabled'):
                            function(*arguments)
                for operation in (pins, read, capture, inventory):
                    operation.assert_not_called()
                self.assertEqual(list(root.iterdir()), [])

    def test_native_entry_fails_before_loading(self):
        emit = mock.Mock()
        with mock.patch.object(I.A, 'Library') as library:
            with self.assertRaisesRegex(ValueError, 'unqualified draft: launch disabled'):
                I.run_inventory('native', [], 'unused', emit)
            library.assert_not_called()
        emit.assert_not_called()


if __name__ == '__main__':
    unittest.main()
