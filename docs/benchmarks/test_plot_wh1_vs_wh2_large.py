import unittest
import xml.etree.ElementTree as ET

import plot_wh1_vs_wh2_large


class LargePlotDataTest(unittest.TestCase):
    def setUp(self):
        self.speed, self.recovery = plot_wh1_vs_wh2_large.read_data(
            plot_wh1_vs_wh2_large.DIRECTORY / "wh1-vs-wh2-large-k-speed.csv",
            plot_wh1_vs_wh2_large.DIRECTORY / "wh1-vs-wh2-large-k-recovery.csv")

    def test_complete_1000_step_grid(self):
        expected = {(k, b) for k in plot_wh1_vs_wh2_large.EXPECTED_K
                    for b in plot_wh1_vs_wh2_large.BLOCKS}
        self.assertEqual(set(self.speed), expected)
        self.assertEqual(set(self.recovery), expected)
        self.assertEqual(max(plot_wh1_vs_wh2_large.LARGE_K), 64000)

    def test_svg_is_reproducible_and_well_formed(self):
        content = plot_wh1_vs_wh2_large.render(self.speed, self.recovery)
        ET.fromstring(content)
        self.assertEqual(
            content,
            (plot_wh1_vs_wh2_large.DIRECTORY / "wh1-vs-wh2-large-k.svg").read_bytes())


if __name__ == "__main__":
    unittest.main()
