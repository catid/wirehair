import unittest
import xml.etree.ElementTree as ET

import plot_wh1_vs_wh2


class PlotDataTest(unittest.TestCase):
    def test_snapshot_cells_and_svg(self):
        data = plot_wh1_vs_wh2.read_data(
            plot_wh1_vs_wh2.DIRECTORY / "wh1-vs-wh2.csv")
        self.assertEqual(len(data), 8)
        self.assertEqual(
            set(data), {(k, b) for k in plot_wh1_vs_wh2.BLOCKS
                        for b in plot_wh1_vs_wh2.WIDTHS})
        ET.fromstring(plot_wh1_vs_wh2.render(data))

    def test_checked_in_svg_is_reproducible(self):
        data = plot_wh1_vs_wh2.read_data(
            plot_wh1_vs_wh2.DIRECTORY / "wh1-vs-wh2.csv")
        self.assertEqual(
            plot_wh1_vs_wh2.render(data),
            (plot_wh1_vs_wh2.DIRECTORY / "wh1-vs-wh2.svg").read_bytes())


if __name__ == "__main__":
    unittest.main()
