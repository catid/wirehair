#!/usr/bin/env python3
"""Render the recorded WH1/WH2 snapshot as a deterministic, standalone SVG."""

import argparse
import csv
import math
from pathlib import Path
import xml.etree.ElementTree as ET


DIRECTORY = Path(__file__).resolve().parent
BLOCKS = (8, 128, 512, 1024)
WIDTHS = (64, 1280)
SERIES = (("wh1_ms", "Wirehair 1", "#0072b2"),
          ("wh2_ms", "Wirehair 2 (default)", "#c45100"))


def read_data(path):
    """Fail closed if this fixed snapshot has missing/duplicate/invalid cells."""
    data = {}
    with path.open(newline="", encoding="utf-8") as source:
        reader = csv.DictReader(source)
        if reader.fieldnames != ["blocks", "block_bytes", "wh1_ms", "wh2_ms"]:
            raise ValueError("unexpected CSV header")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("unexpected CSV row length")
            key = (int(row["blocks"]), int(row["block_bytes"]))
            if key in data:
                raise ValueError("duplicate timing cell")
            values = {field: float(row[field]) for field, _, _ in SERIES}
            if any(not math.isfinite(v) or not 0.002 <= v <= 2.0
                   for v in values.values()):
                raise ValueError("timing outside fixed plot range [0.002, 2] ms")
            data[key] = values
    if set(data) != {(k, b) for k in BLOCKS for b in WIDTHS}:
        raise ValueError("expected four block counts at each of two widths")
    return data


def render(data):
    svg = ET.Element("svg", {
        "xmlns": "http://www.w3.org/2000/svg", "width": "1020", "height": "470",
        "viewBox": "0 0 1020 470", "role": "img",
        "aria-labelledby": "title description",
    })
    ET.SubElement(svg, "title", id="title").text = "Wirehair 1 vs Wirehair 2 lifecycle time"
    ET.SubElement(svg, "desc", id="description").text = (
        "Two panels compare 64-byte and 1280-byte blocks at K=8, 128, 512, 1024. "
        "Lower time is better; both axes are logarithmic. WH2 is slower at all "
        "64-byte points and faster at 1280 bytes except K=8. Rounded one-host "
        "no-loss diagnostic timings, not recovery or universal speed evidence.")

    def element(tag, **attrs):
        return ET.SubElement(svg, tag, {k.replace("_", "-"): str(v)
                                        for k, v in attrs.items()})

    def label(x, y, text, size=14, **attrs):
        node = element("text", x=x, y=y, fill="#243447", font_size=size,
                       font_family="sans-serif", **attrs)
        node.text = text

    element("rect", width=1020, height=470, fill="white")
    label(510, 31, "Wirehair 1 vs Wirehair 2", 23, text_anchor="middle")
    label(510, 55, "No-loss lifecycle · lower is better · logarithmic axes", 15,
          text_anchor="middle")

    for index, (_, name, color) in enumerate(SERIES):
        x = 310 + index * 195
        element("line", x1=x, y1=80, x2=x + 28, y2=80, stroke=color,
                stroke_width=3, stroke_dasharray="none" if index == 0 else "7 4")
        label(x + 36, 85, name)

    top, bottom = 130, 370
    for panel, width in enumerate(WIDTHS):
        left, right = 84 + panel * 510, 470 + panel * 510

        def x_coord(k):
            return left + math.log2(k / 8) / 7 * (right - left)

        def y_coord(ms):
            return bottom - math.log10(ms / 0.002) / 3 * (bottom - top)

        label((left + right) / 2, 115, "B = {:,} bytes per block".format(width),
              17, text_anchor="middle")
        for tick in (0.002, 0.01, 0.1, 1, 2):
            y = y_coord(tick)
            element("line", x1=left, y1=y, x2=right, y2=y, stroke="#dde3e9")
            label(left - 10, y + 5, "{:g}".format(tick), 12, text_anchor="end")
        for k in BLOCKS:
            x = x_coord(k)
            element("line", x1=x, y1=top, x2=x, y2=bottom, stroke="#eef1f4")
            label(x, bottom + 22, str(k), 13, text_anchor="middle")
        label((left + right) / 2, 420, "Original blocks (K)", 14, text_anchor="middle")
        label(left - 58, (top + bottom) / 2, "Lifecycle time (ms)", 14,
              text_anchor="middle",
              transform="rotate(-90 {} {})".format(left - 58, (top + bottom) / 2))

        for index, (field, name, color) in enumerate(SERIES):
            points = [(x_coord(k), y_coord(data[(k, width)][field])) for k in BLOCKS]
            element("polyline", points=" ".join("{:.3f},{:.3f}".format(x, y)
                                                for x, y in points),
                    fill="none", stroke=color, stroke_width=2.5,
                    stroke_dasharray="none" if index == 0 else "7 4")
            for k, (x, y) in zip(BLOCKS, points):
                # Shape and dash pattern distinguish the series without color.
                if index == 0:
                    marker = element("circle", cx=x, cy=y, r=4, fill=color)
                else:
                    marker = element("rect", x=x - 4, y=y - 4, width=8, height=8,
                                     fill=color)
                ET.SubElement(marker, "title").text = (
                    "{}: K={}, B={}, {:.3f} ms".format(name, k, width,
                                                      data[(k, width)][field]))
    label(510, 452, "2026-09-17 · rounded one-host diagnostic · source: wh1-vs-wh2.csv",
          13, text_anchor="middle")
    return ET.tostring(svg, encoding="utf-8", xml_declaration=True) + b"\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="verify SVG is up to date")
    args = parser.parse_args()
    content = render(read_data(DIRECTORY / "wh1-vs-wh2.csv"))
    output = DIRECTORY / "wh1-vs-wh2.svg"
    if args.check:
        if not output.is_file() or output.read_bytes() != content:
            parser.exit(1, "SVG is stale; run plot_wh1_vs_wh2.py\n")
        print("SVG matches source data")
    else:
        output.write_bytes(content)
        print("Wrote {}".format(output))


if __name__ == "__main__":
    main()
