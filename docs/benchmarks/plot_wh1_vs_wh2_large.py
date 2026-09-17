#!/usr/bin/env python3
"""Render the large-K WH1/WH2 speed and recovery snapshots as SVG."""

import argparse
import csv
import math
from pathlib import Path
import xml.etree.ElementTree as ET


DIRECTORY = Path(__file__).resolve().parent
SMALL_K = (8, 128, 512, 1024)
LARGE_K = tuple(range(1000, 64001, 1000))
BLOCKS = (64, 1280)
EXPECTED_K = SMALL_K + LARGE_K
COLORS = ("#0072b2", "#c45100")


def _read(path, fields, trials_field, expected_trials):
    result = {}
    with path.open(newline="", encoding="utf-8") as source:
        reader = csv.DictReader(source)
        if reader.fieldnames != fields:
            raise ValueError("unexpected CSV header in {}".format(path.name))
        for row in reader:
            key = (int(row["blocks"]), int(row["block_bytes"]))
            if key in result:
                raise ValueError("duplicate cell {}".format(key))
            if int(row[trials_field]) != expected_trials:
                raise ValueError("unexpected trial count in {}".format(path.name))
            values = {}
            loss_rate = float(row["loss_rate"])
            if not math.isfinite(loss_rate):
                raise ValueError("non-finite loss_rate in {}".format(path.name))
            values["loss_rate"] = loss_rate
            for field in fields[4:]:
                value = float(row[field])
                if not math.isfinite(value):
                    raise ValueError("non-finite {}".format(field))
                values[field] = value
            result[key] = values
    expected = {(k, b) for k in EXPECTED_K for b in BLOCKS}
    if set(result) != expected:
        raise ValueError("expected {} K/B cells in {}".format(
            len(expected), path.name))
    return result


def read_data(speed_path, recovery_path):
    speed_fields = [
        "blocks", "block_bytes", "speed_trials", "loss_rate",
        "wh1_lifecycle_ms", "wh2_lifecycle_ms"]
    recovery_fields = [
        "blocks", "block_bytes", "trials", "loss_rate",
        "wh1_exact_k_rate", "wh2_exact_k_rate", "wh1_eventual_rate",
        "wh2_eventual_rate", "wh1_terminal_failures",
        "wh2_terminal_failures"]
    speed = _read(speed_path, speed_fields, "speed_trials", 5)
    recovery = _read(recovery_path, recovery_fields, "trials", 16)
    for key in speed:
        if speed[key]["loss_rate"] != 0.0:
            raise ValueError("timing data must be lossless in {}".format(key))
        if recovery[key]["loss_rate"] != 0.10:
            raise ValueError("recovery data must use 10% loss in {}".format(key))
        if not 0.0 <= recovery[key]["wh1_exact_k_rate"] <= 1.0 or \
                not 0.0 <= recovery[key]["wh2_exact_k_rate"] <= 1.0 or \
                not 0.0 <= recovery[key]["wh1_eventual_rate"] <= 1.0 or \
                not 0.0 <= recovery[key]["wh2_eventual_rate"] <= 1.0:
            raise ValueError("recovery rate outside [0, 1]")
        if recovery[key]["wh1_eventual_rate"] < recovery[key]["wh1_exact_k_rate"] or \
                recovery[key]["wh2_eventual_rate"] < recovery[key]["wh2_exact_k_rate"]:
            raise ValueError("eventual rate below exact-K rate in {}".format(key))
    return speed, recovery


def render(speed, recovery):
    width, height = 1260, 730
    svg = ET.Element("svg", {
        "xmlns": "http://www.w3.org/2000/svg", "width": str(width),
        "height": str(height), "viewBox": "0 0 {} {}".format(width, height),
        "role": "img", "aria-labelledby": "title description",
    })
    ET.SubElement(svg, "title", id="title").text = (
        "Wirehair 1 versus Wirehair 2 through K=64000")
    ET.SubElement(svg, "desc", id="description").text = (
        "Four panels show lifecycle milliseconds and exact-K and eventual recovery rates "
        "for 64-byte and 1280-byte blocks. K is sampled every 1000 from "
        "1000 through 64000, with the original small K points included. "
        "Recovery uses 16 paired IID loss trials at ten percent loss and a K+4 horizon. "
        "Lower lifecycle time is better; higher exact-K rate is better.")

    def element(tag, **attrs):
        return ET.SubElement(svg, tag, {k.replace("_", "-"): str(v)
                                        for k, v in attrs.items()})

    def label(x, y, text, size=13, **attrs):
        node = element("text", x=x, y=y, fill="#243447", font_size=size,
                       font_family="sans-serif", **attrs)
        node.text = text

    def line(x1, y1, x2, y2, stroke="#dce3e9", width=1, dash=None):
        element("line", x1=x1, y1=y1, x2=x2, y2=y2, stroke=stroke,
                stroke_width=width, stroke_dasharray=dash or "none")

    element("rect", width=width, height=height, fill="white")
    label(width / 2, 31, "Wirehair 1 vs Wirehair 2: large-K sweep", 24,
          text_anchor="middle")
    label(width / 2, 55,
          "K = 1,000, 2,000, …, 64,000 plus K = 8, 128, 512, 1,024",
          15, text_anchor="middle")
    label(width / 2, 78,
          "Top: no-loss lifecycle time · bottom: exact-K and eventual recovery at 10% IID loss",
          14, text_anchor="middle")

    # Legend: shape and dash pattern distinguish codecs without relying on color.
    for index, name in enumerate(("Wirehair 1", "Wirehair 2 (default)")):
        x = 495 + index * 190
        line(x, 99, x + 28, 99, COLORS[index], 3,
             None if index == 0 else "7 4")
        label(x + 38, 104, name)
    line(875, 99, 903, 99, "#6b7280", 2, "2 4")
    label(913, 104, "fine dotted = eventual")

    top, bottom = 130, 675
    panel_w, panel_h = 525, 235
    panels = ((72, top, "time", 64), (663, top, "time", 1280),
              (72, 405, "recovery", 64), (663, 405, "recovery", 1280))
    x_ticks = (0, 10000, 20000, 30000, 40000, 50000, 60000, 64000)

    for left, panel_top, kind, block_bytes in panels:
        right = left + panel_w
        panel_bottom = panel_top + panel_h
        title = ("Lifecycle time (ms), B={:,}" if kind == "time" else
                 "Recovery rate, B={:,}").format(block_bytes)
        label((left + right) / 2, panel_top - 13, title, 17,
              text_anchor="middle")
        if kind == "time":
            values = [speed[(k, block_bytes)][field]
                      for k in EXPECTED_K
                      for field in ("wh1_lifecycle_ms", "wh2_lifecycle_ms")]
            low = max(0.005, min(values) * 0.75)
            high = max(values) * 1.25

            def y_coord(value):
                return panel_bottom - math.log(value / low) / math.log(high / low) * panel_h
        else:
            low, high = 0.84, 1.01

            def y_coord(value):
                return panel_bottom - (value - low) / (high - low) * panel_h

        def x_coord(k):
            return left + k / 64000.0 * panel_w

        if kind == "time":
            ticks = (low, math.sqrt(low * high), high)
            for value in ticks:
                y = y_coord(value)
                line(left, y, right, y)
                label(left - 9, y + 5, "{:.2g}".format(value), 11,
                      text_anchor="end")
        else:
            for value in (0.85, 0.90, 0.95, 1.0):
                y = y_coord(value)
                line(left, y, right, y)
                label(left - 9, y + 5, "{:g}%".format(value * 100), 11,
                      text_anchor="end")
        for tick in x_ticks:
            x = x_coord(tick)
            line(x, panel_top, x, panel_bottom, "#eef1f4")
            label(x, panel_bottom + 17, "{}k".format(tick // 1000) if tick else "0",
                  11, text_anchor="middle")
        label((left + right) / 2, panel_bottom + 35, "Original blocks (K)", 12,
              text_anchor="middle")

        for index, (codec, color, marker) in enumerate(
                (("wh1", COLORS[0], "circle"), ("wh2", COLORS[1], "rect"))):
            if kind == "recovery":
                eventual_field = codec + "_eventual_rate"
                eventual_points = [(x_coord(k), y_coord(recovery[(k, block_bytes)][
                    eventual_field])) for k in EXPECTED_K]
                element("polyline", points=" ".join(
                    "{:.3f},{:.3f}".format(x, y) for x, y in eventual_points),
                        fill="none", stroke=color, stroke_width=1.3,
                        stroke_dasharray="2 4", stroke_opacity=0.8)
            field = (codec + "_lifecycle_ms" if kind == "time" else
                     codec + "_exact_k_rate")
            points = [(x_coord(k), y_coord((speed if kind == "time" else recovery)[
                (k, block_bytes)][field])) for k in EXPECTED_K]
            element("polyline", points=" ".join("{:.3f},{:.3f}".format(x, y)
                                                for x, y in points),
                    fill="none", stroke=color, stroke_width=2.2,
                    stroke_dasharray="none" if index == 0 else "7 4")
            # Mark only every 1000 point to keep the large sweep readable;
            # preserve all small-K markers so those points remain identifiable.
            for k, (x, y) in zip(EXPECTED_K, points):
                if k not in SMALL_K and k % 5000:
                    continue
                if marker == "circle":
                    node = element("circle", cx=x, cy=y, r=3.2, fill=color)
                else:
                    node = element("rect", x=x - 3.2, y=y - 3.2,
                                   width=6.4, height=6.4, fill=color)
                if kind == "recovery":
                    ET.SubElement(node, "title").text = (
                        "{} K={} B={} exact-K {:.2f}%".format(
                            codec.upper(), k, block_bytes,
                            recovery[(k, block_bytes)][field] * 100.0))
    label(width / 2, 716,
          "Source: checked-in CSV snapshots · current benchmark build · diagnostic, not universal qualification",
          12, text_anchor="middle")
    return ET.tostring(svg, encoding="utf-8", xml_declaration=True) + b"\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="verify SVG is current")
    args = parser.parse_args()
    speed, recovery = read_data(
        DIRECTORY / "wh1-vs-wh2-large-k-speed.csv",
        DIRECTORY / "wh1-vs-wh2-large-k-recovery.csv")
    content = render(speed, recovery)
    output = DIRECTORY / "wh1-vs-wh2-large-k.svg"
    if args.check:
        if not output.is_file() or output.read_bytes() != content:
            parser.exit(1, "SVG is stale; run plot_wh1_vs_wh2_large.py\n")
        print("SVG matches source data")
    else:
        output.write_bytes(content)
        print("Wrote {}".format(output))


if __name__ == "__main__":
    main()
