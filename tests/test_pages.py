"""The HTML pages over a folder of QC figures (qualitymetrics.pages)."""
from __future__ import annotations

import base64
import json
import re
import shutil
import subprocess

import numpy as np
import pytest

matplotlib = pytest.importorskip("matplotlib")
matplotlib.use("Agg")

from qualitymetrics import KilosortResults, pages  # noqa: E402
from qualitymetrics.plots.units import template_sites  # noqa: E402

from .conftest import UV_PER_BIT  # noqa: E402

#: A 3 x 2 PNG: signature, IHDR length and type, width, height, the rest.
PNG = (b"\x89PNG\r\n\x1a\n" + b"\x00\x00\x00\rIHDR" + (3).to_bytes(4, "big")
       + (2).to_bytes(4, "big") + b"\x08\x02\x00\x00\x00")


@pytest.fixture
def ks(sorter_output):
    return KilosortResults.load(sorter_output, uv_per_bit=UV_PER_BIT,
                                label="TEST 2026-01-01 imec0_shank0")


def _folder(tmp_path, names, name="shank"):
    out = tmp_path / name
    out.mkdir()
    for n in names:
        (out / n).write_bytes(PNG if n.endswith(".png") else b"x")
    return out


def _payload(text: str) -> dict:
    m = re.search(r'<script type="application/json" id="tpl-data">(.*?)</script>',
                  text, re.S)
    return json.loads(m.group(1))


def test_figures_appear_in_reading_order_and_unknown_ones_after(tmp_path):
    out = _folder(tmp_path, ["raw_traces.png", "amp_depth_scatter.png",
                             "drift_map.png", "zz_new_figure.png",
                             "dredge_driftmap.pdf", "Thumbs.db"])
    text = pages.write_shank_page(out, title="T").read_text(encoding="utf-8")
    order = [m.group(1) for m in re.finditer(r'<img src="([^"]+)\.png"', text)]
    assert order == ["amp_depth_scatter", "drift_map", "raw_traces", "zz_new_figure"]
    assert 'href="dredge_driftmap.pdf"' in text
    assert "Thumbs.db" not in text


def test_a_missing_figure_is_listed_with_the_reports_reason(tmp_path):
    out = _folder(tmp_path, ["amp_depth_scatter.png"])
    text = pages.write_shank_page(
        out, title="T", skipped={"raw_traces": "no archived recording found",
                                 "phy_tsvs": "microvolts unknown"}
    ).read_text(encoding="utf-8")
    not_drawn = text[text.index('id="not-drawn"'):]
    assert "<td>raw_traces</td><td>no archived recording found</td>" in not_drawn
    assert "<td>phy_tsvs</td><td>microvolts unknown</td>" in not_drawn
    # Listed figures with no recorded reason are still named.
    assert "<td>drift_map</td><td>no reason recorded</td>" in not_drawn
    assert "amp_depth_scatter</td>" not in not_drawn


def test_images_reserve_their_size_from_the_png_header(tmp_path):
    out = _folder(tmp_path, ["drift_map.png"])
    text = pages.write_shank_page(out, title="T").read_text(encoding="utf-8")
    assert '<img src="drift_map.png" width="3" height="2"' in text
    assert pages._png_size(tmp_path / "nope.png") is None


def test_the_title_and_names_are_escaped(tmp_path):
    out = _folder(tmp_path, ["a&b.png"])
    text = pages.write_shank_page(out, title="<x> & y").read_text(encoding="utf-8")
    assert "<title>&lt;x&gt; &amp; y</title>" in text
    assert 'src="a%26b.png"' in text


def test_embedded_templates_decode_to_what_the_grid_draws(tmp_path, ks):
    out = _folder(tmp_path, ["templates_grid.png"])
    text = pages.write_shank_page(out, title="T", ks=ks,
                                  exclude_units=[2]).read_text(encoding="utf-8")
    data = _payload(text)
    assert data["exclude"] == [2]
    assert data["units"] == [int(u) for u in ks.unit_ids]
    nc, nt = data["n_sites"], data["n_samples"]
    waves = np.frombuffer(base64.b64decode(data["waves"]), dtype=np.int8)
    waves = waves.reshape(len(data["units"]), nt, nc)
    for i, u in enumerate(data["units"]):
        chans, wave = template_sites(ks, u, 32)
        assert data["chans"][i * nc:(i + 1) * nc] == chans.tolist()
        peak = np.abs(wave).max()
        # Each unit at int8 relative to its own peak: within half a step.
        np.testing.assert_allclose(waves[i] / 127.0 * peak, wave,
                                   atol=0.5 / 127.0 * peak + 1e-6)
        assert data["peak_site"][i] == int(ks.peak_channel[u])
        assert data["amp_uv"][i] == pytest.approx(ks.unit_amplitude_uv[u])


def test_the_payload_cannot_end_its_script_tag(tmp_path, ks):
    ks.label = "</script><b>"
    data = pages.template_data(ks)
    data["exclude"] = []
    out = _folder(tmp_path, ["templates_grid.png"])
    page = pages.write_shank_page(out, title=ks.label, ks=ks).read_text(encoding="utf-8")
    body = page[page.index("<body>"):]
    assert body.count("</script>") == 2          # the data and the drawing code


def test_a_page_without_templates_says_why_and_keeps_the_figure(tmp_path, ks, monkeypatch):
    def broken(*_a, **_k):
        raise ValueError("no units with templates")
    monkeypatch.setattr(pages, "template_data", broken)
    out = _folder(tmp_path, ["templates_grid.png"])
    text = pages.write_shank_page(out, title="T", ks=ks).read_text(encoding="utf-8")
    assert 'src="templates_grid.png"' in text
    assert "No redraw control: ValueError: no units with templates" in text
    assert "tpl-data" not in text


def test_the_session_page_links_each_shank_page_by_probe_and_shank(tmp_path):
    out = _folder(tmp_path, ["session_raster_passing.png",
                             "session_raster_all_units.png"], name="qc_plots")
    for name in ("p1_g0_imec0_shank0", "p1_g0_imec1_shank3"):
        (out / name).mkdir()
        (out / name / "index.html").write_text("x", encoding="utf-8")
    (out / "p1_g0_imec0_shank1").mkdir()             # figures, no page yet
    (out / "_previous").mkdir()
    (out / "_previous" / "index.html").write_text("x", encoding="utf-8")
    text = pages.write_session_page(
        out, title="S", shank_folders=["p1_g0_imec0_shank0", "p1_g0_imec0_shank2"]
    ).read_text(encoding="utf-8")
    assert 'href="p1_g0_imec0_shank0/index.html"' in text
    assert 'href="p1_g0_imec1_shank3/index.html"' in text   # found on disk
    assert 'href="p1_g0_imec0_shank1/"' in text and "no page yet" in text
    assert "No figures" in text                              # shank2: named, absent
    assert "_previous" not in text
    order = [m.group(1) for m in re.finditer(r'<img src="([^"]+)\.png"', text)]
    assert order == ["session_raster_all_units", "session_raster_passing"]


def test_build_report_writes_the_shank_page(sorter_output, tmp_path):
    from qualitymetrics.report import build_report

    # No recording: the sorter-only figures are enough, and far quicker.
    result = build_report(sorter_output, tmp_path / "qc", uv_per_bit=UV_PER_BIT)
    assert "index.html" in result.made
    text = (tmp_path / "qc" / "index.html").read_text(encoding="utf-8")
    for name in result.made:
        if (tmp_path / "qc" / f"{name}.png").exists():
            assert f'src="{name}.png"' in text
    # The static grid's units are the ones the first redraws leave out.
    assert sorted(_payload(text)["exclude"]) == sorted(
        u for u in (int(x) for x in KilosortResults.load(
            sorter_output, uv_per_bit=UV_PER_BIT).unit_ids))


def test_a_failing_page_is_recorded_not_raised(tmp_path, ks, monkeypatch):
    from qualitymetrics import report

    def broken(*_a, **_k):
        raise OSError("share went away")
    monkeypatch.setattr(pages, "write_shank_page", broken)
    result = report.ReportResult(out_dir=tmp_path)
    result.made["templates_grid"] = "x"
    report.write_page(result, ks)
    assert result.skipped["index.html"] == "OSError: share went away"
    assert result.made == {"templates_grid": "x"}


# ------------------------------------------------- the drawing code itself

NODE_HARNESS = r"""
const fs = require("fs");
const html = fs.readFileSync(process.argv[2], "utf8");
const data = html.match(/<script type="application\/json" id="tpl-data">([\s\S]*?)<\/script>/)[1];
const code = html.match(/<script>\n([\s\S]*?)<\/script>/)[1];
function el(tag) {
  return {tag, attrs: {}, children: [], listeners: {}, textContent: "", value: "",
          setAttribute(k, v) { this.attrs[k] = String(v); },
          appendChild(c) { this.children.push(c); return c; },
          addEventListener(t, f) { this.listeners[t] = f; },
          set innerHTML(_v) { this.children = []; }};
}
const byId = {"tpl-data": Object.assign(el("script"), {textContent: data})};
for (const id of ["tpl-next", "tpl-show", "tpl-ids", "tpl-panels", "tpl-status"]) byId[id] = el(id);
global.document = {getElementById: id => byId[id] || null,
                   createElementNS: (_ns, tag) => el(tag)};
eval(code);
const out = {draws: []};
function panels() {
  return byId["tpl-panels"].children.map(svg => ({
    title: svg.children.filter(c => c.tag === "text").slice(0, 2).map(t => t.textContent),
    lines: svg.children.filter(c => c.tag === "polyline").map(p => ({
      stroke: p.attrs.stroke,
      pts: p.attrs.points.split(" ").map(s => s.split(",").map(Number))}))}));
}
for (let k = 0; k < 3; k++) { byId["tpl-next"].listeners.click(); out.draws.push(panels()); }
byId["tpl-ids"].value = "3, 1, 3, 999";
byId["tpl-show"].listeners.click();
out.typed = panels();
out.status = byId["tpl-status"].textContent;
console.log(JSON.stringify(out));
"""


def _fake_ks(n_units=20, n_sites=40):
    """Enough units for several draws of six; units 0-14 pass QC."""
    from types import SimpleNamespace

    rng = np.random.default_rng(3)
    pos = np.array([[27.0 if i % 2 else 59.0, (i // 2) * 15.0] for i in range(n_sites)])
    t = np.arange(61)
    wave = -np.exp(-((t - 20) ** 2) / 8.0) + 0.4 * np.exp(-((t - 32) ** 2) / 20.0)
    templates = np.zeros((n_units, 61, n_sites), dtype=np.float32)
    peak = rng.integers(0, n_sites, n_units)
    for u in range(n_units):
        templates[u, :, peak[u]] = wave * (10 + u)
        templates[u, :, max(peak[u] - 2, 0)] += 0.4 * wave * (10 + u)
    return SimpleNamespace(
        templates_unwhitened=templates, channel_positions=pos,
        unit_ids=np.arange(n_units), peak_channel=peak, uv_per_bit=2.0, fs=30000.0,
        unit_amplitude_uv={u: float(np.ptp(templates[u]) * 2.0) for u in range(n_units)},
        qc_pass={u: u < 15 for u in range(n_units)})


def _run_page_script(tmp_path, page):
    harness = tmp_path / "harness.js"
    harness.write_text(NODE_HARNESS, encoding="utf-8")
    done = subprocess.run(["node", str(harness), str(page)], capture_output=True,
                          text=True, encoding="utf-8", timeout=60)
    assert done.returncode == 0, done.stderr
    return json.loads(done.stdout)


def _unit_of(panel):
    return int(panel["title"][1].split()[1].rstrip(","))


needs_node = pytest.mark.skipif(shutil.which("node") is None,
                                reason="node is not installed")


@needs_node
def test_the_browser_draws_passing_units_without_repeats_static_ones_last(tmp_path):
    """The page's own script, run in node over a stand-in document."""
    ks = _fake_ks()
    out = _folder(tmp_path, ["templates_grid.png"])
    page = pages.write_shank_page(out, title="T", ks=ks, exclude_units=[0, 1, 2, 3])
    got = _run_page_script(tmp_path, page)
    draws = [[_unit_of(p) for p in d] for d in got["draws"]]
    assert all(p["title"][0] == "Random" for d in got["draws"] for p in d)
    # 15 pass; 4 of them are on the static grid, so the first round is the
    # other 11: six, then five, then the round ends and a new one begins.
    assert all(u < 15 for d in draws for u in d)
    assert all(len(d) == 6 and len(set(d)) == 6 for d in draws)
    assert set(draws[0]) <= set(range(4, 15))
    # Two draws of six cover all eleven: no unit came back inside the round.
    assert set(range(4, 15)) <= set(draws[0]) | set(draws[1])
    # Typed units, any QC verdict, deduplicated, unknown ones named.
    assert [_unit_of(p) for p in got["typed"]] == [3, 1]
    assert [p["title"][0] for p in got["typed"]] == ["Passed QC", "Passed QC"]
    assert "No template for 999" in got["status"]


@needs_node
def test_the_browser_draws_the_peak_trace_as_the_grid_does(tmp_path, ks):
    """Unit 3's peak-site trace in the panel's pixels has the template's shape."""
    out = _folder(tmp_path, ["templates_grid.png"])
    page = pages.write_shank_page(out, title="T", ks=ks)
    panel = _run_page_script(tmp_path, page)["typed"][0]
    assert _unit_of(panel) == 3
    peak_lines = [ln for ln in panel["lines"] if ln["stroke"] == "crimson"]
    assert len(peak_lines) == 1
    assert len(panel["lines"]) == min(32, len(ks.channel_positions))
    pts = np.array(peak_lines[0]["pts"])
    chans, wave = template_sites(ks, 3, 32)
    j = list(chans).index(int(ks.peak_channel[3]))
    expected = wave[:, j] - wave[0, j]
    drawn = -(pts[:, 1] - pts[0, 1])                  # SVG y grows downward
    np.testing.assert_allclose(drawn / np.abs(drawn).max(),
                               expected / np.abs(expected).max(), atol=0.03)
    assert len(pts) == ks.templates_unwhitened.shape[1]
    assert np.all(np.diff(pts[:, 0]) > 0)            # time left to right


# ------------------------------------------------------- review, 2026-11-16

def test_a_page_written_again_keeps_the_reports_reasons(tmp_path):
    out = _folder(tmp_path, ["amp_depth_scatter.png"])
    pages.write_shank_page(out, title="T", skipped={"raw_traces": "withheld: map"})
    again = pages.write_shank_page(out, title="T").read_text(encoding="utf-8")
    assert "<td>raw_traces</td><td>withheld: map</td>" in again
    assert "not_drawn.json" not in again                 # kept, not listed


def test_report_write_page_with_kept_reasons_reads_them(tmp_path, ks):
    from qualitymetrics import report

    out = _folder(tmp_path, ["templates_grid.png"])
    first = report.ReportResult(out_dir=out, skipped={"raw_traces": "withheld: map"})
    report.write_page(first, ks)
    report.write_page(report.ReportResult(out_dir=out), ks, kept_reasons=True)
    text = (out / "index.html").read_text(encoding="utf-8")
    assert "<td>raw_traces</td><td>withheld: map</td>" in text


def test_the_session_page_names_missing_session_figures(tmp_path):
    out = _folder(tmp_path, ["session_raster_all_units.png"], name="qc_plots")
    text = pages.write_session_page(out, title="S").read_text(encoding="utf-8")
    assert "<td>session_raster_passing</td><td>no reason recorded</td>" in text
    assert "<td>session_raster_all_units</td>" not in text


def test_a_figure_cannot_take_the_controls_ids(tmp_path, ks):
    out = _folder(tmp_path, ["templates_grid.png", "tpl-data.png", "units.png"])
    text = pages.write_shank_page(out, title="T", ks=ks).read_text(encoding="utf-8")
    assert text.count('id="tpl-data"') == 1
    assert text.count('id="units"') == 0 or text.count('id="units"') == 1
    assert 'id="fig-tpl-data"' in text


def test_a_template_with_no_finite_value_is_named_not_drawn_as_real(ks, recwarn):
    ks.templates_unwhitened[1][:] = np.nan
    data = pages.template_data(ks)
    assert data["bad"] == [1]
    waves = np.frombuffer(base64.b64decode(data["waves"]), dtype=np.int8)
    assert not [w for w in recwarn if "invalid value" in str(w.message)]
    assert waves.size
