"""HadamardLab – construction of Hadamard matrices, MOLS / orthogonal arrays and lattice designs.

Developed by Dr. Sukanta Dash, Senior Scientist, Division of Statistical Ecology and
Environmental Statistics, ICAR-IASRI, New Delhi-110012 (sukanta.iasri@gmail.com).

Run with::

    shiny run --reload app.py

The numerical work lives in the ``hadamard_lab`` package (pure NumPy, fully unit-tested);
this file only wires inputs to those functions and renders the results.
"""
from __future__ import annotations

import datetime as _dt
import hashlib
import json
import platform
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt                      # noqa: E402
from matplotlib.colors import ListedColormap         # noqa: E402
from matplotlib.patches import Patch                 # noqa: E402
import numpy as np                                   # noqa: E402
import pandas as pd                                  # noqa: E402
from shiny import App, reactive, render, req, ui     # noqa: E402

from hadamard_lab import __version__
from hadamard_lab import analysis as A
from hadamard_lab import constructions as C
from hadamard_lab import content as T
from hadamard_lab import mols as M
from hadamard_lab import planner as P
from hadamard_lab.gf import field as gf_field, modulus_string, prime_power

HERE = Path(__file__).parent
CSS = (HERE / "www" / "style.css").read_text(encoding="utf-8")
DATA_SHA = hashlib.sha256(C.DATA_FILE.read_bytes()).hexdigest()[:16]

# monospaced text output: `output_code` / `render.code` in current Shiny, `text_verbatim` in older releases
_output_code = getattr(ui, "output_code", None) or ui.output_text_verbatim
_render_code = getattr(render, "code", None) or render.text

PREVIEW = 40          # rows/columns shown in the on-screen text preview

# ==========================================================================
# small UI helpers
# ==========================================================================
def stat(label: str, value, sub: str = "", tone: str = ""):
    return ui.div(ui.div(label, class_="lab"), ui.div(str(value), class_="val"),
                  ui.div(sub, class_="sub") if sub else None, class_=f"stat {tone}".strip())


def banner(kind: str, *content):
    return ui.div(*content, class_=f"banner {kind}")


def kv_table(rows):
    return ui.tags.table(*[ui.tags.tr(ui.tags.td(str(k)), ui.tags.td(str(v))) for k, v in rows], class_="kv")


def df_grid(df: pd.DataFrame, filters: bool = False):
    return render.DataGrid(df, width="100%", filters=filters)


METHOD_CHOICES = {"auto": "Automatic (preferred applicable method)"}
METHOD_CHOICES.update({m.key: m.name for m in sorted(C.METHODS.values(), key=lambda m: m.priority)
                       if m.key != "trivial"})

FAMILIES = [
    ("Sylvester / Paley family", ("sylvester", "paley1", "paley2", "paley2_general", "trivial"), "#2b6cb0"),
    ("Block-array constructions", ("baumert_hall", "ehlich", "miyamoto"), "#2f855a"),
    ("Sequence-based (Williamson, Goethals–Seidel, SDS …)",
     ("williamson", "gs_turyn", "gs_base", "cooper_wallis", "sds"), "#d69e2e"),
    ("Kronecker product", ("kronecker",), "#805ad5"),
    ("Open – no matrix known", (), "#c53030"),
    ("Exists in literature, not implemented", (), "#8a94a0"),
]


def family_index(method: str, order: int) -> int:
    if not method:
        return 4 if order in P.KNOWN_OPEN else 5
    for i, (_, keys, _) in enumerate(FAMILIES):
        if method in keys:
            return i
    return 5


# ==========================================================================
# static pages
# ==========================================================================
def _coverage():
    rows = P.existence_table()
    ok = sum(1 for r in rows if r[3])
    return ok, len(rows)


COV_OK, COV_ALL = _coverage()

overview_page = ui.nav_panel(
    "Overview",
    ui.div(
        ui.h2("Combinatorial designs for experimental research"),
        ui.p("Generate and verify Hadamard matrices of order up to 1000 by thirteen classical constructions, "
             "complete sets of mutually orthogonal Latin squares and orthogonal arrays OA(s², s+1, s, 2), "
             "and square-lattice designs. Every object is checked against its defining identity before it is shown."),
        ui.p(ui.tags.b(f"Developed by {T.DEVELOPER}"), ui.tags.br(), T.AFFILIATION_LINE, ui.tags.br(),
             ui.tags.a(T.EMAIL, href=f"mailto:{T.EMAIL}", style="color:#f2c9bd"),
             style="margin-top:.8rem;opacity:.95;font-size:.92rem"),
        class_="hero"),
    ui.layout_columns(
        ui.card(ui.card_header("Hadamard matrices"),
                ui.p(f"{COV_OK} of the {COV_ALL} admissible orders 4, 8, …, 1000 are constructible. "
                     "The only unresolved orders ≤ 1000 are 668, 716 and 892 (no matrix is known to exist) and 876 "
                     "(known to exist, but the required sequences are not bundled)."),
                ui.p("Choose a method or let the planner pick one, normalise, view, verify and export "
                     "(CSV, text, JSON provenance, Plackett–Burman design).")),
        ui.card(ui.card_header("Existence explorer"),
                ui.p("A map of all orders up to 1000 showing which family of constructions reaches each order, with a "
                     "searchable table of all applicable methods.")),
        ui.card(ui.card_header("MOLS and orthogonal arrays"),
                ui.p("For every prime power s < 1000 the complete set of s − 1 MOLS over GF(s) and the corresponding "
                     "orthogonal array, with exhaustive verification for small s and sampled verification for large s.")),
        ui.card(ui.card_header("Lattice designs"),
                ui.p("Simple, triple and — for prime-power k — balanced square-lattice designs for k² treatments "
                     "with field book, concurrence structure and efficiency factor.")),
        col_widths=[3, 3, 3, 3]),
    ui.card(ui.card_header("How correctness is ensured"),
            ui.tags.ul(
                ui.tags.li("Each constructed matrix H is tested for H·Hᵀ = n·I (exact integer arithmetic) and displayed with the result."),
                ui.tags.li("The bundled sequence tables (Williamson, base, Turyn-type sequences, supplementary difference sets) "
                           "are re-verified against their defining equations."),
                ui.tags.li("A test-suite builds every constructible order with every applicable method and checks it independently."),
                ui.tags.li("Four formulas that fail when transcribed literally were re-derived and are documented under About."))),
    ui.div(f"HadamardLab {__version__} · Developed by {T.DEVELOPER}, {T.AFFILIATION_LINE} · {T.EMAIL} · data fingerprint {DATA_SHA}", class_="foot"),
)


def _methods_page():
    table = P.existence_table()
    panels = []
    for m in sorted(C.METHODS.values(), key=lambda m: m.priority):
        if m.key == "trivial":
            continue
        doc = T.METHOD_DOCS[m.key]
        orders = [r[0] for r in table if m.key in r[3]]
        shown = ", ".join(map(str, orders[:40])) + (" …" if len(orders) > 40 else "")
        body = [ui.p(ui.HTML(doc["statement"])),
                ui.p(ui.tags.b("Applies when: "), ui.HTML(doc["condition"])),
                ui.p(ui.tags.b("Reference: "), m.reference or "—"),
                ui.p(ui.tags.b(f"Orders ≤ 1000 reachable ({len(orders)}): "), shown or "none")]
        if doc["notes"]:
            body.append(ui.div(ui.HTML(doc["notes"]), class_="corr"))
        panels.append(ui.accordion_panel(m.name, *body, value=m.key))
    return ui.nav_panel(
        "Methods",
        ui.h3("Constructions"),
        ui.p("Statements follow the classical literature; ‘Orders reachable’ is computed from the actual planner, "
             "so it reflects what this installation can really build."),
        ui.accordion(*panels, id="methods_acc", open=False),
        ui.div(f"HadamardLab {__version__} · Developed by {T.DEVELOPER}, {T.AFFILIATION_LINE} · {T.EMAIL}", class_="foot"))


def _about_page():
    return ui.nav_panel(
        "About",
        ui.h3("About and validation"),
        ui.layout_columns(
            ui.card(ui.card_header("Developer and provenance"), ui.p(ui.HTML(T.DEVELOPER_CREDIT)), ui.p(ui.HTML(T.DATA_CREDIT)),
                    ui.p("Sequence data fingerprint (SHA-256, first 16 hex): ", ui.tags.code(DATA_SHA))),
            ui.card(ui.card_header("Implementation notes"),
                    *[ui.div(ui.tags.b(a + ": "), ui.HTML(b), class_="corr") for a, b in T.CORRECTIONS]),
            col_widths=[6, 6]),
        ui.card(ui.card_header("Scope and limitations"),
                ui.tags.ul(
                    ui.tags.li("Orders are limited to n ≤ 1000. Orders 668, 716 and 892 are open cases of the Hadamard conjecture."),
                    ui.tags.li("Finite fields use a primitive polynomial found by search; for prime powers q = pᵏ with k > 1 the "
                               "labelling of elements (hence the exact matrix) can differ from other software, but the "
                               "matrices are equivalent."),
                    ui.tags.li("Orthogonal arrays are generated from the complete set of MOLS over GF(s). Verification is "
                               "exhaustive for small s and sampled (and labelled as such) when s is large."),
                    ui.tags.li("Williamson-type matrices exist only for the stored orders; the app does not search for new ones."))),
        ui.card(ui.card_header("References"),
                ui.tags.ol(*[ui.tags.li(ui.HTML(r)) for r in T.REFERENCES], class_="ref")),
        ui.card(ui.card_header("Suggested citation"),
                ui.p(f"Dash, S. HadamardLab (v{__version__}): a verified generator of Hadamard matrices, MOLS, "
                     f"orthogonal arrays and lattice designs. {T.DIVISION}, {T.INSTITUTE}. "
                     f"Accessed {_dt.date.today().isoformat()}.")),
        ui.div(f"HadamardLab {__version__} · Developed by {T.DEVELOPER}, {T.AFFILIATION_LINE} · {T.EMAIL}", class_="foot"))


# ==========================================================================
# page assembly
# ==========================================================================
gen_sidebar = ui.sidebar(
    ui.h5("Parameters"),
    ui.input_numeric("order", "Order n", 12, min=1, max=P.ORDER_LIMIT, step=1),
    ui.output_ui("order_hint"),
    ui.input_select("method", "Construction method", METHOD_CHOICES, selected="auto"),
    ui.input_radio_buttons("norm", "Normalisation", A.NORMALISATIONS, selected="none"),
    ui.input_radio_buttons("fmt", "Output notation", A.FORMATS, selected="pm1"),
    ui.input_action_button("generate", "Generate and verify", class_="btn-primary w-100"),
    width=340)

generator_page = ui.nav_panel(
    "Hadamard generator",
    ui.layout_sidebar(
        gen_sidebar,
        ui.output_ui("gen_banner"),
        ui.output_ui("gen_stats"),
        ui.br(),
        ui.navset_card_tab(
            ui.nav_panel("Matrix",
                         ui.layout_columns(
                             ui.card(ui.card_header("Sign pattern"), ui.output_plot("matrix_plot", height="560px"),
                                     full_screen=True),
                             ui.card(ui.card_header("Text preview"), ui.output_ui("preview_note"),
                                     _output_code("matrix_text"), full_screen=True),
                             col_widths=[6, 6])),
            ui.nav_panel("Construction",
                         ui.layout_columns(
                             ui.card(ui.card_header("Construction path"), ui.output_ui("construction_ui")),
                             ui.card(ui.card_header("All methods for this order"),
                                     ui.output_data_frame("alt_table")),
                             col_widths=[5, 7])),
            ui.nav_panel("Verification",
                         ui.layout_columns(
                             ui.card(ui.card_header("Orthogonality check"), ui.output_ui("verify_ui")),
                             ui.card(ui.card_header("Properties"), ui.output_ui("props_ui")),
                             col_widths=[6, 6]),
                         ui.card(ui.card_header("Gram matrix H·Hᵀ and row-sum distribution"),
                                 ui.output_plot("gram_plot", height="380px"))),
            ui.nav_panel("Experimental design",
                         ui.card(ui.card_header("Plackett–Burman design from this matrix"),
                                 ui.output_ui("pb_ui"), ui.output_data_frame("pb_table"))),
            ui.nav_panel("Export",
                         ui.card(ui.card_header("Download"),
                                 ui.p("Files use the notation and normalisation selected in the sidebar."),
                                 ui.layout_columns(
                                     ui.download_button("dl_csv", "Matrix (CSV)", class_="btn-outline-secondary w-100"),
                                     ui.download_button("dl_txt", "Matrix (text)", class_="btn-outline-secondary w-100"),
                                     ui.download_button("dl_json", "Provenance (JSON)", class_="btn-outline-secondary w-100"),
                                     ui.download_button("dl_pb", "PB design (CSV)", class_="btn-outline-secondary w-100"),
                                     col_widths=[3, 3, 3, 3]))),
            id="gen_tabs"),
        fillable=False))

explorer_page = ui.nav_panel(
    "Existence explorer",
    ui.layout_columns(
        ui.card(ui.card_header("Orders 4 … 1000 by construction family"),
                ui.output_plot("exist_plot", height="760px"), full_screen=True),
        ui.div(ui.card(ui.card_header("Summary"), ui.output_ui("exist_summary")),
               ui.card(ui.card_header("Orders and applicable methods"), ui.output_data_frame("exist_table"))),
        col_widths=[7, 5]))

mols_sidebar = ui.sidebar(
    ui.h5("Parameters"),
    ui.input_numeric("m_s", "Order s (prime power ≤ 997)", 7, min=2, max=M.S_LIMIT, step=1),
    ui.input_radio_buttons("m_view", "Output", {"squares": "Individual Latin squares",
                                                "oa": "Orthogonal array OA(s², s+1, s, 2)"}, selected="squares"),
    ui.panel_conditional("input.m_view == 'squares'",
                         ui.input_numeric("m_k", "Square number k (1 … s−1)", 1, min=1, max=M.S_LIMIT - 1, step=1)),
    ui.panel_conditional("input.m_view == 'oa'",
                         ui.input_numeric("m_from", "First row to export (0-based)", 0, min=0, step=1),
                         ui.input_numeric("m_count", "Rows to export", 100000, min=1, step=1000)),
    ui.input_action_button("m_generate", "Generate and verify", class_="btn-primary w-100"),
    width=340)

mols_page = ui.nav_panel(
    "MOLS & orthogonal arrays",
    ui.layout_sidebar(
        mols_sidebar,
        ui.output_ui("m_banner"),
        ui.output_ui("m_stats"),
        ui.br(),
        ui.card(ui.card_header("Result"), ui.output_ui("m_note"), ui.output_data_frame("m_table")),
        ui.card(ui.card_header("Download"), ui.output_ui("m_dl_note"),
                ui.download_button("m_dl", "Download CSV", class_="btn-outline-secondary")),
        fillable=False))

lat_sidebar = ui.sidebar(
    ui.h5("Parameters"),
    ui.input_numeric("l_v", "Number of treatments v = k² (9 … 961)", 25, min=9, max=961, step=1),
    ui.input_numeric("l_r", "Replicates r", 3, min=2, max=32, step=1),
    ui.output_ui("l_hint"),
    ui.input_numeric("l_show", "Replicate to display", 1, min=1, max=32, step=1),
    ui.input_action_button("l_generate", "Generate and verify", class_="btn-primary w-100"),
    width=340)

lattice_page = ui.nav_panel(
    "Lattice designs",
    ui.layout_sidebar(
        lat_sidebar,
        ui.output_ui("l_banner"),
        ui.output_ui("l_stats"),
        ui.br(),
        ui.layout_columns(
            ui.card(ui.card_header("Blocks of the selected replicate"), ui.output_ui("l_blocks"), full_screen=True),
            ui.card(ui.card_header("Structure"), ui.output_ui("l_structure")),
            col_widths=[7, 5]),
        ui.card(ui.card_header("Field book"), ui.output_data_frame("l_book"),
                ui.download_button("l_dl", "Download field book (CSV)", class_="btn-outline-secondary")),
        fillable=False))

MATHJAX = ui.head_content(
    ui.tags.meta(name="viewport", content="width=device-width, initial-scale=1"),
    ui.tags.script("window.MathJax = {tex:{inlineMath:[['\\\\(','\\\\)']]}};"),
    ui.tags.script(src="https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-chtml.js", id="MathJax-script", async_=True),
)

app_ui = ui.page_navbar(
    overview_page, generator_page, explorer_page, mols_page, lattice_page, _methods_page(), _about_page(),
    title=ui.span("HadamardLab", ui.tags.small("combinatorial designs"), class_="brand"),
    id="main_nav",
    header=ui.TagList(ui.tags.style(CSS), MATHJAX),
    window_title="HadamardLab – Hadamard matrices, MOLS and lattice designs | Dr. Sukanta Dash",
    fillable=False,
)


# ==========================================================================
# server
# ==========================================================================
@dataclass
class GenResult:
    ok: bool
    error: str = ""
    construction: Optional[P.Construction] = None
    H: Optional[np.ndarray] = None
    ver: Optional[A.Verification] = None
    props: Optional[A.Properties] = None
    norm: str = "none"
    fmt: str = "pm1"
    meta: dict = field(default_factory=dict)


def _tree_dict(node: P.Node) -> dict:
    return {"order": node.order, "method": node.method, "label": node.label,
            "parameters": {k: v for k, v in node.params.items()},
            "notes": node.notes, "children": [_tree_dict(c) for c in node.children]}


def server(input, output, session):
    # ------------------------------------------------------------------
    # Hadamard generator
    # ------------------------------------------------------------------
    @render.ui
    def order_hint():
        n = input.order()
        if n is None:
            return ui.div("Enter an order.", class_="hint hint-warn")
        ok, msg = P.validate_order(n)
        if not ok:
            return ui.div(msg, class_="hint hint-warn")
        n = int(n)
        ms = P.feasible_methods(n)
        if not ms:
            if n in P.KNOWN_OPEN:
                return ui.div(f"Order {n} is an open case: no Hadamard matrix is known.", class_="hint hint-warn")
            return ui.div(f"No construction for order {n} is implemented.", class_="hint hint-warn")
        return ui.div(ui.tags.b(f"{len(ms)} applicable: "), ", ".join(C.METHODS[m].name for m in ms), class_="hint")

    @reactive.calc
    @reactive.event(input.generate, ignore_none=False)
    def gen() -> GenResult:
        n = input.order()
        ok, msg = P.validate_order(n) if n is not None else (False, "Enter an order.")
        if not ok:
            return GenResult(False, msg)
        method = input.method()
        method = None if method in (None, "auto") else method
        try:
            con = P.construct(int(n), method)
        except ValueError as exc:
            return GenResult(False, str(exc))
        norm, fmt = input.norm(), input.fmt()
        H = A.normalise(con.matrix, norm)
        ver = A.verify(H)
        props = A.properties(H)
        meta = {
            "application": f"HadamardLab {__version__}",
            "developer": {"name": T.DEVELOPER, "designation": T.DESIGNATION, "division": T.DIVISION,
                          "institute": T.INSTITUTE, "email": T.EMAIL},
            "generated_utc": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
            "order": int(n),
            "requested_method": method or "auto",
            "construction": _tree_dict(con.tree),
            "normalisation": norm,
            "notation": fmt,
            "verification": ver.as_dict(),
            "data_fingerprint_sha256_16": DATA_SHA,
            "environment": {"python": platform.python_version(), "numpy": np.__version__},
        }
        return GenResult(True, "", con, H, ver, props, norm, fmt, meta)

    @render.ui
    def gen_banner():
        r = gen()
        if not r.ok:
            return banner("err", ui.tags.b("Cannot generate: "), r.error)
        if r.ver.passed:
            return banner("ok", ui.tags.b("Verified. "),
                          f"H·Hᵀ = {r.ver.order}·I holds exactly for the order-{r.ver.order} matrix built by ",
                          P.describe_node(r.construction.tree), ".")
        return banner("err", ui.tags.b("Verification failed. "), "The matrix is not Hadamard — please report this case.")

    @render.ui
    def gen_stats():
        r = gen()
        req(r.ok)
        tone = "ok" if r.ver.passed else "bad"
        return ui.layout_columns(
            stat("Order", r.ver.order, f"{r.ver.order**2:,} entries"),
            stat("Method", r.construction.tree.label,
                 "recursive" if r.construction.tree.children else "direct"),
            stat("H·Hᵀ = nI", "passed" if r.ver.passed else "FAILED",
                 f"max deviation {max(r.ver.max_row_deviation, r.ver.max_col_deviation)}", tone),
            stat("log₁₀|det H|", f"{r.ver.log10_determinant:.2f}", "Hadamard's bound attained"),
            stat("Time", f"{r.construction.seconds*1000:.0f} ms", "construction"),
            col_widths=[2, 4, 2, 2, 2])

    @render.plot(alt="Black and white image of the Hadamard matrix; light cells are +1, dark cells are −1.")
    def matrix_plot():
        r = gen()
        req(r.ok)
        H, n = r.H, r.H.shape[0]
        fig, ax = plt.subplots(figsize=(6, 6), dpi=110)
        ax.imshow(H, cmap=ListedColormap(["#1d2b3a", "#f4f1ea"]), vmin=-1, vmax=1, interpolation="nearest")
        if n <= 40:
            ax.set_xticks(np.arange(-.5, n, 1), minor=True)
            ax.set_yticks(np.arange(-.5, n, 1), minor=True)
            ax.grid(which="minor", color="#8b97a3", lw=.35)
            ax.tick_params(which="minor", length=0)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f"Hadamard matrix of order {n}  (light = +1, dark = −1)", fontsize=10, color="#17324d")
        fig.tight_layout()
        return fig

    @render.ui
    def preview_note():
        r = gen()
        req(r.ok)
        n = r.H.shape[0]
        txt = (f"First {PREVIEW}×{PREVIEW} entries of {n}×{n}; download for the full matrix."
               if n > PREVIEW else f"Full {n}×{n} matrix.")
        return ui.div(txt, class_="hint")

    @_render_code
    def matrix_text():
        r = gen()
        req(r.ok)
        sub = r.H[:PREVIEW, :PREVIEW]
        return A.to_text(sub, r.fmt)

    @render.ui
    def construction_ui():
        r = gen()
        req(r.ok)

        def render_node(node: P.Node):
            kids = [ui.tags.ul(*[render_node(c) for c in node.children], class_="tree")] if node.children else []
            return ui.tags.li(
                ui.span(f"order {node.order}", class_="ord"), " — ",
                ui.span(P.describe_node(node), class_="nm"),
                *[ui.span(nt, class_="note") for nt in node.notes], *kids)
        ref = C.METHODS[r.construction.tree.method].reference
        return ui.div(ui.tags.ul(render_node(r.construction.tree), class_="tree"),
                      ui.p(ui.tags.b("Reference: "), ref or "—", class_="hint"))

    @render.data_frame
    def alt_table():
        r = gen()
        req(r.ok)
        rows = [{"Method": name, "Applicable": "yes" if ok else "no", "Why": why}
                for _, name, ok, why in P.alternatives(r.H.shape[0])]
        return df_grid(pd.DataFrame(rows))

    @render.ui
    def verify_ui():
        r = gen()
        req(r.ok)
        v = r.ver
        yn = lambda b: ui.span("yes" if b else "no", class_="yes" if b else "no")
        return ui.tags.table(
            ui.tags.tr(ui.tags.td("All entries ±1"), ui.tags.td(yn(v.entries_pm1))),
            ui.tags.tr(ui.tags.td("Rows orthogonal (H·Hᵀ = nI)"), ui.tags.td(yn(v.row_orthogonal))),
            ui.tags.tr(ui.tags.td("Columns orthogonal (Hᵀ·H = nI)"), ui.tags.td(yn(v.column_orthogonal))),
            ui.tags.tr(ui.tags.td("Max |H·Hᵀ − nI|"), ui.tags.td(str(v.max_row_deviation))),
            ui.tags.tr(ui.tags.td("Max |Hᵀ·H − nI|"), ui.tags.td(str(v.max_col_deviation))),
            ui.tags.tr(ui.tags.td("log₁₀|det H| = (n/2)·log₁₀ n"), ui.tags.td(f"{v.log10_determinant:.4f}")),
            class_="kv")

    @render.ui
    def props_ui():
        r = gen()
        req(r.ok)
        return kv_table(r.props.as_rows())

    @render.plot(alt="Left: Gram matrix H times H transpose. Right: histogram of row sums.")
    def gram_plot():
        r = gen()
        req(r.ok)
        H = r.H
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3.6), dpi=110)
        im = a1.imshow(A.gram(H), cmap="Blues", interpolation="nearest")
        a1.set_title("H·Hᵀ (equals n on the diagonal, 0 elsewhere)", fontsize=9)
        a1.set_xticks([])
        a1.set_yticks([])
        fig.colorbar(im, ax=a1, fraction=.046, pad=.03)
        rs = A.row_sums(H)
        a2.hist(rs, bins=min(30, max(5, len(set(rs.tolist())))), color="#24486d", edgecolor="white")
        a2.set_title("Row sums", fontsize=9)
        a2.set_xlabel("sum of a row")
        a2.set_ylabel("rows")
        fig.tight_layout()
        return fig

    @reactive.calc
    def pb():
        r = gen()
        req(r.ok)
        return A.plackett_burman(r.H)

    @render.ui
    def pb_ui():
        D = pb()
        n = D.shape[0]
        return ui.div(
            ui.p(f"{n} runs × {D.shape[1]} two-level factors; the columns are balanced and pairwise orthogonal "
                 f"(checked: {'yes' if A.design_is_orthogonal(D) else 'NO'}). Deleting the constant column of a normalised "
                 "Hadamard matrix gives the classical Plackett–Burman layout."),
            ui.p("Showing the first 24 runs and 12 factors.", class_="hint"))

    @render.data_frame
    def pb_table():
        D = pb().iloc[:24, :12].reset_index()
        return df_grid(D)

    def _fname(ext):
        try:
            r = gen()
            return f"hadamard_{r.H.shape[0]}.{ext}" if r.ok else f"hadamard.{ext}"
        except Exception:
            return f"hadamard.{ext}"

    @render.download(filename=lambda: _fname("csv"))
    def dl_csv():
        r = gen()
        yield A.to_csv(r.H, r.fmt) if r.ok else "No matrix generated\n"

    @render.download(filename=lambda: _fname("txt"))
    def dl_txt():
        r = gen()
        yield (A.to_text(r.H, r.fmt) + "\n") if r.ok else "No matrix generated\n"

    @render.download(filename=lambda: _fname("json"))
    def dl_json():
        r = gen()
        yield A.metadata_json(r.meta) if r.ok else "{}"

    @render.download(filename=lambda: _fname("pb.csv"))
    def dl_pb():
        r = gen()
        yield pb().to_csv() if r.ok else "No matrix generated\n"

    # ------------------------------------------------------------------
    # Existence explorer
    # ------------------------------------------------------------------
    @render.plot(alt="Grid of all orders 4 to 1000 coloured by the construction family that reaches them.")
    def exist_plot():
        rows = P.existence_table()
        cols = 10
        nrow = int(np.ceil(len(rows) / cols))
        grid = np.full((nrow, cols), np.nan)
        fig, ax = plt.subplots(figsize=(7.2, 8.2), dpi=105)
        for idx, (n, status, first, ms) in enumerate(rows):
            grid[idx // cols, idx % cols] = family_index(first, n)
        cmap = ListedColormap([c for _, _, c in FAMILIES])
        ax.imshow(grid, cmap=cmap, vmin=0, vmax=len(FAMILIES) - 1, aspect="auto")
        for idx, (n, *_rest) in enumerate(rows):
            ax.text(idx % cols, idx // cols, str(n), ha="center", va="center", fontsize=6.4, color="white")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title("Orders 4 – 1000 (first applicable construction family)", fontsize=10, color="#17324d")
        fig.legend(handles=[Patch(color=c, label=l) for l, _, c in FAMILIES], loc="lower center",
                   ncol=2, fontsize=7.2, frameon=False)
        fig.tight_layout(rect=[0, 0.11, 1, 1])
        return fig

    @render.ui
    def exist_summary():
        rows = P.existence_table()
        counts = {}
        for n, st, first, ms in rows:
            counts[family_index(first, n)] = counts.get(family_index(first, n), 0) + 1
        items = [ui.tags.li(ui.span(class_="legend-dot", style=f"background:{FAMILIES[i][2]}"),
                            f"{FAMILIES[i][0]}: ", ui.tags.b(str(counts.get(i, 0)))) for i in range(len(FAMILIES))]
        return ui.div(ui.tags.ul(*items, style="list-style:none;padding-left:0;margin:0"),
                      ui.p(f"{COV_OK} of {COV_ALL} orders are constructible. Open cases: "
                           f"{', '.join(map(str, P.KNOWN_OPEN))}.", class_="hint"))

    @render.data_frame
    def exist_table():
        rows = [{"Order": n, "Status": st, "Preferred": C.METHODS[first].name if first else "—",
                 "Applicable methods": len(ms), "All methods": "; ".join(C.METHODS[m].name for m in ms)}
                for n, st, first, ms in P.existence_table()]
        return df_grid(pd.DataFrame(rows), filters=True)

    # ------------------------------------------------------------------
    # MOLS / orthogonal arrays
    # ------------------------------------------------------------------
    @reactive.calc
    @reactive.event(input.m_generate, ignore_none=False)
    def mres():
        s = input.m_s()
        ok, msg = M.validate_s(s) if s is not None else (False, "Enter an order.")
        if not ok:
            return {"ok": False, "error": msg}
        s = int(s)
        f = gf_field(s)
        big = s > 64
        ok_set, msg_set = M.verify_mols_set(s, sample=20 if big else 200)
        ok_oa, msg_oa = M.verify_oa(s, max_pairs=40 if big else None)
        return {"ok": True, "s": s, "p": f.p, "k": f.k, "poly": modulus_string(f),
                "verified": ok_set and ok_oa, "msgs": [msg_set, msg_oa], "view": input.m_view()}

    @render.ui
    def m_banner():
        r = mres()
        if not r["ok"]:
            return banner("err", ui.tags.b("Cannot generate: "), r["error"])
        if r["verified"]:
            return banner("ok", ui.tags.b("Verified. "), " · ".join(r["msgs"]))
        return banner("err", ui.tags.b("Verification failed. "), " · ".join(r["msgs"]))

    @render.ui
    def m_stats():
        r = mres()
        req(r["ok"])
        s = r["s"]
        return ui.layout_columns(
            stat("Order s", s, f"GF({r['p']}{'^' + str(r['k']) if r['k'] > 1 else ''})"),
            stat("MOLS in complete set", s - 1, "upper bound s − 1 attained"),
            stat("Orthogonal array", f"{s*s} × {s+1}", f"{s} levels, strength 2"),
            stat("Field", "prime" if r["k"] == 1 else "extension", r["poly"]),
            col_widths=[3, 3, 3, 3])

    @render.ui
    def m_note():
        r = mres()
        req(r["ok"])
        s = r["s"]
        if input.m_view() == "squares":
            k = int(input.m_k() or 1)
            if not 1 <= k <= s - 1:
                return ui.div(f"Square number must lie between 1 and {s-1}.", class_="hint hint-warn")
            lim = 60
            return ui.div(f"Latin square L{k}[i, j] = {k}·i + j in GF({s}); rows i and columns j run over 0 … {s-1}."
                          + (f" Showing the top-left {lim}×{lim} block." if s > lim else ""), class_="hint")
        return ui.div(f"Rows (i, j) ∈ GF({s})²; columns i, j, then k·i + j for k = 1 … {s-1}. "
                      f"Showing the first 40 rows and up to 14 columns of {s*s:,} × {s+1}.", class_="hint")

    @render.data_frame
    def m_table():
        r = mres()
        req(r["ok"])
        s = r["s"]
        if input.m_view() == "squares":
            k = int(input.m_k() or 1)
            req(1 <= k <= s - 1)
            L = M.latin_square(s, k)[:60, :60]
            return df_grid(pd.DataFrame(L, columns=[str(j) for j in range(L.shape[1])]).reset_index()
                           .rename(columns={"index": "i \\ j"}))
        X = M.oa_rows(s, 0, 40)[:, :14]
        names = ["i", "j"] + [f"k={k}" for k in range(1, s)]
        return df_grid(pd.DataFrame(X, columns=names[:X.shape[1]]))

    def _m_payload():
        r = mres()
        if not r["ok"]:
            return "No data\n", "mols.csv"
        s = r["s"]
        if input.m_view() == "squares":
            k = int(input.m_k() or 1)
            if s * s * (s - 1) <= M.MAX_MATERIALISE:
                parts = []
                for kk in range(1, s):
                    parts.append(f"# L{kk}: i*{kk}+j in GF({s})\n" +
                                 "\n".join(",".join(map(str, row)) for row in M.latin_square(s, kk)))
                return "\n".join(parts) + "\n", f"mols_{s}_all.csv"
            k = min(max(k, 1), s - 1)
            return ("\n".join(",".join(map(str, row)) for row in M.latin_square(s, k)) + "\n"), f"mols_{s}_L{k}.csv"
        start = int(input.m_from() or 0)
        count = int(input.m_count() or 1)
        cells_per_row = s + 1
        count = max(1, min(count, max(1, M.MAX_MATERIALISE // cells_per_row)))
        X = M.oa_rows(s, start, start + count)
        header = ",".join(["i", "j"] + [f"k{k}" for k in range(1, s)])
        body = "\n".join(",".join(map(str, row)) for row in X)
        return header + "\n" + body + "\n", f"oa_{s}_rows_{start}_{start + len(X)}.csv"

    @render.ui
    def m_dl_note():
        r = mres()
        req(r["ok"])
        s = r["s"]
        if input.m_view() == "squares":
            if s * s * (s - 1) <= M.MAX_MATERIALISE:
                return ui.div(f"The download contains all {s-1} squares, one after another.", class_="hint")
            return ui.div("The complete set is too large to export in one file; the download contains the selected square.",
                          class_="hint")
        total = s * s
        return ui.div(f"The full array has {total:,} rows. Use the row range in the sidebar for large s "
                      f"(at most {M.MAX_MATERIALISE // (s + 1):,} rows per file).", class_="hint")

    @render.download(filename=lambda: _m_payload()[1])
    def m_dl():
        yield _m_payload()[0]

    # ------------------------------------------------------------------
    # Lattice designs
    # ------------------------------------------------------------------
    @render.ui
    def l_hint():
        ok, msg = M.validate_lattice(input.l_v()) if input.l_v() is not None else (False, "Enter a number.")
        if not ok:
            return ui.div(msg, class_="hint hint-warn")
        k = int(round(int(input.l_v()) ** 0.5))
        mx = M.max_replicates(k)
        extra = " (balanced lattice when r = k + 1)" if mx == k + 1 else " (cyclic Latin square only)"
        return ui.div(f"k = {k}; replicates 2 … {mx}{extra}.", class_="hint")

    @reactive.calc
    @reactive.event(input.l_generate, ignore_none=False)
    def lres():
        v = input.l_v()
        ok, msg = M.validate_lattice(v) if v is not None else (False, "Enter a number.")
        if not ok:
            return {"ok": False, "error": msg}
        k = int(round(int(v) ** 0.5))
        r = input.l_r()
        if r is None or not (2 <= int(r) <= M.max_replicates(k)):
            return {"ok": False, "error": f"For k = {k} the number of replicates must be between 2 and {M.max_replicates(k)}."}
        r = int(r)
        return {"ok": True, "k": k, "r": r, "summary": M.lattice_summary(k, r),
                "blocks": M.lattice_blocks(k, r), "book": M.lattice_field_book(k, r)}

    @render.ui
    def l_banner():
        x = lres()
        if not x["ok"]:
            return banner("err", ui.tags.b("Cannot generate: "), x["error"])
        s = x["summary"]
        if s["each_treatment_once_per_replicate"]:
            kind = "balanced lattice (every pair of treatments meets exactly once)" if s["balanced"] else \
                {2: "simple lattice", 3: "triple lattice"}.get(x["r"], f"partially balanced lattice with {x['r']} replicates")
            return banner("ok", ui.tags.b("Verified. "),
                          f"Each treatment occurs once per replicate and blocks have size {x['k']}: {kind}.")
        return banner("err", ui.tags.b("Verification failed."), " A replicate is not a partition of the treatments.")

    @render.ui
    def l_stats():
        x = lres()
        req(x["ok"])
        s = x["summary"]
        return ui.layout_columns(
            stat("Treatments", s["treatments"], f"k = {x['k']}"),
            stat("Replicates", s["replicates"], f"{s['blocks']} blocks of {s['block_size']}"),
            stat("Plots", f"{s['plots']:,}"),
            stat("Average efficiency", f"{s['average_efficiency']:.4f}", "harmonic mean of canonical factors"),
            col_widths=[3, 3, 3, 3])

    @render.ui
    def l_blocks():
        x = lres()
        req(x["ok"])
        rep = int(input.l_show() or 1)
        if not 1 <= rep <= x["r"]:
            return ui.div(f"Choose a replicate between 1 and {x['r']}.", class_="hint hint-warn")
        blocks = x["blocks"][rep - 1]
        head = ui.tags.tr(ui.tags.th("Block"), *[ui.tags.th(str(j + 1)) for j in range(x["k"])])
        rows = [ui.tags.tr(ui.tags.th(str(b + 1)), *[ui.tags.td(str(t)) for t in blk]) for b, blk in enumerate(blocks)]
        return ui.div(ui.tags.table(head, *rows, class_="blk"), style="overflow:auto;max-height:520px")

    @render.ui
    def l_structure():
        x = lres()
        req(x["ok"])
        s = x["summary"]
        conc = ", ".join(f"λ = {a}: {b:,} pairs" for a, b in sorted(s["concurrence_counts"].items()))
        eff = ", ".join(f"{a:.4f} (×{b})" for a, b in sorted(s["canonical_efficiencies"].items(), reverse=True))
        return kv_table([("Pair concurrences", conc), ("Balanced (BIBD)", "yes" if s["balanced"] else "no"),
                         ("Canonical efficiency factors", eff),
                         ("Average efficiency", f"{s['average_efficiency']:.4f}")])

    @render.data_frame
    def l_book():
        x = lres()
        req(x["ok"])
        return df_grid(x["book"].head(500), filters=True)

    @render.download(filename=lambda: "lattice_field_book.csv")
    def l_dl():
        x = lres()
        yield x["book"].to_csv(index=False) if x["ok"] else "No design generated\n"


app = App(app_ui, server)


if __name__ == "__main__":
    # Allows `python app.py` (or the IDE "Run" button) to start the server and open the browser.
    # Equivalent to:  shiny run --reload app.py
    from shiny import run_app

    run_app("app:app", host="127.0.0.1", port=8000, launch_browser=True)
