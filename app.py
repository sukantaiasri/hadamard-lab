"""HadamardLab – construction of Hadamard matrices, MOLS / orthogonal arrays and lattice designs.

Developed by Dr. Sukanta Dash, Senior Scientist, Division of Statistical Ecology and
Environmental Statistics, ICAR-IASRI, New Delhi-110012 (sukantadash.iasri@icar.org.in).
Identity, citation and copyright: hadamard_lab/content.py; information pages: hadamard_lab/site.py.

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
from hadamard_lab import seqlab as SL
from hadamard_lab import sequences as SQ
from hadamard_lab import site as S
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

# (long label, method keys, colour, short label for the plot legend)
FAMILIES = [
    ("Sylvester / Paley family", ("sylvester", "paley1", "paley2", "paley2_general", "trivial"), "#2b6cb0",
     "Sylvester / Paley"),
    ("Two-circulant and difference-set cores (Golay pairs, Legendre pairs, Singer / twin-prime / Hall sets)",
     ("golay2c", "legendre", "dscore"), "#0e7c86", "Two-circulant / difference-set cores"),
    ("Block-array constructions (Baumert–Hall, Ehlich, Miyamoto, Scarpis)", ("baumert_hall", "ehlich", "miyamoto", "scarpis"),
     "#2f855a", "Block arrays (Baumert–Hall, Ehlich, Miyamoto, Scarpis)"),
    ("Sequence-based (Williamson, Goethals–Seidel, T-sequences, SDS, Legendre quadruples)",
     ("williamson", "gs_turyn", "gs_base", "gs_golay", "gs_tseq", "cooper_wallis", "sds", "lquad"), "#d69e2e",
     "Sequence-based (Williamson, Goethals–Seidel, SDS …)"),
    ("Kronecker product", ("kronecker",), "#805ad5", "Kronecker product"),
    ("Open – no matrix known", (), "#c53030", "Open – no matrix known"),
]


def family_index(method: str, order: int) -> int:
    if not method:
        return len(FAMILIES) - 1
    for i, (_, keys, _, _) in enumerate(FAMILIES):
        if method in keys:
            return i
    return len(FAMILIES) - 1


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
        ui.p(f"Generate and verify Hadamard matrices of every admissible order up to 1000 by {len(C.METHODS) - 1} constructions "
             "— classical (Sylvester, Paley, Williamson, Goethals–Seidel, Baumert–Hall, Ehlich, Miyamoto, Cooper–Wallis, "
             "supplementary difference sets) and, new in version 2.0, Golay pairs, Legendre pairs and quadruples, cyclic "
             "Hadamard difference sets and the Scarpis–Đoković construction — together with complete sets of mutually "
             "orthogonal Latin squares, orthogonal arrays OA(s², s+1, s, 2) and square-lattice designs. Every object is "
             "checked against its defining identity before it is shown."),
        ui.p(ui.tags.b(f"Developed by {T.DEVELOPER}"), ui.tags.br(), T.AFFILIATION_LINE, ui.tags.br(),
             ui.tags.a(T.EMAIL, href=f"mailto:{T.EMAIL}", style="color:#f2c9bd"),
             style="margin-top:.8rem;opacity:.95;font-size:.92rem"),
        class_="hero"),
    ui.layout_columns(
        ui.card(ui.card_header("Hadamard matrices"),
                ui.p(f"All {COV_OK} of the {COV_ALL} admissible orders 4, 8, …, 1000 are constructible. Orders 668, 716 and 892 — "
                     "the smallest open cases of the Hadamard conjecture from 2005 until August 2026 — are built from the "
                     "sequences announced by Alpöge, Voinov and Reynolds-Haertle (2026); order 876 from T-sequences of length 73."),
                ui.p("Choose a method or let the planner pick one, normalise, view, verify and export "
                     "(CSV, text, JSON provenance, Plackett–Burman design).")),
        ui.card(ui.card_header("Construction lab"),
                ui.p("New in 2.0: build matrices from their ingredients — Golay pairs, T-sequences (with or without "
                     "Williamson matrices), Legendre pairs and quadruples, Singer / twin-prime / Hall difference sets and "
                     "the Scarpis–Đoković map H(q+1) → H(q(q+1)) — and see the autocorrelation identity each one satisfies.")),
        ui.card(ui.card_header("Existence explorer"),
                ui.p("A map of all orders up to 1000 showing which family of constructions reaches each order, with a "
                     "searchable table of all applicable methods.")),
        ui.card(ui.card_header("MOLS and orthogonal arrays"),
                ui.p("For every prime power s < 1000 the complete set of s − 1 MOLS over GF(s) and the corresponding "
                     "orthogonal array, with exhaustive verification for small s and sampled verification for large s.")),
        ui.card(ui.card_header("Lattice designs"),
                ui.p("Simple, triple and — for prime-power k — balanced square-lattice designs for k² treatments "
                     "with field book, concurrence structure and efficiency factor.")),
        col_widths=[3, 3, 2, 2, 2]),
    ui.card(ui.card_header("How correctness is ensured"),
            ui.tags.ul(
                ui.tags.li("Each constructed matrix H is tested for H·Hᵀ = n·I (exact integer arithmetic) and displayed with the result."),
                ui.tags.li("The bundled sequence tables (Williamson, base, Turyn-type sequences, supplementary difference sets) "
                           "are re-verified against their defining equations."),
                ui.tags.li("Golay pairs, Legendre pairs and quadruples, T-sequences and difference sets are checked against their "
                           "autocorrelation identities before they are used; the matrices of orders 668, 716 and 892 recovered from "
                           "the 2026 announcement are re-verified on every run."),
                ui.tags.li("A test-suite builds every constructible order with every applicable method and checks it independently."),
                ui.tags.li("Four formulas that fail when transcribed literally were re-derived and are documented under About."))),
    S.cite_panel(),
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
        title = ui.span(m.name, ui.span("new in 2.0", class_="badge-new") if m.key in C.NEW_IN_V2 else None)
        panels.append(ui.accordion_panel(title, *body, value=m.key))
    return ui.nav_panel(
        "Methods",
        ui.h3("Constructions"),
        ui.p(f"{len(C.METHODS) - 1} constructions; those marked ‘new in 2.0’ can also be run from their ingredients on the "
             "Construction lab tab. Statements follow the classical literature; ‘Orders reachable’ is computed from the "
             "actual planner, so it reflects what this installation can really build."),
        ui.accordion(*panels, id="methods_acc", open=False))


def _about_page():
    return ui.nav_panel(
        "About HadamardLab",
        ui.h3("About and validation"),
        ui.layout_columns(
            ui.card(ui.card_header("Developer and provenance"), ui.p(ui.HTML(T.DEVELOPER_CREDIT)), ui.p(ui.HTML(T.DATA_CREDIT)),
                    ui.p("Sequence data fingerprint (SHA-256, first 16 hex): ", ui.tags.code(DATA_SHA))),
            ui.card(ui.card_header("Implementation notes"),
                    *[ui.div(ui.tags.b(a + ": "), ui.HTML(b), class_="corr") for a, b in T.CORRECTIONS]),
            col_widths=[6, 6]),
        ui.card(ui.card_header("Scope and limitations"),
                ui.tags.ul(
                    ui.tags.li("Orders are limited to n ≤ 1000, where every admissible order is now constructible. A Hadamard matrix "
                               "is known at every order below 2000 (Alpöge, Voinov & Reynolds-Haertle, August 2026; verified by Egan, "
                               f"Kotsireas, London & Ó Catháin, 2026); the smallest open order is {P.SMALLEST_OPEN_ORDER}."),
                    ui.tags.li("Legendre pairs and quadruples of composite length are bundled only where a pair/quadruple is stored "
                               "(found by the built-in search); the app does not search for new ones at run time."),
                    ui.tags.li("Finite fields use a primitive polynomial found by search; for prime powers q = pᵏ with k > 1 the "
                               "labelling of elements (hence the exact matrix) can differ from other software, but the "
                               "matrices are equivalent."),
                    ui.tags.li("Orthogonal arrays are generated from the complete set of MOLS over GF(s). Verification is "
                               "exhaustive for small s and sampled (and labelled as such) when s is large."),
                    ui.tags.li("Williamson-type matrices exist only for the stored orders; the app does not search for new ones."))),
        ui.card(ui.card_header("Recent references (2013–2026) added in version 2.0"),
                ui.p("Surveys, new constructions and computational results that the new methods and the coverage statements rest on.",
                     class_="hint"),
                ui.tags.ul(*[ui.tags.li(ui.HTML(r)) for r in T.recent_references()], class_="ref")),
        ui.card(ui.card_header("References"),
                ui.tags.ol(*[ui.tags.li(ui.HTML(r)) for r in T.REFERENCES], class_="ref")),
        ui.card(ui.card_header("How to cite"),
                ui.p("Please cite HadamardLab as shown on the ", ui.input_action_link("about_to_cite", "How to cite"),
                     " page, where the citation, BibTeX and RIS files are available.")))


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


# ---- Construction lab (new in 2.0) --------------------------------------------------
seq_sidebar = ui.sidebar(
    ui.h5("Construction"),
    ui.input_select("sq_kind", "Build from", SL.KINDS, selected="golay"),
    ui.output_ui("sq_params"),
    ui.input_action_button("sq_generate", "Build and verify", class_="btn-primary w-100"),
    ui.div(ui.p("Each construction starts from its ingredients (a pair, a quadruple, a difference set, T-sequences or "
                "an input matrix). The identity the ingredients must satisfy is checked first, then the assembled matrix.",
                class_="hint"), style="margin-top:.8rem"),
    width=360)

seq_page = ui.nav_panel(
    "Construction lab",
    ui.layout_sidebar(
        seq_sidebar,
        ui.output_ui("sq_banner"),
        ui.output_ui("sq_stats"),
        ui.br(),
        ui.navset_card_tab(
            ui.nav_panel("Ingredients",
                         ui.layout_columns(
                             ui.card(ui.card_header("Sequences"), ui.output_ui("sq_desc"), _output_code("sq_ingredients"),
                                     full_screen=True),
                             ui.card(ui.card_header("Autocorrelation identity"), ui.output_ui("sq_checks"),
                                     ui.output_plot("sq_corr_plot", height="300px")),
                             col_widths=[7, 5])),
            ui.nav_panel("Matrix",
                         ui.layout_columns(
                             ui.card(ui.card_header("Sign pattern"), ui.output_plot("sq_matrix_plot", height="560px"),
                                     full_screen=True),
                             ui.card(ui.card_header("Text preview"), ui.output_ui("sq_preview_note"),
                                     _output_code("sq_matrix_text"), full_screen=True),
                             col_widths=[6, 6])),
            ui.nav_panel("Verification",
                         ui.layout_columns(
                             ui.card(ui.card_header("Orthogonality check"), ui.output_ui("sq_verify_ui")),
                             ui.card(ui.card_header("Properties"), ui.output_ui("sq_props_ui")),
                             col_widths=[6, 6])),
            ui.nav_panel("Export",
                         ui.card(ui.card_header("Download"),
                                 ui.p("The matrix is exported as constructed (±1 notation); the ingredients with their "
                                      "autocorrelation values and the reference as JSON or text."),
                                 ui.layout_columns(
                                     ui.download_button("sq_dl_csv", "Matrix (CSV)", class_="btn-outline-secondary w-100"),
                                     ui.download_button("sq_dl_txt", "Matrix (text)", class_="btn-outline-secondary w-100"),
                                     ui.download_button("sq_dl_seq", "Ingredients (JSON)", class_="btn-outline-secondary w-100"),
                                     ui.download_button("sq_dl_seqtxt", "Ingredients (text)", class_="btn-outline-secondary w-100"),
                                     col_widths=[3, 3, 3, 3]))),
            id="sq_tabs"),
        fillable=False))

MATHJAX = ui.head_content(
    ui.tags.meta(name="viewport", content="width=device-width, initial-scale=1"),
    ui.tags.script("window.MathJax = {tex:{inlineMath:[['\\\\(','\\\\)']]}};"),
    ui.tags.script(src="https://cdn.jsdelivr.net/npm/mathjax@3.2.2/es5/tex-chtml.js", id="MathJax-script", async_=True),
)

app_ui = ui.page_navbar(
    overview_page, generator_page, seq_page, explorer_page, mols_page, lattice_page, _methods_page(),
    ui.nav_menu("About", _about_page(), S.cite_page(), S.legal_page(), value="about_menu"),
    S.contact_page(),
    title=ui.span(ui.tags.img(src="icar-logo.png", class_="brand-logo", alt=""), "HadamardLab",
                  ui.tags.small("combinatorial designs"), class_="brand"),
    id="main_nav",
    header=ui.TagList(S.FONTS, ui.tags.style(CSS), MATHJAX, ui.tags.script(ui.HTML(S.SCRIPT)), S.inst_strip()),
    footer=S.footer(__version__),
    window_title=f"HadamardLab \u2013 Hadamard matrices, sequences, MOLS and lattice designs | {T.INSTITUTE_SHORT}",
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
    # ---- citation (access date = today in India) and information-page links ----
    @render.ui
    def cite_main():
        return ui.HTML(T.citation_html())

    @render.ui
    def cite_home():
        return ui.HTML(T.citation_html())

    @render.text
    def cite_bibtex():
        return T.citation_bibtex()

    @render.download(filename="HadamardLab_citation.bib")
    def dl_bib():
        yield T.citation_bibtex()

    @render.download(filename="HadamardLab_citation.ris")
    def dl_ris():
        yield T.citation_ris()

    _update_nav = getattr(ui, "update_navset", None) or ui.update_navs

    def _goto(target):
        _update_nav("main_nav", selected=target)

    @reactive.effect
    @reactive.event(input.ft_cite, input.home_to_cite, input.about_to_cite)
    def _go_cite():
        _goto("cite")

    @reactive.effect
    @reactive.event(input.ft_copyright, input.ft_disclaimer)
    def _go_legal():
        _goto("legal")

    @reactive.effect
    @reactive.event(input.ft_contact)
    def _go_contact():
        _goto("contact")

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
    # Construction lab (new in 2.0)
    # ------------------------------------------------------------------
    def _sel(name, default=None):
        try:
            v = input[name]()
        except Exception:
            return default
        return default if v is None else v

    @render.ui
    def sq_params():
        kind = input.sq_kind()
        if kind == "golay":
            return ui.input_select("sq_g", "Golay length g = 2^a·10^b·26^c", SL.golay_choices(), selected="10")
        if kind == "tseq":
            L = int(_sel("sq_t", "21") or 21)
            if str(L) not in SL.tseq_choices():
                L = 21
            return ui.TagList(
                ui.input_select("sq_t", "T-sequence length t", SL.tseq_choices(), selected=str(L)),
                ui.input_select("sq_tsrc", "Source of the T-sequences", SL.tseq_source_choices(L)),
                ui.input_select("sq_w", "Combine with Williamson matrices", SL.williamson_choices(L), selected="0"))
        if kind == "legendre":
            L = int(_sel("sq_l", "13") or 13)
            if str(L) not in SL.legendre_choices():
                L = 13
            return ui.TagList(
                ui.input_select("sq_l", "Length ℓ (odd)", SL.legendre_choices(), selected=str(L)),
                ui.input_select("sq_lsrc", "Source of the pair", SL.legendre_source_choices(L)))
        if kind == "lquad":
            return ui.input_select("sq_m", "Length m (even)", SL.lquad_choices(), selected="10")
        if kind == "dscore":
            L = int(_sel("sq_dl", "15") or 15)
            if str(L) not in SL.dscore_choices():
                L = 15
            return ui.TagList(
                ui.input_select("sq_dl", "Length ℓ of the difference set", SL.dscore_choices(), selected=str(L)),
                ui.input_select("sq_dsrc", "Difference set", SL.dscore_source_choices(L)))
        if kind == "scarpis":
            q = int(_sel("sq_q", "7") or 7)
            if str(q) not in SL.scarpis_choices():
                q = 7
            return ui.TagList(
                ui.input_select("sq_q", "Prime power q ≡ 3 (mod 4)", SL.scarpis_choices(), selected=str(q)),
                ui.input_select("sq_sinput", "Input Hadamard matrix of order q + 1", SL.scarpis_input_choices(q), selected="auto"))
        return ui.div()

    @reactive.calc
    @reactive.event(input.sq_generate, ignore_none=False)
    def sq() -> SL.SeqResult:
        kind = input.sq_kind()
        if kind == "golay":
            prm = {"g": _sel("sq_g", "10")}
        elif kind == "tseq":
            t = _sel("sq_t", "21")
            src = _sel("sq_tsrc", "")
            if not src or not src.endswith(":" + str(t)) and src.split(":")[0] not in ("golay2",):
                opts = list(SL.tseq_source_choices(int(t)))
                src = src if src in opts else (opts[0] if opts else "")
            prm = {"t": t, "tsrc": src, "w": _sel("sq_w", "0")}
        elif kind == "legendre":
            L = _sel("sq_l", "13")
            opts = list(SL.legendre_source_choices(int(L)))
            src = _sel("sq_lsrc", "")
            prm = {"l": L, "lsrc": src if src in opts else (opts[0] if opts else "")}
        elif kind == "lquad":
            prm = {"m": _sel("sq_m", "10")}
        elif kind == "dscore":
            L = _sel("sq_dl", "15")
            opts = list(SL.dscore_source_choices(int(L)))
            src = _sel("sq_dsrc", "")
            prm = {"l": L, "dsrc": src if src in opts else (opts[0] if opts else "")}
        else:
            q = _sel("sq_q", "7")
            m = _sel("sq_sinput", "auto")
            prm = {"q": q, "sinput": m if m in SL.scarpis_input_choices(int(q)) else "auto"}
        return SL.run(kind, prm)

    def _sq_method_name(r: SL.SeqResult) -> str:
        if r.kind == "tseq":
            if r.params.get("w", "0") not in ("0", 0):
                return "Cooper–Wallis"
            src = str(r.params.get("tsrc", "")).split(":")[0]
            return C.METHODS[{"base": "gs_base", "turyn": "gs_turyn", "tseq": "gs_tseq"}.get(src, "gs_golay")].name
        return C.METHODS[SL.KIND_TO_METHOD[r.kind]].name

    @render.ui
    def sq_banner():
        r = sq()
        if not r.ok:
            return banner("err", ui.tags.b("Cannot build: "), r.error or "an ingredient check failed.")
        return banner("ok", ui.tags.b("Verified. "), f"{r.title}: H·Hᵀ = {r.order}·I holds exactly for the order-{r.order} matrix.")

    @render.ui
    def sq_stats():
        r = sq()
        req(r.ok)
        return ui.layout_columns(
            stat("Order", r.order, f"{r.order**2:,} entries"),
            stat("Construction", _sq_method_name(r), "from ingredients"),
            stat("Ingredient identity", "passed" if all(c[2] for c in r.checks[:-1]) else "FAILED",
                 r.corr_kind + " autocorrelation" if r.corr_kind not in ("none", "group") else "input matrix", "ok"),
            stat("H·Hᵀ = nI", "passed" if r.checks and r.checks[-1][2] else "FAILED", "exact integer check", "ok"),
            stat("Time", f"{r.seconds*1000:.0f} ms", "construction"),
            col_widths=[2, 4, 2, 2, 2])

    @render.ui
    def sq_desc():
        r = sq()
        req(r.ok)
        return ui.div(ui.p(r.description), ui.p(ui.tags.b("Reference: "), r.reference, class_="hint"))

    @_render_code
    def sq_ingredients():
        r = sq()
        req(r.ok)
        return SL.ingredients_text(r)

    @render.ui
    def sq_checks():
        r = sq()
        req(r.ok)
        yn = lambda ok: "yes" if ok else "no"
        rows = [ui.tags.tr(ui.tags.td(lab), ui.tags.td(ui.span(val, class_="yes" if ok else "no"))) for lab, val, ok in r.checks]
        return ui.tags.table(*rows, class_="kv")

    @render.plot(alt="Total autocorrelation of the ingredient sequences at every shift.")
    def sq_corr_plot():
        r = sq()
        req(r.ok and r.corr is not None)
        c = r.corr
        fig, ax = plt.subplots(dpi=100)
        ax.bar(np.arange(1, len(c)), c[1:], color="#24486d", width=0.8)
        ax.axhline(r.corr_target, color="#b4472e", lw=1, ls="--")
        ax.set_xlabel("shift s")
        ax.set_ylabel("sum of autocorrelations")
        lo, hi = min(int(c[1:].min()), r.corr_target), max(int(c[1:].max()), r.corr_target)
        ax.set_ylim(lo - 1.5, hi + 1.5)
        ax.set_title(f"Σ {r.corr_kind} autocorrelations, s ≠ 0 (peak {int(c[0])} omitted)\ntarget {r.corr_target} at every shift",
                     fontsize=8.5)
        ax.tick_params(labelsize=8)
        fig.tight_layout()
        return fig

    @render.plot(alt="Black and white image of the constructed Hadamard matrix.")
    def sq_matrix_plot():
        r = sq()
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
        ax.set_title(f"{r.title} — order {n}", fontsize=9, color="#17324d")
        fig.tight_layout()
        return fig

    @render.ui
    def sq_preview_note():
        r = sq()
        req(r.ok)
        n = r.H.shape[0]
        return ui.div(f"First {PREVIEW}×{PREVIEW} entries of {n}×{n}; download for the full matrix." if n > PREVIEW
                      else f"Full {n}×{n} matrix.", class_="hint")

    @_render_code
    def sq_matrix_text():
        r = sq()
        req(r.ok)
        return A.to_text(r.H[:PREVIEW, :PREVIEW], "sign")

    @render.ui
    def sq_verify_ui():
        r = sq()
        req(r.ok)
        v = A.verify(r.H)
        yn = lambda b: ui.span("yes" if b else "no", class_="yes" if b else "no")
        return ui.tags.table(
            ui.tags.tr(ui.tags.td("All entries ±1"), ui.tags.td(yn(v.entries_pm1))),
            ui.tags.tr(ui.tags.td("Rows orthogonal (H·Hᵀ = nI)"), ui.tags.td(yn(v.row_orthogonal))),
            ui.tags.tr(ui.tags.td("Columns orthogonal (Hᵀ·H = nI)"), ui.tags.td(yn(v.column_orthogonal))),
            ui.tags.tr(ui.tags.td("Max |H·Hᵀ − nI|"), ui.tags.td(str(v.max_row_deviation))),
            ui.tags.tr(ui.tags.td("log₁₀|det H| = (n/2)·log₁₀ n"), ui.tags.td(f"{v.log10_determinant:.4f}")),
            class_="kv")

    @render.ui
    def sq_props_ui():
        r = sq()
        req(r.ok)
        return kv_table(A.properties(r.H).as_rows())

    def _sq_fname(ext):
        try:
            r = sq()
            return f"hadamard_{r.order}_{r.kind}.{ext}" if r.ok else f"hadamard.{ext}"
        except Exception:
            return f"hadamard.{ext}"

    @render.download(filename=lambda: _sq_fname("csv"))
    def sq_dl_csv():
        r = sq()
        yield A.to_csv(r.H, "pm1") if r.ok else "No matrix generated\n"

    @render.download(filename=lambda: _sq_fname("txt"))
    def sq_dl_txt():
        r = sq()
        yield (A.to_text(r.H, "pm1") + "\n") if r.ok else "No matrix generated\n"

    @render.download(filename=lambda: _sq_fname("ingredients.json"))
    def sq_dl_seq():
        r = sq()
        yield SL.ingredients_json(r) if r.ok else "{}"

    @render.download(filename=lambda: _sq_fname("ingredients.txt"))
    def sq_dl_seqtxt():
        r = sq()
        yield (r.title + "\n" + SL.ingredients_text(r) + "\n") if r.ok else "No construction\n"

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
        cmap = ListedColormap([c for _, _, c, _ in FAMILIES])
        ax.imshow(grid, cmap=cmap, vmin=0, vmax=len(FAMILIES) - 1, aspect="auto")
        for idx, (n, *_rest) in enumerate(rows):
            ax.text(idx % cols, idx // cols, str(n), ha="center", va="center", fontsize=6.4, color="white")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title("Orders 4 – 1000 (first applicable construction family)", fontsize=10, color="#17324d")
        fig.legend(handles=[Patch(color=c, label=sh) for _, _, c, sh in FAMILIES], loc="lower center",
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
                      ui.p(f"{COV_OK} of {COV_ALL} orders are constructible"
                           + (f"; open cases: {', '.join(map(str, P.KNOWN_OPEN))}." if P.KNOWN_OPEN else
                              f". No order ≤ 1000 is open: a Hadamard matrix is known at every order below 2000 "
                              f"(Alpöge, Voinov & Reynolds-Haertle 2026; Egan et al. 2026); the smallest open order is "
                              f"{P.SMALLEST_OPEN_ORDER}."), class_="hint"))

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


app = App(app_ui, server, static_assets=HERE / "www")


if __name__ == "__main__":
    # Allows `python app.py` (or the IDE "Run" button) to start the server and open the browser.
    # Equivalent to:  shiny run --reload app.py
    from shiny import run_app

    run_app("app:app", host="127.0.0.1", port=8000, launch_browser=True)
