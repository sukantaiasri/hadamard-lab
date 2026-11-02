"""Shared information pages for HadamardLab: How to cite, Copyright & disclaimer,
Contact, institutional strip and footer. Identity and citation live in content.py."""
from __future__ import annotations

import html

from shiny import ui

from . import content as T

ICONS = {
    "copy": '<rect x="9" y="9" width="13" height="13" rx="2"/><path d="M5 15H4a2 2 0 0 1-2-2V4a2 2 0 0 1 2-2h9a2 2 0 0 1 2 2v1"/>',
    "mail": '<rect x="2" y="4" width="20" height="16" rx="2"/><path d="m22 6-10 7L2 6"/>',
    "phone": '<path d="M22 16.9v3a2 2 0 0 1-2.2 2 19.8 19.8 0 0 1-8.6-3.1 19.5 19.5 0 0 1-6-6A19.8 19.8 0 0 1 2.1 4.2 2 2 0 0 1 4.1 2h3a2 2 0 0 1 2 1.7c.1.9.4 1.8.7 2.7a2 2 0 0 1-.5 2.1L8 9.8a16 16 0 0 0 6 6l1.3-1.3a2 2 0 0 1 2.1-.4c.9.3 1.8.6 2.7.7a2 2 0 0 1 1.7 2z"/>',
    "pin": '<path d="M20 10c0 6-8 12-8 12S4 16 4 10a8 8 0 0 1 16 0z"/><circle cx="12" cy="10" r="3"/>',
    "globe": '<circle cx="12" cy="12" r="10"/><path d="M2 12h20M12 2a15 15 0 0 1 0 20M12 2a15 15 0 0 0 0 20"/>',
}


def icon(name: str):
    return ui.HTML('<svg class="ico" viewBox="0 0 24 24" width="15" height="15" fill="none" stroke="currentColor" '
                   'stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true">'
                   f"{ICONS[name]}</svg>")


def copy_button(target: str, label: str = "Copy"):
    return ui.tags.button(icon("copy"), ui.span(label, class_="copy-label"), type="button",
                          class_="btn btn-sm btn-outline-secondary copy-btn", **{"data-copy-target": target})


def page_head(title: str, lead: str = "", kicker: str = ""):
    return ui.div(ui.div(kicker, class_="page-kicker") if kicker else None,
                  ui.h2(title, class_="page-title"),
                  ui.p(lead, class_="page-lead") if lead else None, class_="page-head")


def section(title: str, *content):
    return ui.div(ui.h3(title, class_="info-section-title"), *content, class_="info-section")


def mailto(kind: str) -> str:
    from urllib.parse import quote
    subject = f"{T.SHORT_NAME} - " + ("error report" if kind == "error" else "comments/suggestions")
    return f"mailto:{T.EMAIL}?subject={quote(subject)}"


# --------------------------------------------------------------------------
# pages
# --------------------------------------------------------------------------
def cite_page():
    ack = (f"The designs used in this study were obtained using {T.SHORT_NAME} (Dash, {T.YEAR}), "
           f"developed at the {T.INSTITUTE_FULL}, New Delhi.")
    return ui.nav_panel(
        "How to cite",
        ui.div(
            page_head("How to cite", kicker="Citation",
                      lead="If you use this software, or matrices and designs obtained from it, in a publication, "
                           "thesis, report or presentation, please cite it as shown below. The access date is "
                           "filled in automatically with today's date."),
            section("Recommended citation",
                    ui.div(ui.div(ui.output_ui("cite_main", inline=True), id="cite_plain", class_="citation-text"),
                           ui.div(copy_button("cite_plain", "Copy citation"), class_="citation-actions"),
                           class_="citation-box"),
                    ui.p("Format: Author (Year). Title. URL (accessed on date). Edit the style if required to "
                         "match the journal or publisher.", class_="small-note mt-2")),
            section("Reference manager formats",
                    ui.div(ui.div(ui.span("BibTeX"), copy_button("cite_bibtex_pre", "Copy"), class_="code-block-head"),
                           ui.tags.pre(ui.output_text("cite_bibtex", inline=True), id="cite_bibtex_pre",
                                       class_="code-block"),
                           class_="code-block-wrap"),
                    ui.div(ui.download_button("dl_bib", "Download .bib", class_="btn-sm btn-outline-primary"),
                           ui.download_button("dl_ris", "Download .ris (EndNote, Zotero, Mendeley)",
                                              class_="btn-sm btn-outline-primary"),
                           class_="dl-row")),
            section("Suggested acknowledgement",
                    ui.div(ack, id="cite_ack", class_="quote-box"),
                    ui.div(copy_button("cite_ack", "Copy acknowledgement"), class_="citation-actions")),
            section("Citing the methods",
                    ui.p("Where a particular construction is central to your work, please also cite the original "
                         "publication listed for it on the Methods and About tabs.")),
            class_="info-page"),
        value="cite")


def legal_page():
    return ui.nav_panel(
        "Copyright & disclaimer",
        ui.div(
            page_head("Copyright & disclaimer", kicker="Legal",
                      lead="Terms under which this software and its content are made available."),
            section("Copyright",
                    ui.div(ui.tags.strong(T.copyright_line()), class_="notice-box"),
                    ui.p(ui.HTML(f"The software &ldquo;<em>{html.escape(T.TITLE)}</em>&rdquo; ({T.SHORT_NAME}) was "
                                 f"developed by Sukanta Dash at the {T.INSTITUTE_FULL} ({T.INSTITUTE_SHORT}), "
                                 f"{T.PARENT_ORG}, New Delhi.")),
                    ui.p("The application, its source code, user interface, documentation and associated content are "
                         "protected under the Copyright Act, 1957 (India) and applicable international law."),
                    ui.tags.ul(
                        ui.tags.li("Matrices and designs generated with this software may be used freely for research, "
                                   "teaching and extension, including in publications, theses and reports, provided the "
                                   "software is cited as described under How to cite."),
                        ui.tags.li("You may print or download extracts for personal and non-commercial use."),
                        ui.tags.li("Except with the express written permission of the copyright holder, the software and "
                                   "its content may not be redistributed, re-hosted, modified or commercially exploited."),
                        class_="legal-list")),
            section("Disclaimer",
                    ui.tags.ol(
                        ui.tags.li("The software is provided “as is” for research and educational use, without "
                                   "warranty of any kind, express or implied, including warranties of accuracy, "
                                   "completeness, reliability or fitness for a particular purpose."),
                        ui.tags.li("Every matrix shown is checked against its defining identity (H·Hᵀ = nI) and "
                                   "every orthogonal array against its strength; for very large s the check is sampled and "
                                   "labelled as such. Users remain responsible for verifying a design before laying out "
                                   "an experiment or publishing results."),
                        ui.tags.li(f"In no event shall the author, {T.INSTITUTE_SHORT} or the {T.PARENT_ORG} be liable "
                                   "for any loss or damage, including indirect or consequential loss, arising from the use "
                                   "of, or inability to use, this software or its outputs."),
                        ui.tags.li("Links to other websites are provided for convenience only and do not imply endorsement."),
                        ui.tags.li("Every effort is made to keep the service available; no responsibility is accepted for "
                                   "temporary unavailability due to technical reasons beyond our control."),
                        class_="legal-list")),
            section("Data privacy",
                    ui.p("The application does not ask for, store or share any personal information. All matrices and "
                         "designs are computed within your session and are discarded when the session ends.")),
            section("Third-party software",
                    ui.p("This application is written in Python and uses Shiny for Python, NumPy, pandas and Matplotlib; "
                         "formulas are typeset with MathJax. Each package remains the property of its authors and is "
                         "distributed under its own licence.")),
            class_="info-page"),
        value="legal")


def contact_page():
    return ui.nav_panel(
        "Contact",
        ui.div(
            page_head("Contact", kicker="Developer",
                      lead="For further information, queries about the constructions or suggestions for improvement, "
                           "please contact the developer."),
            ui.div(
                ui.div(ui.div(
                    ui.div(ui.span(T.DEVELOPER), ui.span("Developer", class_="badge-corr"), class_="contact-name"),
                    ui.div(T.DESIGNATION, class_="contact-role"),
                    ui.div(T.DIVISION, class_="contact-aff"),
                    ui.div(T.INSTITUTE_FULL, class_="contact-aff"),
                    ui.div("Library Avenue, Pusa, New Delhi 110012", class_="contact-aff"),
                    ui.div(icon("phone"), ui.span(T.PHONE), class_="contact-line"),
                    ui.div(icon("mail"), ui.tags.a(T.EMAIL, href=f"mailto:{T.EMAIL}"), class_="contact-line"),
                    class_="contact-body"), class_="contact-card is-corresponding"),
                class_="contact-grid"),
            ui.layout_columns(
                section("Institute",
                        ui.div(icon("pin"), ui.div(ui.tags.strong(T.INSTITUTE_FULL), ui.br(), T.PARENT_ORG, ui.br(),
                                                   T.ADDRESS), class_="inst-line"),
                        ui.div(icon("phone"), ui.span(f"Phone: {T.PHONE_OFFICE}; Fax: {T.FAX}"), class_="inst-line"),
                        ui.div(icon("globe"), ui.tags.a("iasri.icar.gov.in", href=T.INSTITUTE_URL, target="_blank",
                                                        rel="noopener"), class_="inst-line")),
                section("Feedback and support",
                        ui.p("When reporting a problem, please mention the tab, the order or parameters entered and "
                             "the output received."),
                        ui.div(ui.tags.a("Report an error", href=mailto("error"), class_="btn btn-primary btn-sm"),
                               ui.tags.a("Comments / suggestions", href=mailto("comments"),
                                         class_="btn btn-outline-primary btn-sm"),
                               class_="d-flex flex-wrap gap-2")),
                col_widths=[7, 5], class_="mt-3"),
            class_="info-page"),
        value="contact")


def cite_panel():
    return ui.div(
        ui.div(ui.h3("Cite this software", class_="info-section-title"),
               ui.input_action_link("home_to_cite", "Full citation, BibTeX and RIS →", class_="cite-link"),
               class_="cite-panel-head"),
        ui.div(ui.output_ui("cite_home", inline=True), class_="citation-text"),
        class_="info-section cite-panel")


def inst_strip():
    return ui.div(ui.div(ui.span(T.INSTITUTE_HINDI, class_="inst-hi", lang="hi"), ui.span("|", class_="inst-sep"),
                         ui.span(T.INSTITUTE_FULL, class_="inst-en"), ui.span(T.PARENT_ORG, class_="inst-parent"),
                         class_="inst-strip-inner"),
                  class_="inst-strip", id="inst-strip")


def footer(version: str):
    return ui.tags.footer(
        ui.div(
            ui.div(ui.div(T.SHORT_NAME, class_="footer-title"),
                   ui.div("Hadamard matrices, MOLS, orthogonal arrays and lattice designs", class_="footer-text"),
                   ui.div(f"Version {version} · Developed by Sukanta Dash", class_="footer-text muted"),
                   class_="footer-col"),
            ui.div(ui.div("Information", class_="footer-head"),
                   ui.input_action_link("ft_cite", "How to cite"),
                   ui.input_action_link("ft_copyright", "Copyright"),
                   ui.input_action_link("ft_disclaimer", "Disclaimer"),
                   ui.input_action_link("ft_contact", "Contact"), class_="footer-col"),
            ui.div(ui.div("Feedback", class_="footer-head"),
                   ui.tags.a("Report error", href=mailto("error")),
                   ui.tags.a("Comments / suggestions", href=mailto("comments")), class_="footer-col"),
            ui.div(ui.div(T.INSTITUTE_SHORT, class_="footer-head"),
                   ui.div(T.INSTITUTE_FULL, class_="footer-text"), ui.div(T.ADDRESS, class_="footer-text"),
                   ui.div(f"Phone: {T.PHONE_OFFICE} | Fax: {T.FAX}", class_="footer-text muted"),
                   class_="footer-col"),
            class_="footer-inner"),
        ui.div(ui.span(T.copyright_line(), class_="footer-copy"),
               ui.span("Provided for research and educational use, without warranty."), class_="footer-bottom"),
        class_="app-footer")


# --------------------------------------------------------------------------
# scripts: move the institutional strip above the navbar, copy buttons, and
# link every "Sukanta Dash" to his biodata page (standing instruction)
# --------------------------------------------------------------------------
SCRIPT = """
(function () {
  var SD_URL = '%s';
  var RX = /(Dr\\.?\\s+)?Sukanta Dash/g;
  var SKIP = {A:1, BUTTON:1, PRE:1, CODE:1, INPUT:1, TEXTAREA:1, SELECT:1, SCRIPT:1, STYLE:1, OPTION:1};
  function skip(node) {
    for (var el = node.parentNode; el && el.nodeType === 1; el = el.parentNode) {
      if (SKIP[el.tagName] || (el.classList && el.classList.contains('no-sd-link'))) return true;
    }
    return false;
  }
  function linkify(root) {
    if (!root) return;
    var w = document.createTreeWalker(root, NodeFilter.SHOW_TEXT, null), n, nodes = [];
    while ((n = w.nextNode())) { if (n.nodeValue.indexOf('Sukanta Dash') > -1 && !skip(n)) nodes.push(n); }
    nodes.forEach(function (t) {
      var frag = document.createDocumentFragment(), s = t.nodeValue, last = 0, m;
      RX.lastIndex = 0;
      while ((m = RX.exec(s))) {
        frag.appendChild(document.createTextNode(s.slice(last, m.index)));
        var a = document.createElement('a');
        a.href = SD_URL; a.target = '_blank'; a.rel = 'noopener'; a.className = 'sd-link';
        a.textContent = m[0]; frag.appendChild(a); last = m.index + m[0].length;
      }
      frag.appendChild(document.createTextNode(s.slice(last)));
      t.parentNode.replaceChild(frag, t);
    });
  }
  function copyText(text) {
    if (navigator.clipboard && window.isSecureContext) return navigator.clipboard.writeText(text);
    return new Promise(function (ok, bad) {
      var ta = document.createElement('textarea'); ta.value = text; ta.style.position = 'fixed'; ta.style.opacity = '0';
      document.body.appendChild(ta); ta.select();
      try { document.execCommand('copy') ? ok() : bad(); } catch (e) { bad(e); }
      document.body.removeChild(ta);
    });
  }
  document.addEventListener('click', function (ev) {
    var b = ev.target.closest('.copy-btn'); if (!b) return;
    var el = document.getElementById(b.getAttribute('data-copy-target')); if (!el) return;
    var lab = b.querySelector('.copy-label'), orig = lab ? lab.textContent : '';
    copyText((el.innerText || el.textContent || '').trim()).then(function () {
      b.classList.add('copied'); if (lab) lab.textContent = 'Copied';
      setTimeout(function () { b.classList.remove('copied'); if (lab) lab.textContent = orig; }, 1800);
    });
  });
  function start() {
    var strip = document.getElementById('inst-strip');
    if (strip && document.body.firstChild !== strip) document.body.insertBefore(strip, document.body.firstChild);
    linkify(document.body);
    new MutationObserver(function (ms) {
      ms.forEach(function (m) { m.addedNodes.forEach(function (x) {
        if (x.nodeType === 1) linkify(x); else if (x.nodeType === 3 && x.parentNode) linkify(x.parentNode);
      }); });
    }).observe(document.body, {childList: true, subtree: true});
  }
  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', start); else start();
})();
""" % T.SD_URL

FONTS = ui.TagList(
    ui.tags.link(rel="preconnect", href="https://fonts.googleapis.com"),
    ui.tags.link(rel="preconnect", href="https://fonts.gstatic.com", crossorigin=""),
    ui.tags.link(rel="stylesheet", href="https://fonts.googleapis.com/css2?family=Source+Sans+3:wght@400;500;600;700"
                                        "&family=Source+Serif+4:opsz,wght@8..60,600;8..60,700&display=swap"),
)
