# HadamardLab

A Shiny (Python) application for constructing and **verifying** Hadamard matrices of order ≤ 1000,
complete sets of mutually orthogonal Latin squares (MOLS), orthogonal arrays OA(s², s+1, s, 2) and
square-lattice designs.

**Developed by Dr. Sukanta Dash**, Senior Scientist, Division of Statistical Ecology and Environmental
Statistics, ICAR-IASRI, New Delhi-110012 · sukantadash.iasri@icar.org.in

**How to cite** (the access date is filled in automatically on the How to cite page):
Dash, S. (2026). HadamardLab: A Verified Generator of Hadamard Matrices, Mutually Orthogonal Latin Squares,
Orthogonal Arrays and Lattice Designs. https://iasri-sukanta-hadamardlab.share.connect.posit.cloud/ (accessed on <date>).

© 2026 ICAR-Indian Agricultural Statistics Research Institute, New Delhi. All rights reserved.

## Run

```bash
pip install -r requirements.txt
shiny run --reload app.py            # http://127.0.0.1:8000
# or simply:  python app.py          (starts the server and opens the browser)
```

Docker: `docker build -t hadamardlab . && docker run -p 8000:8000 hadamardlab`.
Posit Connect Cloud / shinyapps.io: `pip install rsconnect-python && rsconnect deploy shiny .`
(deploy without the local `.venv/` folder; `requirements.txt` lists the packages).

## What it does

| Tab | Content |
|---|---|
| Hadamard generator | 13 constructions (Sylvester, Kronecker, Paley I/II, generalised Paley II, Williamson, Goethals–Seidel with Turyn/base sequences, Baumert–Hall, Ehlich, Miyamoto, Cooper–Wallis, supplementary difference sets). Normalisation, three notations, H·Hᵀ = nI verification, construction path, properties, Plackett–Burman design, CSV/TXT/JSON export. |
| Existence explorer | Map and searchable table of all orders 4 … 1000 with every applicable method. 246 of 250 orders are constructible; 668, 716, 892 are open cases and 876 needs data that is not bundled. |
| MOLS & orthogonal arrays | Complete set of s − 1 MOLS over GF(s) for every prime power s ≤ 997, with the orthogonal array. Verification is exhaustive for small s and sampled (and labelled so) for large s. |
| Lattice designs | Simple / triple / balanced square lattices for k² treatments, field book, concurrence structure, efficiency factors. |
| Methods, About | Theorem statements, reachable orders per method, provenance, corrections, references. |
| How to cite, Copyright & disclaimer, Contact | Citation with live access date, BibTeX/RIS downloads, terms of use, developer and institute details. |

## Layout

```
app.py                  Shiny UI + server (thin wrapper)
hadamard_lab/
  gf.py                 GF(p^k) arithmetic, quadratic character, Jacobsthal matrices
  constructions.py      the 13 constructions and their applicability rules
  planner.py            feasibility DP over orders, recursive construction, traces
  analysis.py           verification, normalisation, properties, export, Plackett–Burman
  mols.py               MOLS, orthogonal arrays, lattice designs
  content.py            method statements, corrections log, bibliography; title, URL, citation, copyright (edit here)
  site.py               How to cite, Copyright & disclaimer and Contact pages, institute strip, footer
  data/sequences.json   Williamson / base / Turyn sequences and SDS (verified)
tests/                  115 tests
www/style.css
```

## Tests

```bash
pip install pytest && python -m pytest -q
```

The suite verifies the bundled sequence data against its defining equations, every construction,
**every constructible order with every applicable method** (H·Hᵀ = nI), the finite fields of all prime
powers ≤ 1000, MOLS / OA strength, and lattice efficiency (balanced lattices must give k/(k+1)).

## Sources and implementation notes

Constructions are implemented from the published literature (see the References on the About tab) and
checked numerically. Sequence tables follow those published by C. Koukouvinos (Williamson, base and
Turyn-type sequences) and D. Ž. Djoković (supplementary difference sets); every stored row is
re-verified by the tests. Four formulas fail when transcribed literally and were re-derived
(Ehlich sign, Miyamoto sign table, Cooper–Wallis last term, generalised Paley II congruence); see the
About tab and `hadamard_lab/content.py`. Order 764 is generated from Djoković's supplementary
difference sets for N = 191.

## Limitations

* Orders are capped at 1000. New Williamson/base/Turyn/SDS data can be added to `data/sequences.json`
  (the tests re-verify any new rows).
* Finite-field element labelling for prime powers p^k, k > 1, depends on the primitive polynomial
  chosen; resulting matrices are equivalent to, but not necessarily identical with, those of other software.
