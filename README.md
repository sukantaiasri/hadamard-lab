# HadamardLab

A Shiny (Python) application for constructing and **verifying** Hadamard matrices of order ≤ 1000,
complete sets of mutually orthogonal Latin squares (MOLS), orthogonal arrays OA(s², s+1, s, 2) and
square-lattice designs.

**Version 2.0 (November 2026)** adds seven constructions (Golay pairs, Golay-derived and stored T-sequences,
Legendre pairs, Legendre quadruples, cyclic Hadamard difference sets, Scarpis–Đoković), a *Construction lab*
tab that builds matrices from their ingredient sequences, 24 recent references, and the data behind the orders
668, 716 and 892 (announced by Alpöge, Voinov and Reynolds-Haertle on 12 August 2026 and analysed by Egan,
Kotsireas, London and Ó Catháin, 2026) and 876. **Every admissible order ≤ 1000 is now constructible.**

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
| Hadamard generator | 20 constructions: Sylvester, Kronecker, Paley I/II, generalised Paley II, Williamson, Goethals–Seidel with Turyn/base/Golay/stored T-sequences, Baumert–Hall, Ehlich, Miyamoto, Cooper–Wallis, supplementary difference sets, Golay two-circulant, Legendre pairs (two-circulant core), Legendre quadruples (bordered Goethals–Seidel), cyclic Hadamard difference sets (Singer, twin prime power, Hall), Scarpis–Đoković. Normalisation, three notations, H·Hᵀ = nI verification, construction path, properties, Plackett–Burman design, CSV/TXT/JSON export. |
| Construction lab (new) | Build a matrix from its ingredients: a Golay pair, T-sequences (optionally with Williamson matrices → Cooper–Wallis), a Legendre pair, a Legendre quadruple, a difference set or a Scarpis input matrix. Shows the sequences, the autocorrelation identity they satisfy (with plot), the sign pattern, verification and exports (matrix, ingredients as JSON/text). |
| Existence explorer | Map and searchable table of all orders 4 … 1000 with every applicable method. All 250 orders are constructible; no order ≤ 1000 is open (smallest open order: 2396). |
| MOLS & orthogonal arrays | Complete set of s − 1 MOLS over GF(s) for every prime power s ≤ 997, with the orthogonal array. Verification is exhaustive for small s and sampled (and labelled so) for large s. |
| Lattice designs | Simple / triple / balanced square lattices for k² treatments, field book, concurrence structure, efficiency factors. |
| Methods, About | Theorem statements (new methods flagged), reachable orders per method, provenance, implementation notes, recent references (2013–2026) and the full bibliography. |
| How to cite, Copyright & disclaimer, Contact | Citation with live access date, BibTeX/RIS downloads, terms of use, developer and institute details. |

## Layout

```
app.py                  Shiny UI + server (thin wrapper)
hadamard_lab/
  gf.py                 GF(p^k) arithmetic, quadratic character, Jacobsthal matrices
  constructions.py      the 20 constructions and their applicability rules
  sequences.py          Golay pairs, two-level sequences / difference sets, Legendre pairs and quadruples,
                        two-circulant and bordered Goethals–Seidel arrays, Scarpis–Đoković (new)
  seqlab.py             back-end of the Construction lab tab (new)
  planner.py            feasibility DP over orders, recursive construction, traces
  analysis.py           verification, normalisation, properties, export, Plackett–Burman
  mols.py               MOLS, orthogonal arrays, lattice designs
  content.py            method statements, corrections log, bibliography; title, URL, citation, copyright (edit here)
  site.py               How to cite, Copyright & disclaimer and Contact pages, institute strip, footer
  data/sequences.json   Williamson / base / Turyn sequences, SDS (incl. N = 223), T-sequences (t = 73),
                        Legendre pairs and quadruples (incl. m = 166, 178) – all verified on import
tests/                  179 tests
www/style.css
```

## Tests

```bash
pip install pytest && python -m pytest -q
```

The suite verifies the bundled sequence data against its defining equations, every construction,
**every constructible order with every applicable method** (H·Hᵀ = nI), the finite fields of all prime
powers ≤ 1000, MOLS / OA strength, and lattice efficiency (balanced lattices must give k/(k+1)).

## New constructions in version 2.0

| Method | Order | Reference |
|---|---|---|
| Golay pair, two circulants `[[A, B], [-Bᵀ, Aᵀ]]` | 2g, g = 2^a·10^b·26^c | Golay 1961; Turyn 1974; Borwein & Ferguson 2004 |
| Goethals–Seidel from Golay-derived T-sequences | 4g and 4(2g+1) | Turyn 1974 |
| Goethals–Seidel from stored T-sequences (t = 73) | 4t | Cooper & Wallis 1972; Egan et al. 2026 |
| Legendre pair → two-circulant core | 2ℓ + 2 | Fletcher, Gysin & Seberry 2001; Kotsireas & Koutschan 2021 |
| Cyclic Hadamard difference set → circulant / group-developed core | ℓ + 1 | Singer 1938; Stanton & Sprott 1958; Hall 1956; Whiteman 1962 |
| Legendre quadruple → bordered Goethals–Seidel (orders 668, 716) | 4m + 4 | Egan, Kotsireas, London & Ó Catháin 2026 (Thm 5); Alpöge, Voinov & Reynolds-Haertle 2026 |
| Scarpis–Đoković: H(q+1) → H(q(q+1)), q ≡ 3 (mod 4) | q(q+1) | Scarpis 1898; Đoković 2017 |

The sequences for orders 668, 716 (Legendre quadruples of lengths 166, 178), 892 (supplementary difference
sets in Z₂₂₃) and 876 (T-sequences of length 73) were recovered from the publicly available matrices of the
Hadamard Matrix Explorer tables (https://www.math74.com/Hadamard/, Egan et al. 2026) and are re-verified by
the tests. Legendre pairs of lengths 9, 21, 33, 39 and quadruples of even length m ≤ 44 were found by a small
tabu search written for this release (not shipped with the app).

## Sources and implementation notes

Constructions are implemented from the published literature (see the References on the About tab) and
checked numerically. Sequence tables follow those published by C. Koukouvinos (Williamson, base and
Turyn-type sequences) and D. Ž. Djoković (supplementary difference sets); every stored row is
re-verified by the tests. Four formulas fail when transcribed literally and were re-derived
(Ehlich sign, Miyamoto sign table, Cooper–Wallis last term, generalised Paley II congruence); see the
About tab and `hadamard_lab/content.py`. Order 764 is generated from Djoković's supplementary
difference sets for N = 191.

## Limitations

* Orders are capped at 1000. New Williamson/base/Turyn/SDS/T-sequence/Legendre data can be added to
  `data/sequences.json` (the tests re-verify any new rows).
* Finite-field element labelling for prime powers p^k, k > 1, depends on the primitive polynomial
  chosen; resulting matrices are equivalent to, but not necessarily identical with, those of other software.
