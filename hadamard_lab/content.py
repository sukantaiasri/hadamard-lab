"""Static scientific content (statements, notes, bibliography).

Mathematics is written for MathJax: inline ``\\( ... \\)``.  The prose is original;
the theorem statements are standard results from the cited literature.
"""

# --------------------------------------------------------------------------
# Method documentation
# --------------------------------------------------------------------------
METHOD_DOCS = {
    "sylvester": dict(
        statement=r"""If \(H_a\) and \(H_b\) are Hadamard, so is \(H_a\otimes H_b\). Starting from
        \(H_2=\begin{pmatrix}1&1\\1&-1\end{pmatrix}\), the \(k\)-fold Kronecker power gives an Hadamard
        matrix of order \(2^k\) (the Walsh–Hadamard matrix).""",
        condition=r"\(n = 2^k\)",
        notes=""),
    "kronecker": dict(
        statement=r"""The Kronecker product of Hadamard matrices of orders \(a\) and \(b\) is an Hadamard matrix of order
        \(ab\), because \((A\otimes B)(A\otimes B)^{\top}=AA^{\top}\otimes BB^{\top}=ab\,I\).""",
        condition=r"\(n=ab\) with both orders constructible",
        notes="The planner chooses the factorisation with the shallowest recursion."),
    "paley1": dict(
        statement=r"""Let \(q\equiv 3 \pmod 4\) be a prime power, \(\chi\) the quadratic character of \(GF(q)\) and
        \(Q=(\chi(x-y))\) the Jacobsthal matrix. With \(S=\begin{pmatrix}0&e^{\top}\\-e&Q\end{pmatrix}\) the
        matrix \(H=S+I\) is Hadamard of order \(q+1\) and of skew type (\(H+H^{\top}=2I\)).""",
        condition=r"\(n=q+1\), \(q\equiv3 \pmod 4\) a prime power",
        notes=""),
    "paley2": dict(
        statement=r"""Let \(q\equiv 1 \pmod 4\) be a prime power and \(S=\begin{pmatrix}0&e^{\top}\\e&Q\end{pmatrix}\)
        the (symmetric) conference matrix of order \(q+1\). Then
        \(H=S\otimes\begin{pmatrix}1&1\\1&-1\end{pmatrix}+I\otimes\begin{pmatrix}1&-1\\-1&-1\end{pmatrix}\)
        is Hadamard of order \(2(q+1)\).""",
        condition=r"\(n=2(q+1)\), \(q\equiv1 \pmod 4\) a prime power",
        notes=""),
    "paley2_general": dict(
        statement=r"""Williamson's generalisation: if \(A\) is Hadamard of order \(h\) (even), \(U=I_{h/2}\otimes
        \begin{pmatrix}0&1\\-1&0\end{pmatrix}\), \(B=UA\) and \(S\) is the symmetric conference matrix of order \(q+1\)
        (\(q\equiv 1 \pmod 4\)), then \(H=A\otimes I_{q+1}+B\otimes S\) is Hadamard of order \(h(q+1)\).""",
        condition=r"\(n=h(q+1)\), \(h\ge4\) an Hadamard order, \(q\equiv1 \pmod 4\) a prime power",
        notes=("The theorem is sometimes quoted with the condition “\\(p^{\\alpha}\\equiv 2 \\pmod 4\\)”, which cannot "
               "be met by an odd prime power; the theorem (and this implementation) uses \\(q\\equiv1 \\pmod 4\\).")),
    "williamson": dict(
        statement=r"""If \(A,B,C,D\) are symmetric circulant \(\pm1\) matrices of order \(w\) with
        \(A^2+B^2+C^2+D^2=4wI\), the Williamson array
        \(\begin{pmatrix}A&B&C&D\\-B&A&-D&C\\-C&D&A&-B\\-D&-C&B&A\end{pmatrix}\) is Hadamard of order \(4w\).""",
        condition=r"\(n=4w\) with stored Williamson matrices of order \(w\) (16 orders up to 63)",
        notes="Initial rows originate from tables published by C. Koukouvinos."),
    "gs_base": dict(
        statement=r"""Base sequences \((A,B,C,D)\) of lengths \(m+p,m+p,m,m\) with vanishing total aperiodic
        autocorrelation yield four T-sequences of length \(2m+p\):
        \(((A{+}B)/2,0),\,((A{-}B)/2,0),\,(0,(C{+}D)/2),\,(0,(C{-}D)/2)\). Circulants \(T_i\) of these give
        \(A_1=T_1{+}T_2{+}T_3{+}T_4\), \(A_2=-T_1{+}T_2{+}T_3{-}T_4\), \(A_3=-T_1{-}T_2{+}T_3{+}T_4\),
        \(A_4=-T_1{+}T_2{-}T_3{+}T_4\), which are inserted in the Goethals–Seidel array
        \(\begin{pmatrix}A&BR&CR&DR\\-BR&A&D^{\top}R&-C^{\top}R\\-CR&-D^{\top}R&A&B^{\top}R\\-DR&C^{\top}R&-B^{\top}R&A\end{pmatrix}\)
        (\(R\) the back-diagonal identity). The result has order \(4(2m+p)\).""",
        condition=r"\(n=4(2m+p)\) with stored base sequences",
        notes=""),
    "gs_turyn": dict(
        statement=r"""Turyn-type sequences \((X,Y,Z,W)\) of lengths \(t,t,t,t-1\) with
        \(N_X+N_Y+2N_Z+2N_W=0\) give base sequences \((Z{,}W),(Z{,}-W),X,Y\) of lengths
        \(2t-1,2t-1,t,t\) and hence T-sequences of length \(3t-1\). The Goethals–Seidel array then has order
        \(4(3t-1)\). Order 428 was first obtained this way (Kharaghani and Tayfeh-Rezaie, 2005).""",
        condition=r"\(n=4(3t-1)\), \(t\in\{30,34,36\}\) (orders 356, 404, 428)",
        notes=""),
    "baumert_hall": dict(
        statement=r"""The Baumert–Hall array is a \(12\times12\) array in \(\pm A,\pm B,\pm C,\pm D\) (an orthogonal design
        \(OD(12;3,3,3,3)\)). Replacing the symbols by commuting \(\pm1\) matrices of order \(w\) with
        \(AA^{\top}+BB^{\top}+CC^{\top}+DD^{\top}=4wI\) (e.g. Williamson matrices) gives order \(12w\).""",
        condition=r"\(n=12w\) with stored Williamson matrices of order \(w\) (or \(w=1\))",
        notes="The array is verified in the test-suite to be an exact orthogonal design OD(12; 3,3,3,3)."),
    "ehlich": dict(
        statement=r"""Let \(k\equiv3 \pmod 4\) and \(m=k-2\equiv1 \pmod 4\) be prime powers, \(A\) the Jacobsthal matrix of
        \(GF(m)\) (so \(AA^{\top}=mI-J\)) and \(G\) the Jacobsthal matrix of \(GF(k)\) (skew, from Paley I). Then
        \(H=\begin{pmatrix}1&e^{\top}\\e&K\end{pmatrix}\) with
        \(K=A\otimes G-I_m\otimes(J_k-I_k)+J_m\otimes I_k\) is Hadamard of order \((k-1)^2\).""",
        condition=r"\(n=(k-1)^2\), \(k\in\{7,11,19,27,31\}\) (orders 36, 100, 324, 676, 900)",
        notes=("Implementation note: the formula is sometimes transcribed with \\(+I_m\\otimes(J_k-I_k)\\). With that sign the matrix is not "
               "orthogonal (the row sums of \\(K\\) must be \\(-1\\)); the sign above is required and is checked by the test-suite.")),
    "miyamoto": dict(
        statement=r"""For a prime power \(q\equiv1 \pmod 4\) and an Hadamard matrix \(K\) of order \(q-1\), Miyamoto (1991)
        builds from the Paley conference matrix (blocks \(C_1,C_2,C_4\) of size \(m=(q-1)/2\)) and the four
        quadrants of \(K\) two \(4m\times4m\) block matrices \(U\) and \(V\). The blocks
        \(T_{ij}=U_{ij}\otimes J_2+V_{ij}\otimes\begin{pmatrix}1&-1\\-1&1\end{pmatrix}\), bordered by signed rows and
        columns of ones, assemble an Hadamard matrix of order \(4q\).""",
        condition=r"\(n=4q\), \(q\equiv1 \pmod 4\) a prime power, Hadamard matrix of order \(q-1\) available",
        notes=("Implementation note: the border/sign pattern as commonly transcribed does not give an orthogonal matrix. "
               "The block skeleton is unchanged; the sign table in the source (<code>MIYAMOTO_SIGNS</code>) was found by exhaustive "
               "search at \\(q=5\\) and verified for every admissible \\(q\\le 249\\) with several choices of \\(K\\).")),
    "cooper_wallis": dict(
        statement=r"""With Williamson matrices \(A,B,C,D\) of order \(w\) and circulants \(T_1,\dots,T_4\) of T-sequences of
        length \(L\), put \(X=\sum T_i\otimes M_i\) with \((M_i)=(A,B,C,D)\), \(Y=(-B,A,D,-C)\), \(Z=(-C,-D,A,B)\),
        \(W=(-D,C,-B,A)\) (same pairing). Using \(X,Y,Z,W\) in the Goethals–Seidel array gives order \(4wL\).""",
        condition=r"\(n=4wL\), stored Williamson order \(w\ge3\), T-sequence length \(L\in\{7,11,19,47,\dots,107\}\)",
        notes=("Implementation note: a transcription with \\(-T_4\\otimes A\\) in \\(W\\) fails. The coefficient array must have orthogonal "
               "<em>columns</em>, which forces \\(+T_4\\otimes A\\).")),
    "sds": dict(
        statement=r"""Let \(G=\mathbb Z_N\) and \(S_1,\dots,S_4\subset G\) be supplementary difference sets with
        \(\sum|S_i|=N+\lambda\) and every non-zero element represented \(\lambda\) times in total. With
        \(A_i(x,y)=-1\) if \(y-x\in S_i\) and \(+1\) otherwise, the Goethals–Seidel array has order \(4N\). Djoković's
        sets in \(\mathbb Z_{103},\mathbb Z_{127},\mathbb Z_{151},\mathbb Z_{163},\mathbb Z_{181},\mathbb Z_{191},\mathbb Z_{239}\)
        give orders 412, 508, 604, 652, 724, 764 and 956.""",
        condition=r"\(n=4N\), \(N\in\{103,127,151,163,181,191,239\}\)",
        notes="Order 764, open for many years, follows from Djoković's sets for \\(N=191\\) (2008)."),
}

# --------------------------------------------------------------------------
# Implementation notes on formulas that fail when transcribed literally (About tab)
# --------------------------------------------------------------------------
CORRECTIONS = [
    ("Ehlich", "Sign of the term \\(I_m\\otimes(J_k-I_k)\\) in \\(K\\) must be negative."),
    ("Miyamoto", "Border and block sign pattern replaced by a computationally verified table."),
    ("Cooper–Wallis", "Last term of \\(W\\) is \\(+T_4\\otimes A\\)."),
    ("Paley II (generalised)", "Condition is \\(q\\equiv1 \\pmod 4\\), not \\(p^\\alpha\\equiv2 \\pmod 4\\)."),
    ("Coverage", "Orders 764 (supplementary difference sets) and 964 (Miyamoto) can be generated; only 668, 716, 892 (open) and 876 are unavailable."),
]

# --------------------------------------------------------------------------
# Bibliography
# --------------------------------------------------------------------------
REFERENCES = [
    "Baumert, L. D. and Hall, M. Jr. (1965). A new construction for Hadamard matrices. <i>Bull. Amer. Math. Soc.</i> 71, 169–170.",
    "Bose, R. C., Shrikhande, S. S. and Parker, E. T. (1960). Further results on the construction of mutually orthogonal Latin squares and the falsity of Euler's conjecture. <i>Canad. J. Math.</i> 12, 189–203.",
    "Bush, K. A. (1952). Orthogonal arrays of index unity. <i>Ann. Math. Statist.</i> 23, 426–434.",
    "Cooper, J. and Wallis, J. (1972). A construction for Hadamard arrays. <i>Bull. Austral. Math. Soc.</i> 7, 269–278.",
    "Craigen, R. (1996). Hadamard matrices and designs. In: <i>CRC Handbook of Combinatorial Designs</i> (Colbourn and Dinitz, eds.), CRC Press, 370–377.",
    "Djoković, D. Ž. (1992–1994). Constructions of Hadamard matrices of Goethals–Seidel type (several papers in <i>Bull. Austral. Math. Soc.</i>, <i>J. Combin. Theory A</i>, <i>Europ. J. Combin.</i>, <i>Combinatorica</i> and <i>Austral. J. Combin.</i>).",
    "Djoković, D. Ž. (2008). Hadamard matrices of order 764 exist. <i>Combinatorica</i> 28, 487–489.",
    "Ehlich, H. (1965). Neue Hadamard-Matrizen. <i>Arch. Math.</i> 16, 34–36.",
    "Goethals, J. M. and Seidel, J. J. (1967). Orthogonal matrices with zero diagonal. <i>Canad. J. Math.</i> 19, 1001–1010.",
    "Hall, M. Jr. (1986). <i>Combinatorial Theory</i>, 2nd ed. Wiley.",
    "Hedayat, A. and Wallis, W. D. (1978). Hadamard matrices and their applications. <i>Ann. Statist.</i> 6, 1184–1238.",
    "Kharaghani, H. and Tayfeh-Rezaie, B. (2005). A Hadamard matrix of order 428. <i>J. Combin. Designs</i> 13, 435–440.",
    "Miyamoto, M. (1991). A construction of Hadamard matrices. <i>J. Combin. Theory Ser. A</i> 57, 86–108.",
    "Paley, R. E. A. C. (1933). On orthogonal matrices. <i>J. Math. Phys.</i> 12, 311–320.",
    "Plackett, R. L. and Burman, J. P. (1946). The design of optimum multifactorial experiments. <i>Biometrika</i> 33, 305–325.",
    "Seberry, J. and Yamada, M. (1992). Hadamard matrices, sequences and block designs. In: <i>Contemporary Design Theory</i> (Dinitz and Stinson, eds.), Wiley, 431–560.",
    "Sylvester, J. J. (1867). Thoughts on inverse orthogonal matrices, simultaneous sign successions, and tessellated pavements. <i>Phil. Mag.</i> 34, 461–475.",
    "Wallis, J. and Whiteman, A. L. (1972). Some classes of Hadamard matrices with constant diagonal. <i>Bull. Austral. Math. Soc.</i> 7, 233–249.",
    "Williamson, J. (1944). Hadamard's determinant theorem and the sum of four squares. <i>Duke Math. J.</i> 11, 65–81.",
    "Yates, F. (1936). A new method of arranging variety trials involving a large number of varieties. <i>J. Agric. Sci.</i> 26, 424–455.",
]

DEVELOPER = "Dr. Sukanta Dash"
DESIGNATION = "Senior Scientist"
DIVISION = "Division of Statistical Ecology and Environmental Statistics"
INSTITUTE = "ICAR-IASRI, New Delhi-110012"
EMAIL = "sukantadash.iasri@icar.org.in"
PHONE = "011-25847127-4158"
AFFILIATION_LINE = f"{DESIGNATION}, {DIVISION}, {INSTITUTE}"

# --------------------------------------------------------------------------
# Software identity, citation and copyright (edit here)
# --------------------------------------------------------------------------
SHORT_NAME = "HadamardLab"
TITLE = ("HadamardLab: A Verified Generator of Hadamard Matrices, Mutually Orthogonal Latin Squares, "
         "Orthogonal Arrays and Lattice Designs")
YEAR = 2026
URL = "https://iasri-sukanta-hadamardlab.share.connect.posit.cloud/"
AUTHORS = [{"cite": "Dash, S.", "bib": "Dash, S."}]
INSTITUTE_FULL = "ICAR-Indian Agricultural Statistics Research Institute"
INSTITUTE_SHORT = "ICAR-IASRI"
INSTITUTE_HINDI = "\u092d\u093e\u0915\u0943\u0905\u0928\u0941\u092a-\u092d\u093e\u0930\u0924\u0940\u092f \u0915\u0943\u0937\u093f \u0938\u093e\u0902\u0916\u094d\u092f\u093f\u0915\u0940 \u0905\u0928\u0941\u0938\u0902\u0927\u093e\u0928 \u0938\u0902\u0938\u094d\u0925\u093e\u0928"
PARENT_ORG = "Indian Council of Agricultural Research"
ADDRESS = "Library Avenue, Pusa, New Delhi 110012, India"
INSTITUTE_URL = "https://iasri.icar.gov.in/"
PHONE_OFFICE = "91-11-25847121-24, 25841254 (PBX)"
FAX = "91-11-25841564"
COPYRIGHT_HOLDER = "ICAR-Indian Agricultural Statistics Research Institute, New Delhi"
SD_URL = "https://farmerfirst.halforfarmers.com/sukanta-dash-biodata.html"   # standing instruction


def access_date():
    """Today's date in India (UTC+05:30), used as the access date of the citation."""
    import datetime as _d
    return _d.datetime.now(_d.timezone(_d.timedelta(hours=5, minutes=30))).date()


def long_date(d=None) -> str:
    d = d or access_date()
    return f"{d.day} {d.strftime('%B')} {d.year}"


def citation_authors() -> str:
    x = [a["cite"] for a in AUTHORS]
    if len(x) == 1:
        return x[0]
    if len(x) == 2:
        return f"{x[0]} and {x[1]}"
    return ", ".join(x[:-1]) + ", and " + x[-1]


def citation_text(d=None) -> str:
    return f"{citation_authors()} ({YEAR}). {TITLE}. {URL} (accessed on {long_date(d)})."


def citation_html(d=None) -> str:
    import html
    return (f"{html.escape(citation_authors())} ({YEAR}). <em>{html.escape(TITLE)}</em>. "
            f"<a href='{URL}' target='_blank' rel='noopener'>{URL}</a> (accessed on {long_date(d)}).")


def citation_bibtex(d=None) -> str:
    authors = " and ".join(a["bib"] for a in AUTHORS)
    key = AUTHORS[0]["cite"].split(",")[0].lower() + str(YEAR) + SHORT_NAME.lower()
    return ("@misc{" + key + ",\n"
            f"  author       = {{{authors}}},\n"
            f"  title        = {{{TITLE}}},\n"
            f"  year         = {{{YEAR}}},\n"
            f"  howpublished = {{\\url{{{URL}}}}},\n"
            f"  note         = {{Accessed on {long_date(d)}}}\n"
            "}")


def citation_ris(d=None) -> str:
    d = d or access_date()
    lines = ["TY  - COMP"] + [f"AU  - {a['cite']}" for a in AUTHORS] + [
        f"TI  - {TITLE}", f"PY  - {YEAR}", f"UR  - {URL}", f"Y2  - {d.strftime('%Y/%m/%d')}",
        f"PB  - {INSTITUTE_FULL}, New Delhi", "ER  - ", ""]
    return "\r\n".join(lines)


def copyright_line() -> str:
    return f"\u00a9 {YEAR} {COPYRIGHT_HOLDER}. All rights reserved."


DEVELOPER_CREDIT = (
    f"HadamardLab was developed by <b>{DEVELOPER}</b>, {DESIGNATION}, {DIVISION}, {INSTITUTE}. "
    f"Contact: <a href='mailto:{EMAIL}'>{EMAIL}</a>. All constructions are implemented from the published "
    "literature listed below and verified numerically; no third-party generator code is used."
)

DATA_CREDIT = (
    "Initial rows for Williamson, base and Turyn-type sequences follow tables published by C. Koukouvinos; "
    "the supplementary difference sets are due to D. Ž. Djoković. Every stored row is re-verified against its "
    "defining equations by the test-suite, and Galois-field polynomials are found by direct search."
)
