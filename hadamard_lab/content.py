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
    "golay2c": dict(
        statement=r"""A Golay complementary pair is a pair \(a, b\) of \(\pm1\) sequences of length \(g\) whose
        non-periodic autocorrelations cancel: \(N_a(s)+N_b(s)=0\) for every \(s\ge1\). Periodic autocorrelations then
        cancel as well, so the circulants satisfy \(AA^{\top}+BB^{\top}=2gI\) and
        \(H=\begin{pmatrix}A&B\\-B^{\top}&A^{\top}\end{pmatrix}\) is Hadamard of order \(2g\). Pairs exist for every
        \(g=2^a10^b26^c\): the seeds of lengths 2, 10 and 26 are combined by Turyn's product
        \(x=a\otimes\tfrac{c+d}{2}+\tilde b\otimes\tfrac{c-d}{2},\; y=b\otimes\tfrac{c+d}{2}-\tilde a\otimes\tfrac{c-d}{2}\)
        (\(\tilde{}\) = reversal). No other length \(g\le105\) admits a pair (Borwein & Ferguson 2004; the smallest
        undecided length is 106).""",
        condition=r"\(n=2g\), \(g=2^a10^b26^c\) (23 orders \(\le1000\))",
        notes="The seeds were verified by exhaustive search; Turyn's product is checked for every combination of seeds in the test-suite."),
    "gs_golay": dict(
        statement=r"""A Golay pair \((a,b)\) of length \(g\) gives T-sequences in two ways: directly,
        \(\big(\tfrac{a+b}{2},\tfrac{a-b}{2},0,0\big)\) of length \(g\), and through Turyn's base sequences
        \(BS(g+1,g)=\big((a,1),(a,-1),b,b\big)\), which give T-sequences of length \(2g+1\). The four circulants
        \(A_1=T_1{+}T_2{+}T_3{+}T_4,\dots\) (as for base sequences) in the Goethals–Seidel array give orders
        \(4g\) and \(4(2g+1)\); the same T-sequences also feed the Cooper–Wallis construction.""",
        condition=r"\(n=4g\) or \(n=4(2g+1)\) with \(g=2^a10^b26^c\)",
        notes=""),
    "gs_tseq": dict(
        statement=r"""T-sequences stored directly (rather than derived from base or Turyn-type sequences) are
        placed in the Goethals–Seidel array exactly as above. The bundled set of length 73 was recovered from the
        order-876 matrix \(G(T(73),W(3))\) of the Hadamard Matrix Explorer tables (Egan, Kotsireas, London &
        Ó Catháin 2026); it is verified to have disjoint supports and zero total non-periodic autocorrelation,
        gives order 292 directly and, with Williamson matrices of order 3, order 876 by Cooper–Wallis.""",
        condition=r"\(n=4t\) with stored T-sequences of length \(t\) (\(t=73\))",
        notes="Order 876 was listed as ‘exists in the literature, not implemented’ in version 1; it is now built from these T-sequences."),
    "legendre": dict(
        statement=r"""A Legendre pair is a pair \(a,b\) of \(\pm1\) sequences of odd length \(\ell\) with
        \(P_a(s)+P_b(s)=-2\) for all \(s\ne0\) (periodic autocorrelation), equivalently
        \(AA^{\top}+BB^{\top}=(2\ell+2)I-2J\) for the circulants, normalised to row sum \(+1\). The two-circulant core
        with a \(2\times2\) border (Fletcher, Gysin & Seberry 2001)
        \(\begin{pmatrix}-&-&j&j\\-&+&j&-j\\j^{\top}&j^{\top}&A&B\\j^{\top}&-j^{\top}&B^{\top}&-A^{\top}\end{pmatrix}\)
        is Hadamard of order \(2\ell+2\). Sources used here: the Legendre symbol sequences modulo a prime \(\ell\)
        (two copies when \(\ell\equiv3\), the choices \(a_0=\pm1\) when \(\ell\equiv1 \pmod 4\)); two copies of any
        sequence with two-level periodic autocorrelation (m-sequences, twin primes, Hall sextic residues); and pairs of
        composite length found by computer search. It is conjectured that a Legendre pair exists for every odd
        \(\ell\) (Đoković & Kotsireas 2015); this would prove the Hadamard conjecture for all \(n\equiv0 \pmod 4\).""",
        condition=r"\(n=2\ell+2\), \(\ell\) odd: prime, \(2^k-1\), twin primes \(q(q+2)\), Hall prime, or a stored pair",
        notes=("Recent lengths: 77 (Turner, Kotsireas, Bulutoglu & Geyer 2021), 91, 93, 123 (Balonin & Đoković 2021), "
               "85 and 87 (Kotsireas, Koutschan, Bulutoglu et al. 2023); lengths \\(\\equiv0 \\pmod 3\\) by Kotsireas & Koutschan (2021). "
               "The pairs of lengths 9, 21, 33 and 39 bundled here were found by HadamardLab's own tabu search.")),
    "dscore": dict(
        statement=r"""Let \(D\) be a Hadamard difference set in an abelian group \(G\) of order \(\ell\equiv3 \pmod 4\)
        (parameters \((\ell,\tfrac{\ell-1}{2},\tfrac{\ell-3}{4})\)) and \(f=+1\) on \(D\), \(-1\) elsewhere, signed so
        that \(\sum f=-1\). The group-developed matrix \(M=(f(y-x))\) satisfies \(MM^{\top}=(\ell+1)I-J\) and
        \(MJ=-J\), so the bordered matrix \(\begin{pmatrix}1&e^{\top}\\e&M\end{pmatrix}\) is Hadamard of order
        \(\ell+1\) (circulant core when \(G\) is cyclic). Besides Paley's quadratic residues (= Paley I) the classical
        families are: Singer sets / m-sequences, \(\ell=2^k-1\) (maximal-length shift-register sequences from a
        primitive polynomial over GF(2)); twin prime powers, \(\ell=q(q+2)\) with
        \(D=\{(x,y):\chi_1(x)\chi_2(y)=1\}\cup\{(x,0)\}\) (Stanton & Sprott 1958; Whiteman 1962); and Hall's sextic
        residues, \(\ell=4x^2+27\) prime, \(D=C_0\cup C_1\cup C_3\) (Hall 1956).""",
        condition=r"\(n=\ell+1\), \(\ell\in\{2^k-1\}\cup\{q(q+2)\}\cup\{4x^2+27\ \text{prime}\}\)",
        notes=("For twin prime powers that are not both prime (7·9, 9·11, 23·25, 25·27, 27·29) the group "
               "GF(q)⁺ × GF(q+2)⁺ is not cyclic and the core is group-developed rather than circulant; the cyclic cases "
               "(3·5, 5·7, 11·13, 17·19, 29·31) also give Legendre pairs. Hall's union of classes is selected by testing the "
               "twenty candidate unions.")),
    "lquad": dict(
        statement=r"""Four \(\pm1\) sequences of even length \(m\) whose periodic autocorrelations sum to \(-4\) at
        every non-zero shift (a <em>Legendre quadruple</em>, Egan, Kotsireas, London & Ó Catháin 2026) have row sums
        \((\pm2,0,0,0)\); with circulants normalised to row sums \(2,0,0,0\),
        \(AA^{\top}+BB^{\top}+CC^{\top}+DD^{\top}=(4m+4)I-4J\), and the Goethals–Seidel array with a border of
        width four,
        \(H=\begin{pmatrix}K&E\\F&M\end{pmatrix}\) with \(K=2I_4-J_4\), \(E,F\) signed blocks of ones and \(M\) the
        Goethals–Seidel array, is Hadamard of order \(4m+4\) (their Theorem 5). The matrices of orders 668 and 716
        announced by Alpöge, Voinov and Reynolds-Haertle on 12 August 2026 — the two smallest orders that had been
        open since 2005 — are of exactly this type, with \(m=166\) and \(m=178\); their sequences are bundled.""",
        condition=r"\(n=4m+4\), \(m\) even with a stored quadruple (\(m\le44\), 166, 178)",
        notes=("Alpöge, Voinov and Reynolds-Haertle published the matrices but not the search method; Egan et al. (2026) "
               "identified the structure. Quadruples of length \\(m\\le44\\) were found by HadamardLab's own search.")),
    "scarpis": dict(
        statement=r"""Let \(q\equiv3 \pmod 4\) be a prime power and \(A\) any Hadamard matrix of order \(n=q+1\),
        normalised so that its first row and column are \(+1\); let \(C\) be its core with rows \(c(x)\) indexed by
        \(x\in GF(q)\), so that \(c(x)\cdot c(y)=-1\) for \(x\ne y\) and every row sum is \(-1\). Put
        \(B_0=A'\otimes j\) (\(A'\) = \(A\) without its first row) and, for \(r\in GF(q)\),
        \(B_r=\big[\,j^{\top}\otimes c(r)\;\big|\;B_{r,x}\ (x\in GF(q))\,\big]\) where row \(k\) of \(B_{r,x}\) is
        \(c(xr+k)\). Then \(B=(B_0;B_1;\dots;B_q)\) is Hadamard of order \(qn=q(q+1)\). Scarpis (1898) proved this
        for primes; Đoković (2017) observed that only the field structure of the row labels is used, so any prime
        power \(q\) and any input matrix of order \(q+1\) work.""",
        condition=r"\(n=q(q+1)\), \(q\in\{3,7,11,19,23,27,31\}\) (orders 12, 56, 132, 380, 552, 756, 992)",
        notes="The input matrix of order q + 1 is built by the planner (Paley I by default); q = 27 gives order 756, the smallest case beyond Scarpis's original theorem."),
    "sds": dict(
        statement=r"""Let \(G=\mathbb Z_N\) and \(S_1,\dots,S_4\subset G\) be supplementary difference sets with
        \(\sum|S_i|=N+\lambda\) and every non-zero element represented \(\lambda\) times in total. With
        \(A_i(x,y)=-1\) if \(y-x\in S_i\) and \(+1\) otherwise, the Goethals–Seidel array has order \(4N\). Djoković's
        sets in \(\mathbb Z_{103},\mathbb Z_{127},\mathbb Z_{151},\mathbb Z_{163},\mathbb Z_{181},\mathbb Z_{191},\mathbb Z_{239}\)
        give orders 412, 508, 604, 652, 724, 764 and 956; the sets in \(\mathbb Z_{223}\) (sizes 106, 106, 106, 100,
        \(\lambda=195\)) recovered from the order-892 matrix of Alpöge, Voinov and Reynolds-Haertle (2026) give order 892.""",
        condition=r"\(n=4N\), \(N\in\{103,127,151,163,181,191,223,239\}\)",
        notes="Order 764, open for many years, follows from Đoković's sets for \\(N=191\\) (2008); order 892 was open until August 2026."),
}

# --------------------------------------------------------------------------
# Implementation notes on formulas that fail when transcribed literally (About tab)
# --------------------------------------------------------------------------
CORRECTIONS = [
    ("Ehlich", "Sign of the term \\(I_m\\otimes(J_k-I_k)\\) in \\(K\\) must be negative."),
    ("Miyamoto", "Border and block sign pattern replaced by a computationally verified table."),
    ("Cooper–Wallis", "Last term of \\(W\\) is \\(+T_4\\otimes A\\)."),
    ("Paley II (generalised)", "Condition is \\(q\\equiv1 \\pmod 4\\), not \\(p^\\alpha\\equiv2 \\pmod 4\\)."),
    ("Coverage (v2.0)", "Every admissible order ≤ 1000 can now be generated. Orders 668, 716 and 892 — the smallest open cases of the "
                        "Hadamard conjecture from 2005 until August 2026 — are built from the sequences of Alpöge, Voinov and "
                        "Reynolds-Haertle (2026) as analysed by Egan, Kotsireas, London and Ó Catháin (2026); order 876 from the "
                        "T-sequences of length 73 (Cooper–Wallis with Williamson order 3). The smallest order with no known Hadamard "
                        "matrix is now 2396."),
    ("Legendre quadruples / pairs", "Searched sequences (m ≤ 44, ℓ ∈ {9, 21, 33, 39}) are checked against their autocorrelation "
                                    "identities on import; the two-circulant and bordered arrays are verified for every stored length."),
]

# --------------------------------------------------------------------------
# Bibliography
# --------------------------------------------------------------------------
REFERENCES = [
    "Alpöge, L., Voinov, P. and Reynolds-Haertle, S. (2026). Hadamard matrices of orders 668, 716, 892, 1132, 1244, 1388, 1436, 1676, 1772, 1916, 1948 and 1964. Announcement, 12 August 2026; see <a href='https://epoch.ai/frontiermath/open-problems/hadamard' target='_blank' rel='noopener'>epoch.ai/frontiermath/open-problems/hadamard</a>.",
    "Balonin, N. A. and Đoković, D. Ž. (2021). Three new lengths for cyclic Legendre pairs. <i>Information and Control Systems</i> 1, 2–7.",
    "Baumert, L. D. and Hall, M. Jr. (1965). A new construction for Hadamard matrices. <i>Bull. Amer. Math. Soc.</i> 71, 169–170.",
    "Best, D., Đoković, D. Ž., Kharaghani, H. and Ramp, H. (2013). Turyn-type sequences: classification, enumeration, and construction. <i>J. Combin. Des.</i> 21, 24–35.",
    "Borwein, P. and Ferguson, R. A. (2004). A complete description of Golay pairs for lengths up to 100. <i>Math. Comp.</i> 73, 967–985.",
    "Bose, R. C., Shrikhande, S. S. and Parker, E. T. (1960). Further results on the construction of mutually orthogonal Latin squares and the falsity of Euler's conjecture. <i>Canad. J. Math.</i> 12, 189–203.",
    "Bright, C. (2018). A doubling construction for Williamson matrices. arXiv:1803.01480.",
    "Bright, C., Kotsireas, I. and Ganesh, V. (2020). Applying computer algebra systems with SAT solvers to the Williamson conjecture. <i>J. Symbolic Comput.</i> 100, 187–209.",
    "Bulutoglu, D., Baczkowski, D. and Yauney, J. (2026). Determining the group that sends each Legendre pair to an equivalent Legendre pair. arXiv:2604.22423.",
    "Bush, K. A. (1952). Orthogonal arrays of index unity. <i>Ann. Math. Statist.</i> 23, 426–434.",
    "Cati, M. and Pasechnik, D. V. (2024). A database of constructions of Hadamard matrices (implemented in SageMath). arXiv:2411.18897.",
    "Cooper, J. and Wallis, J. (1972). A construction for Hadamard arrays. <i>Bull. Austral. Math. Soc.</i> 7, 269–278.",
    "Craigen, R. (1996). Hadamard matrices and designs. In: <i>CRC Handbook of Combinatorial Designs</i> (Colbourn and Dinitz, eds.), CRC Press, 370–377.",
    "Craigen, R. and Kharaghani, H. (2007). Hadamard matrices and Hadamard designs. In: <i>Handbook of Combinatorial Designs</i>, 2nd ed. (Colbourn and Dinitz, eds.), Chapman &amp; Hall/CRC, 273–280.",
    "Đoković, D. Ž. (1992–1994). Constructions of Hadamard matrices of Goethals–Seidel type (several papers in <i>Bull. Austral. Math. Soc.</i>, <i>J. Combin. Theory A</i>, <i>Europ. J. Combin.</i>, <i>Combinatorica</i> and <i>Austral. J. Combin.</i>).",
    "Đoković, D. Ž. (2008). Hadamard matrices of order 764 exist. <i>Combinatorica</i> 28, 487–489.",
    "Đoković, D. Ž. (2011). Small orders of Hadamard matrices and base sequences. <i>Int. Math. Forum</i> 6(62), 3061–3067 (arXiv:1008.2043).",
    "Đoković, D. Ž. (2017). Generalization of Scarpis's theorem on Hadamard matrices. <i>Linear Multilinear Algebra</i> 65(10) (arXiv:1601.00635).",
    "Đoković, D. Ž. (2023). Hadamard matrices: skew of orders 276, 292 and symmetric 372. arXiv:2301.02751.",
    "Đoković, D. Ž., Golubitsky, O. and Kotsireas, I. S. (2014). Some new orders of Hadamard and skew-Hadamard matrices. <i>J. Combin. Des.</i> 22, 270–277.",
    "Đoković, D. Ž. and Kotsireas, I. S. (2015). Compression of periodic complementary sequences and applications. <i>Des. Codes Cryptogr.</i> 74, 365–377.",
    "Egan, R., Kotsireas, I. S., London, S. and Ó Catháin, P. (2026). Recent progress in the construction of Hadamard matrices. arXiv:2610.27175; matrices at <a href='https://www.math74.com/Hadamard/' target='_blank' rel='noopener'>math74.com/Hadamard</a>.",
    "Ehlich, H. (1965). Neue Hadamard-Matrizen. <i>Arch. Math.</i> 16, 34–36.",
    "Fletcher, R. J., Gysin, M. and Seberry, J. (2001). Application of the discrete Fourier transform to the search for generalised Legendre pairs and Hadamard matrices. <i>Australas. J. Combin.</i> 23, 75–86.",
    "Goethals, J. M. and Seidel, J. J. (1967). Orthogonal matrices with zero diagonal. <i>Canad. J. Math.</i> 19, 1001–1010.",
    "Goethals, J. M. and Seidel, J. J. (1970). A skew Hadamard matrix of order 36. <i>J. Austral. Math. Soc.</i> 11, 343–344.",
    "Golay, M. J. E. (1961). Complementary series. <i>IRE Trans. Inform. Theory</i> IT-7, 82–87.",
    "Hall, M. Jr. (1956). A survey of difference sets. <i>Proc. Amer. Math. Soc.</i> 7, 975–986.",
    "Hall, M. Jr. (1986). <i>Combinatorial Theory</i>, 2nd ed. Wiley.",
    "Hedayat, A., Sloane, N. J. A. and Stufken, J. (1999). <i>Orthogonal Arrays: Theory and Applications</i>. Springer.",
    "Hedayat, A. and Wallis, W. D. (1978). Hadamard matrices and their applications. <i>Ann. Statist.</i> 6, 1184–1238.",
    "Holzmann, W. H., Kharaghani, H. and Tayfeh-Rezaie, B. (2008). Williamson matrices up to order 59. <i>Des. Codes Cryptogr.</i> 46, 343–352.",
    "Horadam, K. J. (2007). <i>Hadamard Matrices and Their Applications</i>. Princeton University Press.",
    "Kharaghani, H., Mohammadian, A. and Tayfeh-Rezaie, B. (2026). A search for Hadamard matrices of Williamson type. arXiv:2605.08661.",
    "Kharaghani, H. and Tayfeh-Rezaie, B. (2005). A Hadamard matrix of order 428. <i>J. Combin. Designs</i> 13, 435–440.",
    "Kotsireas, I. and Koutschan, C. (2021). Legendre pairs of lengths ℓ ≡ 0 (mod 3). <i>J. Combin. Des.</i> 29(12), 870–887.",
    "Kotsireas, I., Koutschan, C., Bulutoglu, D., Arquette, D., Turner, J. and Ryan, K. (2023). Legendre pairs of lengths ℓ ≡ 0 (mod 5). <i>Special Matrices</i> 11, 20230105.",
    "Kotsireas, I. S., Koutschan, C. and Winterhof, A. (2025). Quaternary Legendre pairs II. <i>Discrete Math.</i> 348(9), 114501.",
    "Leung, K. H., Momihara, K. and Xiang, Q. (2019). A new family of Hadamard matrices of order 4(2q² + 1). arXiv:1907.02623.",
    "London, S. and Kotsireas, I. (2025). New Turyn-type sequences. <i>Cryptogr. Commun.</i> 17(5), 1601–1610.",
    "Lumsden, T., Kotsireas, I. and Bright, C. (2025). New results on periodic Golay pairs. <i>Math. Comp.</i>, doi:10.1090/mcom/4096.",
    "Miyamoto, M. (1991). A construction of Hadamard matrices. <i>J. Combin. Theory Ser. A</i> 57, 86–108.",
    "Paley, R. E. A. C. (1933). On orthogonal matrices. <i>J. Math. Phys.</i> 12, 311–320.",
    "Plackett, R. L. and Burman, J. P. (1946). The design of optimum multifactorial experiments. <i>Biometrika</i> 33, 305–325.",
    "Ramos, A. F., Hulak, D. B. and de Queiroz, R. J. G. B. (2026). Multiplier obstructions for Legendre pairs of length 333. <i>J. Combin. Des.</i>, to appear (arXiv:2607.20765).",
    "Scarpis, U. (1898). Sui determinanti di valore massimo. <i>Rend. R. Ist. Lombardo Sci. Lett.</i> (2) 31, 1441–1446.",
    "Seberry, J. and Balonin, N. A. (2015). The propus construction for symmetric Hadamard matrices. arXiv:1512.01732.",
    "Seberry, J. and Yamada, M. (1992). Hadamard matrices, sequences and block designs. In: <i>Contemporary Design Theory</i> (Dinitz and Stinson, eds.), Wiley, 431–560.",
    "Seberry, J. and Yamada, M. (2020). <i>Hadamard Matrices: Constructions using Number Theory and Linear Algebra</i>. Wiley.",
    "Singer, J. (1938). A theorem in finite projective geometry and some applications to number theory. <i>Trans. Amer. Math. Soc.</i> 43, 377–385.",
    "Stanton, R. G. and Sprott, D. A. (1958). A family of difference sets. <i>Canad. J. Math.</i> 10, 73–77.",
    "Sylvester, J. J. (1867). Thoughts on inverse orthogonal matrices, simultaneous sign successions, and tessellated pavements. <i>Phil. Mag.</i> 34, 461–475.",
    "Turner, J., Kotsireas, I., Bulutoglu, D. and Geyer, A. (2021). Legendre pair of length 77 using complementary binary matrices with fixed marginals. <i>Des. Codes Cryptogr.</i> (arXiv:2101.10918).",
    "Turyn, R. J. (1972). An infinite class of Williamson matrices. <i>J. Combin. Theory Ser. A</i> 12, 319–321.",
    "Turyn, R. J. (1974). Hadamard matrices, Baumert–Hall units, four-symbol sequences, pulse compression, and surface wave encodings. <i>J. Combin. Theory Ser. A</i> 16, 313–333.",
    "Wallis, J. and Whiteman, A. L. (1972). Some classes of Hadamard matrices with constant diagonal. <i>Bull. Austral. Math. Soc.</i> 7, 233–249.",
    "Whiteman, A. L. (1962). A family of difference sets. <i>Illinois J. Math.</i> 6, 107–121.",
    "Williamson, G., Yacobi, O. and Zinn-Justin, P. (2026). Generating Hadamard matrices with transformers. arXiv:2604.11101.",
    "Williamson, J. (1944). Hadamard's determinant theorem and the sum of four squares. <i>Duke Math. J.</i> 11, 65–81.",
    "Yates, F. (1936). A new method of arranging variety trials involving a large number of varieties. <i>J. Agric. Sci.</i> 26, 424–455.",
]

# References that are new in version 2.0 (2013–2026), listed separately on the About tab
RECENT_REFERENCE_KEYS = ("Alpöge", "Balonin", "Best, D.", "Bright, C. (2018)", "Bright, C., Kotsireas", "Bulutoglu", "Cati",
                         "Đoković, D. Ž. (2017)", "Đoković, D. Ž. (2023)", "Đoković, D. Ž., Golubitsky", "Đoković, D. Ž. and Kotsireas",
                         "Egan", "Kharaghani, H., Mohammadian", "Kotsireas, I. and Koutschan", "Kotsireas, I., Koutschan, C., Bulutoglu",
                         "Kotsireas, I. S., Koutschan, C. and Winterhof", "Leung", "London", "Lumsden", "Ramos", "Seberry, J. and Balonin",
                         "Seberry, J. and Yamada, M. (2020)", "Turner", "Williamson, G.")


def recent_references():
    return [r for r in REFERENCES if any(r.startswith(k) for k in RECENT_REFERENCE_KEYS)]

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
    "the supplementary difference sets for N ≤ 191 and 239 are due to D. Ž. Đoković. The sets in Z<sub>223</sub> "
    "(order 892), the Legendre quadruples of lengths 166 and 178 (orders 668, 716) and the T-sequences of length 73 "
    "(orders 292, 876) were recovered from the publicly available matrices of the Hadamard Matrix Explorer tables "
    "(Egan, Kotsireas, London and Ó Catháin, 2026); the first three are the matrices announced by Alpöge, Voinov and "
    "Reynolds-Haertle (2026). Golay pairs are generated from the seeds of lengths 2, 10 and 26; the remaining Legendre "
    "pairs and quadruples were found by HadamardLab's own search. Every stored row is re-verified against its "
    "defining equations by the test-suite, and Galois-field polynomials are found by direct search."
)
