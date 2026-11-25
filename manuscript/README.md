# Granvik's real-part conjecture in OEIS A292273

The note proves Granvik's conjecture in A292273, revision 24, FORMULA, 6 July 2024, for every positive integer. It proves the adjacent parity formulas, support and density, generalized finite trajectory sums, Dirichlet inverse and convolution square-root identities. It corrects the odd-prime-power formula and two reciprocal constants in A087003. The terminal-1 interpretation applies whenever that terminal value is reached; universal Collatz termination is outside the statements.

## Files

- `note.tex`: complete mathematical text and bibliography.
- `note.pdf`: typeset note, distributed separately from the source archive.
- `metadata.json`: note identity, authors, license and mathematical evidence link.
- `verify.py`: standard-library Python 3 finite integer checks.
- `verification.json`: deterministic output of those finite checks.
- `source-values.json`: displayed sequence values with exact scientific source locators.
- `LICENSE.txt`: license for the note and accompanying original code.

## Build

Use Tectonic 0.17.0:

```sh
tectonic note.tex
```

Alternatively, with a LaTeX installation containing `geometry`, `amsmath`, `amssymb`, `amsthm`, `booktabs`, `array`, `longtable`, `fontenc`, `lmodern`, `xurl`, and `hyperref`:

```sh
pdflatex -interaction=nonstopmode -halt-on-error note.tex
pdflatex -interaction=nonstopmode -halt-on-error note.tex
```

The bibliography is contained in `note.tex`; no BibTeX step or external scientific data file is required for the PDF. Tectonic obtains its standard typesetting resources on the first build.

## Finite checks

```sh
python3 verify.py > verification-new.json
diff -u verification.json verification-new.json
```

The bounds and counts are included in the output. The implementation compares a sieve-based Möbius function with trial factorization, checks the four residue classes and sequence identities through 65537, compares displayed values with the scientific sources, and checks the finite Collatz identity for starting values 1 through 1000. It checks divisor convolutions through 2048, the rational square root through 512, generalized odd branches at finite prefixes, prime-parameter squarefree counts, and the source corrections by exact arithmetic. The constant checks use rational intervals with explicit tails, not floating-point comparisons. These are finite checks; the note supplies the all-index and asymptotic proofs.

## Formal statement

The declaration `D5.S3.Arith.Mobius.GranvikMobiusRealPart.result` in the [pinned Lean file](https://github.com/the-omega-institute/trureturing/blob/3bedf2870aacc6e7f17353444aadaa49fb73157c/D5/S3/Arith/Mobius/GranvikMobiusRealPart.lean) proves the real-part equality for every positive natural index using the standard integer-valued Möbius function. Its companion `a` is the parity-based definition. The note's general real-valued criterion and additional corollaries are written proofs, not further declarations in that file.

## Attribution and license

Haobo Ma (ChronoAI Pte Ltd, Singapore), followed by Wenlin Zhang (National University of Singapore; ChronoAI Pte Ltd, Singapore). CC BY 4.0 (`cc-by-4.0`). The source attribution identifies Granvik's 2024 conjecture, Karttunen's 2017 entry, LeBrun's 2004 cancellation argument, Michon's 2007 inverse relation, Jameson's 2010 density proof, and the dated later OEIS formulas. The coverage table distinguishes proved results, corrected formulas, prior proofs and statements outside the coverage.
