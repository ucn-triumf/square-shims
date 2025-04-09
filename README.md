# square-shims

To install:

1.  git clone this repo
2.  cd square-shims
3.  git clone https://github.com/jmartin1454/patchlib.git
4.  git clone https://github.com/jmartin1454/Pis

Then run particular-coils-onesheet.py.

What the code does:

- sets up an array of coils inscribed on a cube in the design geometry
  of the TUCAN MSR in early 2025 (the "fake shims")

- sets up an array of fluxgates distributed throughout a smaller cube

- calculates the matrix in B=M*I by setting each coil current I in
  turn and measuring B at each sensor position.  This calculation is
  currently done assuming free space using the patch library.  A
  future upgrade will be to use COMSOL to do the calculation.

- inverts the matrix.  Carefully excises the zero mode (if desired).

- selects a target field B_target (uses the Pis library)

- calculates I_set = M^{-1} B_target and sets those currents on the
  cube

- makes graphs of the resultant true field and graphs them in a few
  views (presently along axes), comparing with the target field.
  Again, this is presently done in free space, using the patch
  library.  A future upgrade will be to use a linear combination of
  COMSOL results to do the calculation.

- a lot of the different functionalities are implemented through
  command-line switches.  Run the program using the --help switch for
  more information.

