Heegaard builder
----------------

This repository contains source code for taking a one-vertex triangulation of
a 3-manifold with boundary, and filling the boundary with a handlebody to
obtain a one-vertex triangulation of a closed 3-manifold.

The filling is specified by endowing the boundary of the input triangulation
with an object called a filling bouquet, which gives a combinatorial encoding
of a system of attaching circles. In the special case where the input
triangulation is a *handlebody*, the system of attaching circles effectively
describes a Heegaard splitting, and hence the filling algorithm gives a way
to turn a Heegaard splitting of a closed 3-manifold *M* into a one-vertex
triangulation of *M*.

The filling algorithm was designed in joint work with *James Morgan* and *Em
Thompson*. Details can be found in our preprint at:
    https://arxiv.org/abs/2312.17556.

The main implementation for this algorithm is contained in the script
``heegaardbuilder.py``. The best ways to use this script are to either:
- run it with Regina's regina-python command-line program in interactive
    mode; or
- import it into the Python console provided by Regina's graphical user
    interface.

See ``heegaardExample.txt`` for a transcript of a regina-python session
using ``heegaardbuilder.py``.

For details on how to install Regina, visit:
    https://regina-normal.github.io/.

— *Alex He (a.he@uqconnect.edu.au)*
