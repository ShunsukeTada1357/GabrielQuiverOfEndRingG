###############################################################################
# Test script for `VisualizeGquiverOfEnd`
#
# Input:
#   A "poset matrix" P (transitive closure), i.e.
#     P[i, j] = 1  ⇔  there exists a directed path i → ... → j in the Hasse graph,
#     P[i, j] = 0  otherwise.
#   (In particular, P has 1's on the diagonal.)
#
# Example:
#   If the underlying graph is 1 → 2 → 3, then the transitive-closure matrix is
#     P = [1 1 1;
#          0 1 1;
#          0 0 1]
###############################################################################

include("Gquiver.jl")
using .Gquiver

###############################################################################
# Example 1: Total order on 5 elements (A5)
###############################################################################
A5 = [1 1 1 1 1;
      0 1 1 1 1;
      0 0 1 1 1;
      0 0 0 1 1;
      0 0 0 0 1]

# Visualize the Gabriel quiver of End(G)
Gquiver.VisualizeGquiverOfEnd(A5)

# Compute the adjacency matrix ARQ and list all directed cycles (if any)
ARQ, AllInt = Gquiver.GquiverOfEndRing(A5)
cycles = Gquiver.all_directed_cycles(ARQ)

# Save the Gabriel quiver using Graphviz (DOT/PNG, etc.)
Gquiver.output(A5, "A5")


###############################################################################
# Example 2: A-tilde type poset on 6 elements (expected to contain cycles)
###############################################################################
Atilda = [1 1 0 0 1 1;
          0 1 0 0 0 0;
          0 1 1 0 0 0;
          0 1 1 1 1 0;
          0 0 0 0 1 0;
          0 0 0 0 1 1]

# Visualize the Gabriel quiver of End(G)
Gquiver.VisualizeGquiverOfEnd(Atilda)

# Detect a directed cycle and (optionally) extract one cycle
ok, cyc = Gquiver.Isdirected_with_cycle(Atilda)

# Compute ARQ and list all directed cycles (if any)
ARQ, AllInt = Gquiver.GquiverOfEndRing(Atilda)
cycles = Gquiver.all_directed_cycles(ARQ)

# Save the Gabriel quiver using Graphviz (DOT/PNG, etc.)
Gquiver.output(Atilda, "Atilda")


###############################################################################
# Other basic posets (transitive-closure matrices)
###############################################################################
ComGrid25 = [
    1 1 1 1 1 1 1 1 1 1;
    0 1 1 1 1 0 1 1 1 1;
    0 0 1 1 1 0 0 1 1 1;
    0 0 0 1 1 0 0 0 1 1;
    0 0 0 0 1 0 0 0 0 1;
    0 0 0 0 0 1 1 1 1 1;
    0 0 0 0 0 0 1 1 1 1;
    0 0 0 0 0 0 0 1 1 1;
    0 0 0 0 0 0 0 0 1 1;
    0 0 0 0 0 0 0 0 0 1
]

ComGrid33 = [
    1 1 1 1 1 1 1 1 1;
    0 1 1 0 1 1 0 1 1;
    0 0 1 0 0 1 0 0 1;
    0 0 0 1 1 1 1 1 1;
    0 0 0 0 1 1 0 1 1;
    0 0 0 0 0 1 0 0 1;
    0 0 0 0 0 0 1 1 1;
    0 0 0 0 0 0 0 1 1;
    0 0 0 0 0 0 0 0 1
]

B12 = [1 1 1 1 1;
       0 1 1 0 0;
       0 0 1 0 0;
       0 0 1 1 1;
       0 0 1 0 1]

Atilda4  = [1 1 0 1;
            0 1 0 0;
            0 1 1 1;
            0 0 0 1]

Atilda5  = [1 1 0 0 1;
            0 1 0 0 0;
            0 1 1 0 0;
            0 1 1 1 1;
            0 0 0 0 1]

Atilda6  = [1 1 0 0 0 1;
            0 1 0 0 0 0;
            0 1 1 1 0 0;
            0 0 0 1 0 0;
            0 0 0 1 1 1;
            0 0 0 0 0 1]

Atilda10 = [1 1 0 0 0 0 0 0 0 1;
            0 1 0 0 0 0 0 0 0 0;
            0 1 1 1 0 0 0 0 0 0;
            0 0 0 1 0 0 0 0 0 0;
            0 0 0 1 1 1 0 0 0 0;
            0 0 0 0 0 1 0 0 0 0;
            0 0 0 0 0 1 1 1 0 0;
            0 0 0 0 0 0 0 1 0 0;
            0 0 0 0 0 0 0 1 1 1;
            0 0 0 0 0 0 0 0 0 1]
