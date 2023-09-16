##This is a test code for the function VisualizeGquiverOfEnd
## input is "Poset Matrix", 
## i.e., if  there is a path from i to j on the graph, then the (i,j) component of the "Poset Matrix" is 1, otherwise 0.
## For example if the graph is 1→2→3, then
#poset=[1 1 1;0 1 1; 0 0 1 ]    
using GraphRecipes
using Plots,PyPlot
using PyCall
using Combinatorics
using LightGraphs, GraphPlot, LinearAlgebra
using AbstractAlgebra, Nemo
import Base.Iterators: flatten
using .GabrieQuiverOfEndG: VisualizeGquiverOfEnd


Pos = [1 0 0; 1 1 1; 0 0 1 ]
ZZ = AbstractAlgebra.ZZ
Atype=[1 1 1; 0 1 1; 0 0 1]
A4=[1 1 1 1 ; 0 1 1 1 ; 0 0 1 1; 0 0 0 1]
C12=[1 1 1 1 1; 0 1 1 0 0 ; 0 0 1 0 0 ; 0 0 1 1 1; 0 0 1 0 1]
A5=[1 1 1 1 1; 0 1 1 1 1; 0 0 1 1 1; 0 0 0 1 1; 0 0 0 0 1]
ComGrid=
[1 1 1 1 1 1
0 1 1 0 1 1 
0 0 1 0 0 1
0 0 0 1 1 1 
0 0 0 0 1 1 
0 0 0 0 0 1]

Atilda4=[1 1 0 1; 0 1 0 0; 0 1 1 1; 0 0 0 1]
zig3=[1 0 0; 1 1 1 ; 0 0 1]
Atilda6=[1 1 0 0 0 1 ; 0 1 0 0 0 0 ; 0 1 1 1 0 0; 0 0 0 1 0 0 ;0 0 0 1 1 1; 0 0 0 0 0 1]
Atilda10=[1 1 0 0 0 0 0 0 0 1; 0 1 0 0 0 0 0 0 0 0; 0 1 1 1 0 0 0 0 0 0; 0 0 0 1 0 0 0 0 0 0 
 ;0 0 0 1 1 1 0 0 0 0; 0 0 0 0 0 1 0 0 0 0 ; 0 0 0 0 0 1 1 1 0 0 ; 0 0 0 0 0 0 0 1 0 0;
 0 0 0 0 0 0 0 1 1 1; 0 0 0 0 0 0 0 0 0 1]
Dtype=[1 1 1 1; 0 1 1 1; 0 0 1 0; 0 0 0 1]
C11=[1 1 1 1; 0 1 0 1 ; 0 0 1 1 ; 0 0 0 1]
star3=[1 1 1 1; 0 1 0 0; 0 0 1 0; 0 0 0 1 ]
star4=[1 1 1 1 1 ; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1]
star5=[1 1 1 1 1 1; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1]


#________________
VisualizeGquiverOfEnd(zig3)
VisualizeGquiverOfEnd(A5)
VisualizeGquiverOfEnd(C11)
VisualizeGquiverOfEnd(A4)
VisualizeGquiverOfEnd(Atype)
VisualizeGquiverOfEnd(Dtype)
VisualizeGquiverOfEnd(ComGrid)
VisualizeGquiverOfEnd(Atilda4)
VisualizeGquiverOfEnd(Atilda6)
VisualizeGquiverOfEnd(Atilda10)
VisualizeGquiverOfEnd(star4)
VisualizeGquiverOfEnd(star5)
VisualizeGquiverOfEnd(star3)

