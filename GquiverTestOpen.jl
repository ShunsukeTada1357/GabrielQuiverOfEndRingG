##This is a test code for the function VisualizeGquiverOfEnd
## input is "Poset Matrix", 
## i.e., if  there is a path from i to j on the graph, then the (i,j) component of the "Poset Matrix" is 1, otherwise 0.
## For example if the graph is 1→2→3, then
#poset=[1 1 1;0 1 1; 0 0 1 ]    
include("Gquiver.jl")
using .Gquiver


Pos = [1 0 0; 1 1 1; 0 0 1 ]
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
ComGrid25=[  1 1 1 1 1 1 1 1 1 1 
  0 1 1 1 1 0 1 1 1 1 
 0 0 1 1 1 0 0 1 1 1 
 0 0 0 1 1 0 0 0 1 1 
 0 0 0 0 1 0 0 0 0 1  
0 0 0 0 0 1 1 1 1 1  
0 0 0 0 0 0 1 1 1 1 
0 0 0 0 0 0 0 1 1 1  
0 0 0 0 0 0 0 0 1 1 
0 0 0 0 0 0 0 0 0 1 ]




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
star6=[1 1 1 1 1 1 1; 0 1 0 0 0 0 0; 0 0 1 0 0 0 0; 0 0 0 1 0 0 0; 0 0 0 0 1 0 0; 0 0 0 0 0 1 0; 0 0 0 0 0 0 1]
flower3 = [1 1 1 1 1 1 1 1 1 1; 0 1 0 1 0 0 0 0 0 0; 0 0 1 1 0 0 0 0 0 0; 0 0 0 1 0 0 0 0 0 0; 
           0 0 0 0 1 0 1 0 0 0; 0 0 0 0 0 1 1 0 0 0; 0 0 0 0 0 0 1 0 0 0;
           0 0 0 0 0 0 0 1 0 1; 0 0 0 0 0 0 0 0 1 1; 0 0 0 0 0 0 0 0 0 1]
#________________
Gquiver.VisualizeGquiverOfEnd(zig3)
Gquiver.VisualizeGquiverOfEnd(A5)
Gquiver.VisualizeGquiverOfEnd(C11)
Gquiver.VisualizeGquiverOfEnd(A4)
Gquiver.VisualizeGquiverOfEnd(Atype)
Gquiver.VisualizeGquiverOfEnd(Dtype)
Gquiver.VisualizeGquiverOfEnd(ComGrid)
Gquiver.VisualizeGquiverOfEnd(Atilda4)
Gquiver.VisualizeGquiverOfEnd(Atilda6)
Gquiver.VisualizeGquiverOfEnd(Atilda10)
Gquiver.VisualizeGquiverOfEnd(star4)
Gquiver.VisualizeGquiverOfEnd(star5)
Gquiver.VisualizeGquiverOfEnd(star3)

