#This is a code to visualize Gabriel quiver of END(G), where G is all the direct sum
#of the non-iso interval modules. 
# This code is the implementation of Proposition 6.8. of the paper 
#"EXACT STRUCTURES FOR PERSISTECE MODULES" by B. Blanchette, T. BrüStle, and E.J.Hanson.       

module GabrieQuiverOfEndG


export VisualizeGquiverOfEnd


using GraphRecipes
using Plots,PyPlot
using PyCall
using Combinatorics
using LightGraphs, GraphPlot, LinearAlgebra
using AbstractAlgebra, Nemo
import Base.Iterators: flatten

#a<b iff a→b
ZZ = AbstractAlgebra.ZZ
function standardbasis(dimension,field)
    R = field
    n=dimension 
    Id = identity_matrix(R, n)
    listmat= [Id[:,i] for i in 1:n  ]
    return listmat
end


#_______________________________
function FullPosetMat(AdjMat)
    n = length(AdjMat[:,1])
    AdjOne = AdjMat + Matrix(I,n,n)
    AdjMatPower = Matrix(I,n,n)
    for i in 1:n
        AdjMatPower = AdjMatPower * AdjOne
    end
    for j in 1:n
        for k in 1:n
              if AdjMatPower[j,k] !=0
                AdjMatPower[j,k]=1
              end
        end
    end              
    return AdjMatPower
end

#_______________________________
#poset mat to posetmat with zero diagonal 
function FullPosetToHasse(Poset) ##poset to no-direction path
    n= length(Poset[1,:])
    Idd=MatrixSpace(ZZ,n,n)(identity_matrix(ZZ,n))
    Poset=MatrixSpace(ZZ,n,n)(Poset)
    return  (Poset -Idd)
end    
#_______________________________

#_______________________________
#有向グラフの隣接行列
function PosetToAdjMat(Poset)
    n= length(Poset[1,:])
    rad = FullPosetToHasse(Poset)
    rad=MatrixSpace(ZZ,n,n)(rad)
    radsq=rad^2
    Adj=MatrixSpace(ZZ,n,n)(identity_matrix(ZZ,n))
    for i in 1:n 
        for j in 1:n
            if radsq[i,j]==0 && rad[i,j]!=0
                Adj[i,j]=1
            else 
                Adj[i,j]=0
            end
        end          
    end    
return Adj 
end    
#_______________________________

#_______________________________
function MakeConvexhullofSub(Poset,Sub)
    n =length(Poset[1,:])
    vec=[0 for i in 1:n ]   
    for i in Sub
        vec[i]=1
    end   
    Down =  Poset *vec
    Down = [i for i in 1:n if Down[i]!=0 ]
    Up =    transpose(vec) * Poset
    Up = [i for i in 1:n if Up[i]!=0 ]

    conv=intersect(Up,Down)
    return conv
end
#________________
function CheckGraphConvex(Poset,Subset)
    conv = MakeConvexhullofSub(Poset,Subset)
    if Set(conv)== Set(Subset)
        return true
    else
        return false 
    end       
end

#________________
function PathIsThere(PosetMat,n,m)
    size= length(PosetMat[1,:])
    PosetMat=MatrixSpace(ZZ,size,size)(PosetMat)
    if  size < n || size < m 
        return "wrong size input"
    else
        if ((transpose(PosetMat)^size)*standardbasis(size,ZZ)[n])[m]==0
           return false 
        else 
            print(string("there is a path ",string(n) , " to " , string(m),"\n" ))
            return  true
        end
    end 
end    
#_______________________________

#_______________________________
##check n and m are comparable in a poset 
function Iscomparable(PosetMat,n,m)
    size= length(PosetMat[1,:])
    PosetMat=MatrixSpace(ZZ,size,size)(PosetMat)
    if  size < n || size <m 
        return "wrong size input"
    else
        if PathIsThere(PosetMat,n,m) == true || PathIsThere(PosetMat,m,n)
            return true
        else
            return false
        end
    end            
end    
#_______________________________

#_______________________________
# check Subset is antichain in a Poset.
function IsAntiChain(Subset,Poset)
    n=length(Subset)
    if n==1 
        return true

    else
        for i in 1:n-1
            for  j in (i+1):n 
                if Iscomparable(Poset,Subset[i],Subset[j])
                    return false  
                end
            end
        end 
        return true
    end   
end    
#_______________________________
#take minimal antichain from subset
function TakeMinAntiChain(PosetMat, Sub)
    PosetMat=FullPosetToHasse(PosetMat)
    new=[]
    for i in Sub
            if iszero([PosetMat[:,i][k] for k in Sub])
                push!(new,i)
            end 
    end    
return new
end

#_______________________________
function TakeMaxAntiChain(PosetMat, Sub)
    PosetMat=FullPosetToHasse(PosetMat)
    new=[]
    for i in Sub
            if iszero([PosetMat[i,:][k] for k in Sub])
                push!(new,i)
            end 
    end    
return new
end

#_______________________________
# Check elem <= a for all  a in A 
function IsPointSmallerThanSet(Poset,elem,A)
    Up=Poset[elem,:]
    cap=[Up[i] for i in A]
    if 0 in cap
        return false 
    else 
        return true
    end 
end
#_______________________________
# Check elem < a for all  a in A 
function IsPointStrictlySmallerThanSet(Poset,elem,A)
    n= length(Poset[1,:])
    Poset=MatrixSpace(ZZ,n,n)(Poset)
    Id=identity_matrix(ZZ,n) 
    Poset = Poset -Id

    Up=Poset[elem,:]
    cap=[Up[i] for i in A]
    if 0 in cap
        return false 
    else 
        return true
    end 
end

#_______________________________
# Check elem => a for all  a in A 
function IsPointLargerThanSet(Poset,elem,A)
    Down=Poset[:,elem]
    cap=[Down[i] for i in A]
    if 0 in cap
        return false 
    else 
        return true
    end 
end  

#________________
#check A<=B as antichain
function IsASmallerAnichainThanB(Poset,A,B)
    if  IsAntiChain(A,Poset) == false
        print( string(A, " is not antichain \n"))
        return false 
    end    
    if  IsAntiChain(B,Poset) == false
        print( string(B, " is not antichain \n"))
        return false 
    end    
    c=0

    for a in A 
        for b in B 
            if PathIsThere(Poset, a, b) 
                c+=1
                break
            end
        end
    end
    
    if c != length(A)
            println("there is a s.t. a!<= b for all b in B  ")
            return false
    end
    
    c=0

    for b in B 
        for a in A 
            if PathIsThere(Poset, a, b) 
                c+=1
                break
            end
        end
    end
    if c != length(B)
        println("there is b s.t. b!<= a for all a in A")
        return false
    end

            
    return true     
end
#_______________________________

#_______________________________
function IsALargerAnichainThanB(Poset,A,B)
    return IsASmallerAnichainThanB(transpose(Poset),A,B)
end

##________________

function CompareAntiChains(Poset,A,B)    
    if  IsAntiChain(A,Poset) == false
        println(A, " is not antichain \n")
        return false 
    end    
    if  IsAntiChain(B,Poset) == false
        println(B, " is not antichain \n" )
        return  false 
    end    

    if  IsALargerAnichainThanB(Poset,A,B) 
        print(string(A, " < or = ", B, "\n"))
        return true
    elseif   IsALargerAnichainThanB(Poset,B,A) 
        print(string(B, " < or = " ,A,"\n" ))
        return true
    else
        print("no relation \n")
        return false    
    end 
end    
#________________
# Take the down set of a subset of a Poset
# This is not used.
#Down(D) = <-infty, D>:={x <=D}:={x in P : x <=s for all s in subset  }
function DownSet(PosetMat,Subset)
    n = length(PosetMat[1,:])
    Subset = TakeMinAntiChain(PosetMat, Subset)
    Want = [i for i in 1:n if IsPointSmallerThanSet(PosetMat,i,Subset) ]
   
    return Want
end
##
#_______________________________
#Take the up set of a subset of Poset
function UpSet(PosetMat,Subset)
    PosetMat=transpose(PosetMat)
return DownSet(PosetMat,Subset)
end
#_______________________________

#Take a spread set <A,B> = {x in P: A<=x, x<=B }
# for A<=B
function SpreadsSet(Poset,A,B)#a<b iff a\to b
    n = length(Poset[1,:])
    if !IsASmallerAnichainThanB(Poset,A,B)
        return false
    end
    Want=[x for x in 1:n if IsPointLargerThanSet(Poset,x,A) & IsPointSmallerThanSet(Poset,x,B)]
return Want
end        
#_______________________________
#
function SpreadPairOfInterval(Poset,Interval)
    A=TakeMinAntiChain(Poset,Interval)
    B=TakeMaxAntiChain(Poset,Interval)
    return [A,B]
end    

#_______________________________
function RightFook(PosetMat,SubSet)
    All=[i for i in 1: length(PosetMat[1,:]) ]
    X=UpSet(PosetMat,SubSet)
    Hook=[i for i in All if (i in X) == false]    
    return Hook        
end
#________________

#________________
function RightFookTest(PosetMat,SubSet)
    All=[i for i in 1: length(PosetMat[1,:]) ]
    X=upperset(PosetMat,SubSet)
    Hook=[i for i in All if (i in X) == false]    
    return Hook        
end
#________________
function LeftFook(PosetMat,SubSet)
    n=length(PosetMat[1,:])
    All=[i for i in 1:n ]
    X=DownSet(PosetMat,SubSet)
    Hook=[i for i in All if (i in X) ==false  ]    
    return Hook        
end
#________________
function LeftFookTest(PosetMat,SubSet)
    n=length(PosetMat[1,:])
    All=[i for i in 1:n ]
    X=lowerset(PosetMat,SubSet)
    Hook=[i for i in All if (i in X ) == false  ]    
    return Hook        
end
#________________
function IsMinimal(PosetMat, elem)
    PosetMat=FullPosetToHasse(PosetMat)
    if iszero(PosetMat[:,elem])
        return true
    else
        return false    
    end   
end
#________________
function IsMaximal(PosetMat, elem)
    return IsMinimal(transpose(PosetMat),elem )
end
#________________
#take intervals
function IntervalsForMe(Poset)
    n=length(Poset[1,:]) 
    combis = collect(combinations(1:n))
    Intervals=[]
    for inter in combis
        if  CheckInterval(Poset,inter)
            push!(Intervals,inter)
        end    
    end
    return Intervals   
end    
#_______________________________

#________________
#<A,B<
#This is not needed
function LeftUpRightFook(Poset,AntiChainA,AntiChainB)
    if IsAntiChain(AntiChainA,Poset) & IsAntiChain(AntiChainB,Poset)
        XX=UpSet(Poset,AntiChainA)
        YY=RightFook(Poset,AntiChainB)
        ZZ=intersect(XX,YY)
        return ZZ
    end
    print("not antichain")
    return false
end
#________________
#<A,B<
#This is we want
function LeftUpRightFookTest(Poset,AntiChainA,AntiChainB)
    if IsAntiChain(AntiChainA,Poset) & IsAntiChain(AntiChainB,Poset)
        XX=upperset(Poset,AntiChainA)
        YY=RightFookTest(Poset,AntiChainB)
        ZZ=intersect(XX,YY)
        return ZZ
    end
    print("not antichain")
    return false
end
#________________
#________________
#>A,B>
#this code is not needed
function LeftFookRightDown(Poset,AntiChainA,AntiChainB)
    if IsAntiChain(AntiChainA,Poset) & IsAntiChain(AntiChainB,Poset)
        #if IsASmallerAnichainThanB(Poset,AntiChainA,AntiChainB )
            XX=DownSet(Poset,AntiChainB)
            YY=LeftFook(Poset,AntiChainA)
            ZZ=intersect(XX,YY)
            return ZZ
        #end

    end
    println("not A<=B")
    return false
end
#________________
#>A,B>
#This is we need
function LeftFookRightDownTest(Poset,AntiChainA,AntiChainB)
    if IsAntiChain(AntiChainA,Poset) & IsAntiChain(AntiChainB,Poset)
        #if IsASmallerAnichainThanB(Poset,AntiChainA,AntiChainB )
            XX=lowerset(Poset,AntiChainB)
            YY=LeftFookTest(Poset,AntiChainA)
            ZZ=intersect(XX,YY)
            return ZZ
        #end

    end
    println("not A<=B")
    return false
end
#________________
##copy AdjMat is important
function CheckGraphConnectedness(AdjMat)
    AdjMatTrans = transpose(AdjMat)
    n = length(AdjMat[:,1])
    visited = 0 * Vector(1:n)
    visit_queue = [1]
    while length(visit_queue) > 0 
        cur_i = visit_queue[1] 
        deleteat!(visit_queue, 1)
        visited[cur_i] = 1
        neighbors_i = AdjMat[:,cur_i] + AdjMatTrans[:,cur_i]
        for j in 1:n 
            if neighbors_i[j] > 0 && 0 == visited[j] 
                append!(visit_queue, j)
            end
        end
    end
    return sum(visited) == n
end
#________________
function DiGraphFromPoset(Poset)
    n = length(Poset[1,:])

    Adj = PosetToAdjMat(Poset)

    ID = identity_matrix(ZZ, n)
    return Adj + ID

end
#DiGraphFromPoset(Pos)
#________________

function SubDiGraphFromPoset(Poset,Subset)
    n = length(Poset[1,:])
    zerovec= [0 for i in 1:n ]
    Digraph = DiGraphFromPoset(Poset)
    for i in 1:n
        if i in Subset
           #do nothig
        else
            Digraph[i, :] = zerovec
            Digraph[:, i] = zerovec
        end
    end
    return Digraph

end
#________________
function ConnectedCompOfSubDiGraph(DiGraphh, SubDiGraph,elem)
    if SubDiGraph[elem, elem] == 0
        println("elem is not in SubGraph")
        return false 
    end 
    n = length(DiGraphh[:,1])
    SubDiGraph = MatrixSpace(ZZ,n,n)(SubDiGraph)
    SubGraph =  SubDiGraph + transpose(SubDiGraph)
    components=[]
    visited =[]
    visit_queue = [elem]

    while length(visit_queue) > 0 
        cur_i = visit_queue[1] 
        deleteat!(visit_queue, 1)
        push!(visited,cur_i)
        neighbors_i = SubGraph[cur_i,:]
        for j in 1:n 
            if ((neighbors_i[j] != 0) && ((j in visited) == false) )&&  ((j in visit_queue)  == false) 
                push!(visit_queue, j)
            end
        end
    end
    visited = Set(visited)
    visited=[i for i in visited]
    return visited 

end       
#________________

function ConnectedCompsOfSubDiGraph(DiGraph, SubDiGraph)
    k=length(DiGraph[1,:])
    subvertex=[i for i in 1:k if SubDiGraph[i,i] !=0 ]
    n= length(subvertex)
    m = 0
    elem=subvertex[1]
    i=1
    visited=[]
    components = Vector{Vector{Int64}}()
    while  m < n
        comp = ConnectedCompOfSubDiGraph(DiGraph,SubDiGraph,elem)
        push!(components,comp)
        visited= union!(visited,comp) 
        visited= Set(visited)
        m = length(visited)
        if m != n
            elem = [i for i in subvertex if !(i in visited) ][1]
        end
    end
     
    return components
end      
#________________

function IsConnectedCompnentOfSubDiGraph(DiGraph,SubDiGraph, subset )
    if  iszero(SubDiGraph)
        return false
    end    

    if Set(subset) in [Set(comp) for comp in ConnectedCompsOfSubDiGraph(DiGraph,SubDiGraph)]
        return true
    else 
        return false    
    end    
end
#________________
function lowerset(Poset,subset)
    n = length(Poset[1,:])
    www=[0 for i in 1:n]
    for i in subset
        www[i] = 1 
    end
    www =Poset*www 
    return [i for i in 1:n if www[i] != 0 ] 
end
#________________
function upperset(Poset, subset)
    return lowerset(transpose(Poset), subset )
end  
#________________
function DoesElemCoverConvexSet(Poset,conv, elem)
    if elem in conv
        return false
    end   
    Want = copy(conv)
    Want = push!(Want,elem) 
    Want = MakeConvexhullofSub(Poset,Want)
     if elem in lowerset(Poset,conv)
        return false
     end
    #if elem in UpSet(Poset,conv)
        if length(Want) == length(conv) + 1
           return true
        else 
           return false
        end
    #end 
    return false
end    
#________________

function CheckInterval(Poset,Interval)
    if CheckGraphConvex(Poset,Interval) & CheckGraphConnectedness(Poset[Interval,Interval])
        return true
    else 
        return false
    end        
end    
#________________
#________________
#This code is needed
function IntInjectiveIrrTest(Poset,Interval)
    n= length(Poset[1,:])
    AllIntervals= IntervalsForMe(Poset)
    A = TakeMinAntiChain(Poset,Interval)
    B = TakeMaxAntiChain(Poset,Interval)
    Want=[]
    #Want = Vector{Vector{Int64}}()
    CandidateIntervals = AllIntervals
    m = length(CandidateIntervals)
    zerovec=[0 for i in 1:n]
    for i in 1:m     
        C= TakeMinAntiChain(Poset, CandidateIntervals[i])
        DiGraphh=DiGraphFromPoset(Poset)
        SubDiGraph= SubDiGraphFromPoset(Poset,CandidateIntervals[i])

        for c in C
            SubDiGraphminusc=copy(SubDiGraph)
            SubDiGraphminusc[c,:] = zerovec
            SubDiGraphminusc[:,c] = zerovec

            if IsConnectedCompnentOfSubDiGraph(DiGraphh,SubDiGraphminusc,Interval)
                if  Set(CandidateIntervals[i])==Set(union(Interval,LeftUpRightFookTest(Poset,[c],A)))
                    push!(Want,CandidateIntervals[i] )
                end    
            end
        end
    
    end   

    return Want

end 

#________________
#________________
#This code is needed
function IntSurjectiveIrrTest(Poset, Interval)
    Intervals = IntervalsForMe(Poset)
    m = length(Poset[1,:])
    n = length(Intervals)
    Want=[]
    for Inte in Intervals
        D = TakeMaxAntiChain(Poset,Inte)
        for x in 1 :m
            if DoesElemCoverConvexSet(Poset,Inte, x)
                if Set(Interval) == Set(union(Inte, LeftFookRightDownTest(Poset,D,[x]) ) ) 
                    push!(Want,Inte )
                end
            end
        end
    end
    return Want

end

#_______________________
#This code is needed
function OutFlowTest(Poset, Interval)
    return union(IntSurjectiveIrrTest(Poset,Interval), IntInjectiveIrrTest(Poset,Interval))
end
#________________

function nthoflist(Alist,a)
    n =length(Alist)
    if a in Alist
        j=1
        for i in Alist
            if a == i
                return j
            end    

            j=j+1
        end
        return j-1
    else
        println("is not in the list")
        return false
    end
end    
#________________

#________________
#This code is needed
function GquiverOfEndRing(Poset)
    AllIntervals=IntervalsForMe(Poset)
    n= length(AllIntervals)
    ARQ=zero_matrix(ZZ,n,n)
    i=1
    for Int in AllIntervals 
        proj=[0 for j in 1:n]
        if iszero(OutFlowTest(Poset,Int)) == false
            for int in OutFlowTest(Poset,Int)
                proj[nthoflist(AllIntervals,int)]=1
            end
            ARQ[i,:]=proj    
        end
        i=i+1
    end
    println(AllIntervals)
    return [ARQ,AllIntervals]   
end
#________________

#________________
#This code is needed
function VisualizeGquiverOfEnd(Poset)
    AR = GquiverOfEndRing(Poset)
    ARQ=AR[1]
    AllInt=AR[2]
    ARQ = Matrix(ARQ)
    n= length(ARQ[1,:])
    
 return graphplot(ARQ, 
 names=[string(i) for i in AllInt ], 
 curvature_scalar=0.0,
 #size=(2000,2000), #if we want big size graph. 
 nodeshape=:rect,
 markercolor = range(colorant"white", stop=colorant"white", length=n),
 linecolor = :darkgrey,
 fontsize = 6)
end
#________________
#end


end #module end






##futurework
#(1)graphvizを用いて記述する．
#(2)グラフが一定になるようにする．
#(2)色を付ける(インプットが加群、resolutionの１項目)


