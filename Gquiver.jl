#This is a code to visualize Gabriel quiver of END(G), where G is all the direct sum
#of the non-iso interval modules. 
# This code is the implementation of Proposition 4.11 of the paper
#"STABILIZATION OF THE SPREAD-GLOBAL DIMENSION" by BENJAMIN BLANCHETTE, JUSTIN DESROCHERS, ERIC J. HANSON, AND LUIS SCOCCOLA
#see also "EXACT STRUCTURES FOR PERSISTECE MODULES" by B. Blanchette, T. BrüStle, and E.J.Hanson.   

module GabrieQuiverOfEndG

#_______________________________
#poset mat to posetmat with zero diagonal   
function FullPosetToHasse(Poset)
    n = size(Poset, 1)
    Idd = Matrix(I, n, n)
    return Poset - Idd
end
#_______________________________
#有向グラフの隣接行列
function PosetToAdjMat(Poset)
    n = size(Poset, 1)
    rad = FullPosetToHasse(Poset)
    radsq = rad * rad
    Adj = Matrix(I, n, n)
    for i in 1:n
        for j in 1:n
            if radsq[i, j] == 0 && rad[i, j] != 0
                Adj[i, j] = 1
            else
                Adj[i, j] = 0
            end
        end
    end
    return Adj
end

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


function UpSet(PosetMat, Subset)
    n = size(PosetMat, 1)
    M = PosetMat^n
    M = M[Subset, :]
    Want = []
    for i in 1:n
        col = M[:, i]
        if any(!=(0), col)
            push!(Want, i)
        end
    end
    return Want
end

#_______________________________
#Take the up set of a subset of Poset
function DownSet(PosetMat,Subset)
    PosetMat=transpose(PosetMat)
return UpSet(PosetMat,Subset)
end
#_______________________________

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

function CheckInterval(Poset,Interval)
    if CheckGraphConvex(Poset,Interval) & CheckGraphConnectedness(Poset[Interval,Interval])
        return true
    else 
        return false
    end        
end    
#________________

function cocov(poset, interval)
    want = []
    n = size(poset)[1]
    for i in 1:n
        upper = setdiff(UpSet(poset, [i]), [i])
        lower = setdiff(DownSet(poset, interval),interval)
        if isempty(intersect(upper, lower)) && (i in lower )
            push!(want, i)
        end
    end
    return want
end


function cov(poset, interval)
    want = []
    n = size(poset)[1]
    for i in 1:n
        upper = setdiff(UpSet(poset,interval),interval)
        lower = setdiff(DownSet(poset,[i]),[i])
        if isempty(intersect(upper, lower)) && (i in upper )
            push!(want, i)
        end
    end
    return want
end


function IntInjectiveIrrTest(Poset, Interval)
    AllIntervals = IntervalsForMe(Poset)
    Want = []
    for cand in AllIntervals
        
        for s in cocov(Poset, cand)
            if  Set(Interval) ==  union(Set(cand), setdiff(Set(UpSet(Poset, [s])), Set(UpSet(Poset, cand))))
                push!(Want, sort(collect(cand)))
            end
        end
    end

    return Want
end

#________________

function IntSurjectiveIrr(Poset, Interval)
    Intervals = IntervalsForMe(Poset)
    m = length(Poset[1,:])
    n = length(Intervals)
    Want=[]
    for Inte in Intervals
        D = TakeMaxAntiChain(Poset,Inte)
        for x in 1 :m
            if DoesElemCoverConvexSet(Poset,Inte, x)
                if Set(Interval) == Set(union(Inte, LeftFookRightDown(Poset,D,[x]) ) ) 
                    push!(Want,Inte )
                end
            end
        end
    end
    return Want
end



function IntSurjectiveIrrTest(Poset, Interval) 
    Intervals = IntervalsForMe(Poset)
    m = size(Poset, 1)
    Want = []
    co = cov(Poset, Interval)
    for x in co 
        newset = union(Interval, setdiff(Set(DownSet(Poset, [x])), Set(DownSet(Poset, Interval))))
        push!(Want, sort(collect(newset)))  # Set を Vector に変換
    end
    return Want
end

function InFlow(Poset, Interval)
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
#_____
function GquiverOfEndRing(Poset)
    AllIntervals=IntervalsForMe(Poset)
    n= length(AllIntervals)
    ARQ=zeros(Int, n, n)
    i = 1
    for a in AllIntervals 
        proj=[0 for j in 1:n]
        if iszero(InFlow(Poset,a)) == false
            for int in InFlow(Poset,a)
                proj[nthoflist(AllIntervals,int)]=1
            end
            ARQ[:,i]=proj    
        end
        i=i+1
    end
    println(AllIntervals)
    return [ARQ,AllIntervals]   
end
#________________

function VisualizeGquiverOfEnd(Poset)
    AR = GquiverOfEndRing(Poset)
    ARQ=AR[1]
    AllInt=AR[2]
    ARQ = Matrix(ARQ)
    n= length(ARQ[1,:])
 return graphplot(ARQ, 
 names=[string(i) for i in AllInt ], 
 curvature_scalar=0.0,
 size=(2000,2000),
 nodeshape=:rect,
 markercolor = range(colorant"white", stop=colorant"white", length=n),
 linecolor = :darkgrey,
 fontsize = 15)
end
#________________

end #module end
