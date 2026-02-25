#This is a code to visualize Gabriel quiver of the ring automorphism ring END(G), where G is all the direct sum of the non-iso interval modules. 
# This code is the implementation of Proposition 4.11 of the paper
#"STABILIZATION OF THE SPREAD-GLOBAL DIMENSION" by BENJAMIN BLANCHETTE, JUSTIN DESROCHERS, ERIC J. HANSON, AND LUIS SCOCCOLA
#see also "EXACT STRUCTURES FOR PERSISTECE MODULES" by B. Blanchette, T. BrüStle, and E.J.Hanson.   

module Gquiver 
using Plots,PyPlot
using GraphRecipes #simple graph 
using Combinatorics
using LinearAlgebra

#Input poset is a matrix M s.t. M(i,j) = 1 iff i<= j,  and  M(i,j) = 0 iff i !<= j.
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

function CheckGraphConnectedness(AdjMat)
    AdjMatTrans = transpose(AdjMat)              
    n = size(AdjMat, 1)                           
    visited = zeros(Int, n)                       
    visit_queue = [1]                             

    while !isempty(visit_queue)
        cur_i = popfirst!(visit_queue)    
    
        
        visited[cur_i] = 1                        
        neighbors_i = AdjMat[:,cur_i] + AdjMatTrans[:,cur_i]  

        for j in 1:n
            if neighbors_i[j] > 0 && visited[j] == 0
                push!(visit_queue, j)             
            end
        end
    end

    return all(visited .== 1)                     
end



function UpSet(Poset, Subset)
    n = size(Poset, 1)
    M = Poset[Subset, :] 
    Want = Int[]
    for i in 1:n
        if any(!=(0),  M[:, i])# Subset のどれかから i に到達できる
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
    combis = filter(!isempty, collect(combinations(1:n)))
    #combis = collect(combinations(1:n))
    Intervals=[]
    for inter in combis
        if  CheckInterval(Poset,inter)
            push!(Intervals,inter)
        end    
    end
    #@show Intervals
    #println(Intervals)
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

function cocov(poset, interval)#cocov(S) := max(↓S \ S).
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


function cov(Poset, interval)#cov(S):=min(↑S∖S)
    want = []
    upper = setdiff(UpSet(Poset,interval),interval)
    n = size(Poset,1)
    for i in upper
        lower_i = DownSet(Poset,[i])
        if intersect(upper, lower_i) == [i]
            push!(want, i)
        end
    end
    return want
end

function cocov(Poset, interval)#cocov(S) := max(↓S \ S).
    return cov(transpose(Poset), interval)
end


function IntInjectiveIrr(Poset, Interval,AllIntervals)# collect intervals J s.t. J ->Interval is irr injective.
    #AllIntervals = IntervalsForMe(Poset)
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




####改善できる？
function IntSurjectiveIrr(Poset, Interval) 
    Want = []
    co = cov(Poset, Interval)
    for x in co 
        newset = union(Interval, setdiff(Set(DownSet(Poset, [x])), Set(DownSet(Poset, Interval))))
        push!(Want, sort(collect(newset)))  # Set を Vector に変換
    end
    return Want
end

function InFlow(Poset, Interval,AllIntervals)
    return union(IntSurjectiveIrr(Poset,Interval), IntInjectiveIrr(Poset,Interval,AllIntervals))
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
    AllIntervals = IntervalsForMe(Poset)
    n = length(AllIntervals)
    ARQ = zeros(Int, n, n)

    for (i, a) in enumerate(AllIntervals)         # i = 終点のindex
        for int in InFlow(Poset, a, AllIntervals)               # int -> a
            u = nthoflist(AllIntervals, int)      # u = 始点のindex
            ARQ[u, i] = 1                         # u -> i
        end
    end

    return [ARQ, AllIntervals]
end
#________________


# ARQ: n×n adjacency matrix (0/1). edge i -> j if ARQ[i,j] != 0
# We have edge information
function find_directed_cycle(ARQ)
    n = size(ARQ, 1)
    @assert size(ARQ, 2) == n

    color  = zeros(Int, n)   # 0=unvisited, 1=visiting(in stack), 2=done
    parent = fill(0, n)

    # 1つのサイクルを復元する補助
    function build_cycle(u, v)
        # back-edge u -> v (v is ancestor in current DFS stack)
        cyc = Int[v]
        cur = u
        while cur != v && cur != 0
            push!(cyc, cur)
            cur = parent[cur]
        end
        push!(cyc, v)        # 閉じる
        reverse!(cyc)         # v ... v の順に
        return cyc
    end

    function dfs(u)
        color[u] = 1
        for v in 1:n
            ARQ[u, v] == 0 && continue
            if color[v] == 0
                parent[v] = u
                cyc = dfs(v)
                cyc === nothing || return cyc
            elseif color[v] == 1
                # back-edge found -> cycle exists
                return build_cycle(u, v)
            end
        end
        color[u] = 2
        return nothing
    end

    for s in 1:n
        if color[s] == 0
            cyc = dfs(s)
            cyc === nothing || return cyc
        end
    end

    return nothing
end

# ARQ[u,v] != 0 means edge u -> v
# returns Vector{Vector{Int}} of simple directed cycles
# We lose edge information
function all_directed_cycles(ARQ)
    n = size(ARQ, 1)
    @assert size(ARQ, 2) == n

    # adjacency list
    adj = [Int[] for _ in 1:n]
    for u in 1:n
        for v in 1:n
            ARQ[u, v] != 0 && push!(adj[u], v)
        end
    end

    cycles  = Vector{Vector{Int}}()
    blocked = falses(n)
    B       = [Set{Int}() for _ in 1:n]
    stack   = Int[]

    function unblock(u)
        blocked[u] = false
        for w in collect(B[u])
            delete!(B[u], w)
            blocked[w] && unblock(w)
        end
    end

    # search cycles starting/ending at s, within "allowed" vertex set
    function circuit(v, s, allowed::BitVector)
        found = false
        push!(stack, v)
        blocked[v] = true

        for w in adj[v]
            allowed[w] || continue
            if w == s
                push!(cycles, copy(stack))  # one cycle found
                found = true
            elseif !blocked[w]
                found |= circuit(w, s, allowed)
            end
        end

        if found
            unblock(v)
        else
            for w in adj[v]
                allowed[w] || continue
                push!(B[w], v)
            end
        end

        pop!(stack)
        return found
    end

    # To avoid duplicates, restrict to vertices >= s (by index)
    for s in 1:n
        allowed = falses(n)
        allowed[s:n] .= true

        blocked .= false
        for i in 1:n
            empty!(B[i])
        end
        empty!(stack)

        circuit(s, s, allowed)
    end

    # Canonicalize each cycle: rotate so the smallest vertex comes first
    function canonical(cyc::Vector{Int})
        k = length(cyc)
        mpos = argmin(cyc)
        return vcat(cyc[mpos:end], cyc[1:mpos-1])
    end

    canon = [canonical(c) for c in cycles]
    # unique (string key is simple and robust)
    keys = Set{String}()
    out  = Vector{Vector{Int}}()
    for c in canon
        key = join(c, ",")
        if !(key in keys)
            push!(keys, key)
            push!(out, c)
        end
    end

    return out
end


function Isdirected_with_cycle(Poset)
    ARQ = GquiverOfEndRing(Poset)[1]
    cyc = find_directed_cycle(ARQ)
    if cyc === nothing
        return true, nothing
    else
        println("It is cyclic. One cycle is: ", cyc)
        return false, cyc
    end
end

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
 linecolor = :darkgrey,
 fontsize = 15,
 markercolor = fill(colorant"white", n)
  )
end



# cyc: [1,2,3,1] のように閉じた列でも、[1,2,3] のように閉じてなくてもOK
function cycle_to_vertices_edges(cyc::Vector{Int})
    isempty(cyc) && return (Int[], Tuple{Int,Int}[])

    # 末尾が先頭と同じなら落として扱う（閉じている表現を正規化）
    c = (length(cyc) >= 2 && cyc[end] == cyc[1]) ? cyc[1:end-1] : cyc

    # 頂点（順序を保って重複除去）
    seen = falses(maximum(c))
    verts = Int[]
    for v in c
        if v > length(seen)
            resize!(seen, v); seen[v] = false
        end
        if !seen[v]
            push!(verts, v)
            seen[v] = true
        end
    end

    # 辺（連続ペア + 最後→最初 で閉じる）
    edges = Tuple{Int,Int}[]
    for i in 1:length(c)-1
        push!(edges, (c[i], c[i+1]))
    end
    push!(edges, (c[end], c[1]))

    return verts, edges
end

function VisualizeGquiverOfEnd_cycle(Poset, cycle::Vector{Int};
                                    vcolor_cycle=colorant"gold",
                                    ecolor_cycle=colorant"red")
    ARQ, AllInt = GquiverOfEndRing(Poset)
    ARQ = Matrix(ARQ)
    n = size(ARQ, 1)

    verts, hedges = cycle_to_vertices_edges(cycle)

    # 頂点色
    vcolor = fill(colorant"white", n)
    for v in verts
        1 <= v <= n || continue
        vcolor[v] = vcolor_cycle
    end

    # 辺色：ARQ から辺リストを作り、その順で色ベクトルを作る
    edges_u = Int[]
    edges_v = Int[]
    for u in 1:n, v in 1:n
        if ARQ[u,v] != 0
            push!(edges_u, u)
            push!(edges_v, v)
        end
    end

    ecolor = fill(colorant"darkgrey", length(edges_u))
    hset = Set(hedges)
    for k in eachindex(edges_u)
        if (edges_u[k], edges_v[k]) in hset
            ecolor[k] = ecolor_cycle
        end
    end

    return graphplot(
        ARQ,
        names=[string(i) for i in AllInt],
        curvature_scalar=0.0,
        size=(2000,2000),
        nodeshape=:rect,
        fontsize=15,
        markercolor=vcolor,
        linecolor=ecolor,
    )
end

#########For Graph viz
function arq_to_dot(ARQ, labels::Vector{String};
                    highlight_vertices=Int[],
                    highlight_edges=Tuple{Int,Int}[])
    n = size(ARQ, 1)
    @assert size(ARQ, 2) == n
    @assert length(labels) == n

    hv = Set(highlight_vertices)
    he = Set(highlight_edges)

    io = IOBuffer()
    println(io, "digraph G {")
    println(io, "  rankdir=LR;")
    println(io, "  node [shape=box, style=filled, fontname=\"Helvetica\"];")
    println(io, "  edge [color=\"#666666\"];")

    # nodes
    for v in 1:n
        fill = (v in hv) ? "#FFD966" : "white"   # gold / white
        # label はダブルクォートを避ける
        lab = replace(labels[v], "\"" => "\\\"")
        println(io, "  $v [label=\"$lab\", fillcolor=\"$fill\"];")
    end

    # edges
    for u in 1:n, v in 1:n
        ARQ[u,v] == 0 && continue
        if (u,v) in he
            println(io, "  $u -> $v [color=\"red\", penwidth=3.0];")
        else
            println(io, "  $u -> $v;")
        end
    end

    println(io, "}")
    return String(take!(io))
end

function render_dot(dot::String; outpath="gquiver.png")
    dotpath = replace(outpath, r"\.[^\.]+$" => ".dot")
    open(dotpath, "w") do f
        write(f, dot)
    end
    # PNGなら -Tpng, PDFなら -Tpdf
    run(`dot -Tpng $dotpath -o $outpath`)
    return outpath
end


###

function output(Poset, name::String)
    ARQ, AllInt = GquiverOfEndRing(Poset)
    labels = [string(x) for x in AllInt]
    ok, cyc = Gquiver.Isdirected_with_cycle(Poset)
    if !ok && cyc !== nothing
        verts, edges = cycle_to_vertices_edges(cyc)
        dot = arq_to_dot(ARQ, labels; highlight_vertices=verts, highlight_edges=edges)
        render_dot(dot; outpath="gquiver_cycle_"*name*".png")
    else
        dot = arq_to_dot(ARQ, labels)
        render_dot(dot; outpath="gquiver_"*name*".png")
    end
end
end # module Gquiver
