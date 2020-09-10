#=
IBS_Master:
- Julia version: 1.5.0
- Author: vloge
- Date: 2020-09-09
=#


#Example genome inputs
#TAATGCCATGGGATGTT
#TAGTGCTAAGGGATGTC
#GATTTTCCCACGGCCAC
#GCCCGGATATGAGGTAA
#TTTTGGGCGGTTGCAAA
#TAAAATTAGGACATGGTGGCG

# note genome length more than 15 is difficult to compute

#-----------------------------
# PART -1 CODE: INITIALIZATON
#-----------------------------

mutable struct node
    kmer::String
    suffix::String
    prefix::String
    count::Int
end
node(x) = node(x,x[2:end],x[1:end-1],1)

printdata(self::node) = println("\nkmer = $(self.kmer); \nsuffix = $(self.suffix); \nprefix = $(self.prefix);  \ncount = $(self.count);  ")


println("Enter Genome = ")
Genome = readline()
println("Enter length of kmer = ")
klen = parse(Int,readline())

list = []
for i = 1:length(Genome) - klen + 1
    append!(list,[Genome[i:i+klen-1]])
end


function find_node(kmers,str,len)
    for i = 1:len
        if kmers[i].kmer == str
            return i
        end
    end
    return 0
end

using OrderedCollections
lexder = sort(list)

kmers = Array{node, 1}(undef, length(Set(lexder)))
no = 0
for i in lexder
    global no
    j = find_node(kmers,i,no)
    if j==0
        no = no + 1
        kmers[no] = node(i)
    else
        kmers[j].count = kmers[j].count + 1
    end
end

lexderset = sort(collect(Set(lexder)))

println("\n\n")
for i =1:length(kmers)
    printdata(kmers[i])
end


#-----------------------------
# PART -2 CODE: POSSIBILITY ANALYSIS
#-----------------------------

function Adjmat(kmers)
    l = length(kmers)
    A = zeros(Int64,(l,l))
    for i in 1:l
        curnode = kmers[i]
        #println("suffixes")
        for j in 1:l
            if curnode.suffix == kmers[j].prefix
                A[i,j] = 1
                #println("A[$i,$j]  $(curnode.kmer)  $(kmers[j].kmer)")
            end
        end
        #println("prefixes")
        for j in 1:l
            if curnode.prefix == kmers[j].suffix
                A[j,i] = 1
                #println("A[$i,$j]  $(kmers[j].kmer)  $(curnode.kmer)")
            end
        end
    end
    return A
end
A=Adjmat(kmers)


function degree(A)
    l = length(kmers)
    deg = Dict()
    for i in 1:l
        in_ = sum(A[:,i])
        out = sum(A[i,:])
        deg[kmers[i].kmer] = [in_,out]
    end
    return deg
end
function degree(A,X::String)              # additional methods
    i = find_node(kmers,X,length(kmer))
    in_ = sum(A[:,i])
    out = sum(A[i,:])
    deg = [in_,out]
    return deg
end
function degree(A,X::Int)                # additional methods
    i = X
    in_ = sum(A[:,i])
    out = sum(A[i,:])
    deg = [in_,out]
    return deg
end
degA = sort!(OrderedDict(degree(A)))

function checkstend(deg)
    # true means it is possible to form hamiltonian path
    arr = collect(values(deg))
    res = !(length(findall(x->x[1]==0,arr))>1 || length(findall(x->x[2]==0,arr))>1)
    return res
end
if !checkstend(degA)
    println("This case is not possible to solve as it more than one start or end nodes")
    println("Please dont execute the rest code, even if executed, they ")
end

using GraphRecipes, Plots
hamiltonian = graphplot(A, names=lexderset,arrow=arrow(:closed,:head,1,1), curvature_scalar=0.0, nodeshape=:circle, nodesize=0.15, nodecolor=:white )
png(hamiltonian,"hamiltonian")


#--------------------------------------------------------------------
# PART -3 CODE: PREDICTION OF START NODES, COMPLEX POSSIBLE RESULTS
#--------------------------------------------------------------------

function start_finder(deg)
    list = []
    arr = collect(values(deg))
    predefined = findall(x->x[1]==0,arr)
    if length(predefined) != 0
        for i in predefined
            append!(list,[kmers[i].kmer])
        end
    else
        for i = 1:length(arr)
            D = arr[i]
            if D[2] != 0    # if node has more out's than in's
                append!(list,[kmers[i].kmer])
            end
        end
    end
    return list
end
start = start_finder(degA)
#println(start)



# Analysis solving functions

function matches(str::String,list)
    mat = []
    list = collect(Set(list))
    for i in list
        if str[2:end] == i[1:end-1]
            append!(mat,[i])
        end
    end
    return mat
end

function del!(a, item)
    b = copy(a)
    q = findall(x->x==item, a)
    deleteat!(b,q[1])
    return b
end

function string_findall(str::String,sub::String)
    indexes = [0]
    for i = 1:length(str)
        if i == indexes[end] +1
            # println(str[i:end])
            ind = findfirst(sub, str[i:end])
            if ind == nothing
                break
            end
            append!(indexes,indexes[end]+ind[1])
        end
    end
    deleteat!(indexes,1)
    return indexes
end

function formatter(str::String)
    colen = sort(string_findall(str[1:end-1],";"),rev=true)
    if length(colen) != 0
        for i in colen
            if str[i+1] == ']' || str[i+1] == ';'
                str = string(str[1:i-1],str[i+1:end])
            end
        end
    end
    return str
end

# Main Recressive Function
function path_rec(start::String,stack)
    res = ""
    outs = matches(start,stack)
    #println(start," - ",outs, " - ",length(stack))

    if length(stack) == 0
        return res,true

    elseif length(outs) > 0 && length(stack)>0
        pos = false
        count = length(outs)
        for i in outs
            next = i
            sub = next[end]
            childstr, possibility = path_rec(next,del!(stack,next))
            #println(childstr, "   ", possibility)
            if possibility == false
                sub = ""
                count = count -1
                #println("this child is not possible")
            else
                pos = true
                sub = string(sub,childstr)
                #println("this child is possible")
            end
            res = string(res,sub)
            if count > 1 && pos
                res = string(res,";")
            end
        end
        if count > 1
            res = string("[",res,"]")
        end
        res = formatter(res)
        #println(res, "   ", pos)
        return res,pos
    else
        return res,false
    end
end



complex_sols = []
for i in start
    result,TF = path_rec(i,del!(lexder,i))
    if TF
        append!(complex_sols,[string(i,result)])
    end
    #println("Let start be $i -> result = $TF ")
end

println("\n")
println("Given genome  = $Genome \nComplex-Results ($(length(complex_sols))) are")
for i =1:length(complex_sols)
    println(complex_sols[i])
end

#--------------------------------------------------------------------
# PART -4 CODE: INTEPRETER FOR COMPLEX RESULTS TO SIMPLE RESULTS
#--------------------------------------------------------------------

# default ARGUMENTS:
# a is the array with simple solution string -> initially array with null string,
# nesting is the count for no of start brackets - "["  -> initially zero
# string given is stored in-memory IO buffer
function interpreter(io::IOBuffer, a = [""], nesting = 0)
    s = a[end]
    while !eof(io)
        c = read(io, Char)
        # println(a,"    ",s,"    ",c,"    ",nesting)
        if c == ']'
            if nesting == 0
                error("Unexpected \"]\"")
            else
                break
            end
        elseif c == ';' # if colen
            # it added new string element to array of string before '['
            push!(a, s)
        elseif c == '[' # if bracket starts
            #all previous string gets copied to var s in next iteration
            interpreter(io, a, nesting + 1)
        else
            # * means concatination
            a[end] *= c
        end
    end
    return a
end

interpreter(complex::String) = interpreter(IOBuffer(complex))

sols = []
for i =1:length(complex_sols)
    append!(sols,interpreter(complex_sols[i]))
end

println("\n")
println("Given genome  = $Genome \nFinal-Results ($(length(sols))) are")
for i =1:length(sols)
    println(sols[i])
end

# notes:
# both degree in,out should be equal
# remove start and end point as it is fixed
# repeation == Total degree / 2
# any node cant have more than 4 in and 4 out
# as our genome space is limited to [A G C T]

