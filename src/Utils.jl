import DataStructures

"""
    strongly_connected_components(V::Array, E::Array{Tuple})

Strongly connected components of a directed graph in reverse topological
order. Translated from sympy/utilities/iteratbles.py/strongly_connected_components.


"""
function strongly_connected_components(V, E)
    Gmap = Dict(vi => [] for vi in V)
    for (v1, v2) in E
        push!(Gmap[v1], v2)
    end
    # Non-recursive Tarjan's algorithm:
    lowlink = Dict()
    indices = Dict()
    stack = DataStructures.OrderedDict() #########
    callstack = []
    components = []
    next_index = 1

    function start(v)
        index = length(stack) + 1
        indices[v] = lowlink[v] = index
        stack[v] = missing
        push!(callstack, (v, Gmap[v]))  #####
    end

    function finish(v1)
        # Finished a component?
        if lowlink[v1] == indices[v1]
            component = [pop!(stack)[1]]
            
            while component[end] != v1
                push!(component, pop!(stack)[1])
            end
            push!(components, reverse(component))
        end
        v2, _ = pop!(callstack)
        if length(callstack) > 0
            v1, _ = callstack[end]
            lowlink[v1] = min(lowlink[v1], lowlink[v2])
        end
    end

    for v in V
        if v in keys(indices)
            continue
        end
        start(v)
        while length(callstack) > 0
            v1, it1 = callstack[end]

            if length(it1) == 0
                v2 = missing
            else
                v2 = popfirst!(it1)  ######
            end
            # Finished children of v1?
            if v2 === missing
                finish(v1)
            # Recurse on v2
            elseif !(v2 in keys(indices))
                start(v2)
            elseif v2 in keys(stack)
                lowlink[v1] = min(lowlink[v1], indices[v2])
            end
        end
    end
    # Reverse topological sort order:
    return components
end

"""
Computes the period of an irreducible transition matrix
T must be a 2D array
"""
function breadth_first_search(T)
    n = size(T)[1]

    non_tree_edge_values = Set()
    visited = Set(1)
    
    overall_newly_visited = Set(1)
    level = Dict(1 => 0)
    current_level = 0
    done = false  # imitate a do-while loop
    while !done  # runs at most len(class_) times
        done = length(visited) == n
        current_level += 1
        new_overall_newly_visited = Set()
        
        # this loop and the while loop above run a combined len(class_) number of times.
        # so this triple nested loop runs through each of the n states once.
        for i in overall_newly_visited

            # the loop below runs len(class_) number of times
            # complexity is around about O(n * avg(len(class_)))
            newly_visited = Set(j for j in 1:n if T[i, j] != 0)
            new_overall_newly_visited = union(new_overall_newly_visited, newly_visited)

            new_tree_edges = setdiff(newly_visited, visited)
            for j in new_tree_edges
                level[j] = current_level
            end

            new_non_tree_edges = intersect(newly_visited, visited)
            new_non_tree_edge_values = Set(level[i]-level[j]+1 for j in new_non_tree_edges)

            non_tree_edge_values = union(non_tree_edge_values, new_non_tree_edge_values)
            visited = union(visited, new_tree_edges)
        end
        overall_newly_visited = new_overall_newly_visited
    end

    positive_ntev = Set(val_e for val_e in non_tree_edge_values if val_e > 0)
    if length(positive_ntev) == 0
        return n
    else
        return gcd(positive_ntev...)
    end
end
