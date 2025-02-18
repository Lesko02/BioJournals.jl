function horspool_find(jss::JournaledString, needle::LongDNA{4})
    # Yet to be coded
end

function horspool_find(jst::JSTree, needle::LongDNA{4})
    # Yet to be coded
end

function myers_ukkoken_find(jss::JournaledString, needle::LongDNA{4})
    # Yet to be coded
end

function myers_ukkoken_find(jst::JSTree, needle::LongDNA{4})
    # Yet to be coded
end


function approximate_findall(query, tolerance::Int64, seq::LongDNA{4})
    results = UnitRange{Int64}[]
    pos = findfirst(query, tolerance, seq)
    while pos !== nothing
        push!(results, pos)
        pos = findnext(query, tolerance, seq, last(pos)+tolerance+1)
    end
    return results
end

function approximate_search(jss::JournaledString, needle::LongDNA{4})
    query = ApproximateSearchQuery(needle)
    tolerance = ceil(Int64, (length(needle) / 100) * 5)  
    indexMatrix = Dict{Int64, Vector{UnitRange{Int64}}}()
    vector = approximate_findall(query, tolerance, jss.reference)
    to_remove = Set{UnitRange{Int64}}()
    to_add = Set{UnitRange{Int64}}()

    for i in 1:length(jss.deltaMap)
        indexMatrix[i] = vector
    end

    for i in 1:length(jss.deltaMap)
        empty!(to_add)
        empty!(to_remove)
        for range in indexMatrix[i]
            for ( _, entry) in jss.deltaMap[i]
                
                if entry.position in range
                    seq = apply_delta(jss.reference, entry)

                    for element in approximate_findall(query, tolerance, seq)
                        push!(to_add, element)
                    end
                    
                    push!(to_remove, range)
                end

            end
        end 
        indexMatrix[i]= filter(x -> all(y -> x != y, to_remove), indexMatrix[i])
        append!(indexMatrix[i], to_add)
        indexMatrix[i] = collect(Set(indexMatrix[i]))
    end
    return indexMatrix
end

function approximate_search(jss::JournaledString, needle::LongDNA{4},
    tol::Int64)

    if tol <= 0 || tol >= 100
        error("Tolerance cannot less or 0% or more than 100%")
    end

    tolerance = ceil(Int64, (length(needle) / 100) * tol) 
    query = ApproximateSearchQuery(needle)
    indexMatrix = Dict{Int64, Vector{UnitRange{Int64}}}()
    vector = approximate_findall(query, tolerance, jss.reference)
    to_remove = Set{UnitRange{Int64}}()
    to_add = Set{UnitRange{Int64}}()

    for i in 1:length(jss.deltaMap)
        indexMatrix[i] = vector
    end

    for i in 1:length(jss.deltaMap)
        empty!(to_add)
        empty!(to_remove)
        for range in indexMatrix[i]
            
            for ( _, entry) in jss.deltaMap[i]
                    
                if entry.position in range
                    seq = apply_delta(jss.reference, entry)

                    for element in approximate_findall(query, tolerance, seq)
                        push!(to_add, element)
                    end
                        
                    push!(to_remove, range)
                end

            end

        end 
        indexMatrix[i]= filter(x -> all(y -> x != y, to_remove), indexMatrix[i])
        append!(indexMatrix[i], to_add)
        indexMatrix[i] = collect(Set(indexMatrix[i]))
    end
    return indexMatrix
end

function approximate_search(jst::JSTree, needle::LongDNA{4})

    seq = jst.root
    query = ApproximateSearchQuery(needle)
    vector = UnitRange{Int}[]
    indexMatrix = Dict{String, Vector{UnitRange{Int64}}}(
        name => UnitRange{Int64}[] for name in keys(jst.children))
    tolerance = ceil(Int64, (length(needle) / 100) * 5) 

    vector = approximate_findall(query, tolerance, seq)
    to_remove = Set{UnitRange{Int64}}()
    to_add = Set{UnitRange{Int64}}()
    range2 = Set{UnitRange{Int64}}()

    for (name, _) in indexMatrix
        indexMatrix[name] = vector
    end

    println(vector)

    for (name, node) in jst.children

        for range in indexMatrix[name]

            empty!(to_add)
            empty!(to_remove)
            if !isnothing(node.deltaMap)
                empty!(range2)
                seq = flatten(jst, node.parent.name)
                range2 = push!(range2, range)
                range2 = union!(range2, approximate_findall(query, 
                tolerance, seq))
                union!(to_add, range2)
                for (time, entry) in node.deltaMap
                    for range3 in range2
                        if entry.position in range3
                            
                            seq = apply_delta(seq, entry)
                                
                            for element in approximate_findall(query,
                                tolerance, seq)
                                push!(to_add, element)
                            end
                            push!(to_remove, range3)   
                        end
                    end
                    
                end
            append!(indexMatrix[name], to_add)
            indexMatrix[name]= filter(x -> all(y -> x != y, to_remove), 
            indexMatrix[name])
                
        end
        indexMatrix[name] = collect(Set(indexMatrix[name]))
        end
    end
    return indexMatrix
end

function exact_search(jss::JournaledString, needle::LongDNA )
    results = Dict(i => UnitRange{Int64}[] for i in 1:length(jss.deltaMap))
    query = ExactSearchQuery(needle)
    vector = UnitRange{Int64}[]
    for i in 1:length(jss.deltaMap)

        seq = apply_delta(jss.reference, jss.deltaMap[i])
        vector = BioSequences.findall(query, seq)
        append!(results[i], vector)

    end
    return results
end

function exact_search(jst::JSTree, needle::LongDNA{4})
    indexMatrix = Dict{String, Vector{UnitRange{Int64}}}()
    query = ExactSearchQuery(needle)
    vector = UnitRange{Int}[]

    indexMatrix = Dict{String, Vector{UnitRange{Int64}}}(
        name => UnitRange{Int64}[] for name in keys(jst.children))
    
    for (name, child) in jst.children

        empty!(vector)
        if (!isnothing(child.deltaMap))
        seq = flatten(jst, name)
        vector = BioSequences.findall(query, seq)
        end

        indexMatrix[name] = append!(indexMatrix[name], vector)
        
    end
    return indexMatrix
end

