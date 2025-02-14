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

function approximate_search(jst::JSTree, needle::LongDNA{4})
    # Yet to be coded
end

function slow_search(jss::JournaledString, needle::LongDNA )
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

function slow_search(jst::JSTree, needle::LongDNA{4})
    query = ExactSearchQuery(needle)
    vector = UnitRange{Int}[]
    for (name, child) in jst.children

        empty!(vector)
        if (!isnothing(child.deltaMap))
        seq = flatten(jst, name)
        vector = BioSequences.findall(query, seq)
        end
        
        if isempty(vector)
            println("No match at child: $name")
        else
            println("Match at child: $name")
            println("Ranges: ", vector)
        end
    end
end
