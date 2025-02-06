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
        pos = findnext(query, tolerance, seq, last(pos)+1)
    end
    return results
end

function approximate_search(jss::JournaledString, needle::LongDNA{4})
    query = ApproximateSearchQuery(needle)
    tolerance = ceil(Int64, (length(needle) / 100) * 5)  
    indices = UnitRange{Int64}[]
    vector = UnitRange{Int64}[]
    indices = approximate_findall(query, tolerance, jss.reference)
    to_remove = UnitRange{Int64}[]

    for range in indices
        for i in length(jss.deltaMap)
            for ( _, entry) in jss.deltaMap[i]
                empty!(vector)

                if entry.position in range
                    seq = apply_delta(jss.reference, entry)
                    vector = push!(approximate_findall(query, tolerance, seq))
                    
                    push!(to_remove, range)

                    if isempty(vector)
                        println("No match at Deltamap N° $i")
                    else
                        println("Match at Deltamap N° $i")
                        println("Ranges: ", vector)
                    end
                end
            end 
            filter!(x -> x ∉ to_remove, indices)
            if !isempty(indices)
                println("Match at Deltamap N° $i")
                println("Ranges: ", indices) 
            end
        
        end
    end
end

function approximate_search(jst::JSTree, needle::LongDNA{4})
    # Yet to be coded
end

function slow_search(jss::JournaledString, needle::LongDNA )

    query = ExactSearchQuery(needle)
    vector = UnitRange{Int64}[]
    for i in 1:length(jss.deltaMap)
        empty!(vector)
        seq = apply_delta(jss.reference, jss.deltaMap[i])
        vector = BioSequences.findall(query, seq)

        if isempty(vector)
            println("No match at Deltamap N° $i")
        else
            println("Match at Deltamap N° $i")
            println("Ranges: ", vector)
        end

    end

end

function slow_search(jst::JSTree, needle::LongDNA{4})
    query = ExactSearchQuery(needle)
    vector = UnitRange{Int64}[]

    for (name, child) in jst.children

        empty!(vector)

        if (!isnothing(child.deltaMap))
            seq = apply_delta(flatten(jst, name), child.deltaMap)
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
