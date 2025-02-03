function horspool_find(jss::JournaledString, needle::LongDNA{4})
    # Yet to be coded
end

function horspool_findfind(jst::JSTree, needle::LongDNA{4})
    # Yet to be coded
end

function myers_ukkoken_find(jss::JournaledString, needle::LongDNA{4})
    # Yet to be coded
end

function myers_ukkoken_find(jst::JSTree, needle::LongDNA{4})
    # Yet to be coded
end

function naive_search(jss::JournaledString, needle::LongDNA )

    query = ExactSearchQuery(needle)
    vector = UnitRange{Int}[]
    for i in 1:length(jss.deltaMap)
        empty!(vector)
        seq = apply_delta(jss.reference, jss.deltaMap[i])
        vector = BioSequences.findall(query, seq)

        if isempty(vector)
            println("No match at Deltamap N° $i")
        else
            println("Match at Deltamap N° $i")
            print(" Ranges: ", vector)
        end

    end

end

function naive_find(jst::JSTree, needle::LongDNA{4})
    # Yet to be coded
end