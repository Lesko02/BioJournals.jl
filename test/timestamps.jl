### -*- Mode: Julia -*-

### timestamps.jl

@testset "Timestamp" begin
    ts = Timestamp(3)

    @test ts == Timestamp(UInt64(3))
    @test ts == Timestamp(UInt8(3))
    @test ts == 3
    @test 3 == ts
    @test ts != -1
    @test -1 != ts
    @test ts < Timestamp(4)
    @test Timestamp(4) > ts
    @test ts < 4
    @test !(ts < -1)
    @test 2 < ts
    @test -1 < ts
    @test ts <= 3
    @test 3 <= ts
    @test -5 <= ts
    @test !(ts <= -1)
    @test ts <= Timestamp(3)
    @test ts <= Timestamp(4)
    @test !(Timestamp(4) <= ts)
    @test ts == big(3)
    @test big(3) == ts
    @test ts < big(4)
    @test big(2) < ts
    @test ts <= big(3)
    @test big(-1) < ts
    @test ts < big(typemax(UInt64)) + 1
    @test UInt64(ts) === UInt64(3)
    @test convert(Timestamp, 9) == Timestamp(9)
    @test convert(Int64, ts) == 3
    @test BioJournals.addto(ts, UInt64(2)) == Timestamp(5)
    @test BioJournals.addto(ts, Timestamp(2)) == Timestamp(5)
    @test ts + UInt32(2) == Timestamp(5)
    @test ts + UInt64(2) == Timestamp(5)
    @test hash(ts) == hash(Timestamp(3))
    @test isequal(ts, Timestamp(3))
    @test sprint(show, ts) == sprint(show, UInt64(ts))
    @test sprint(print, ts) == repr(UInt64(ts))
    @test_throws DomainError Timestamp(-1)

    @test isempty(DeltaMap())
    @test length(DeltaMap(2)) == 2
    @test length(JournaledString(dna"ACGT").deltaMap) == 10
end

### timestamps.jl ends here.
