# tests for the triggering module 

@testset "Classic STA/LTA" begin 
    # use random SeisChannel with fs > 1.0 
    S = randSeisData(s=1.0, c=0.0, nx=10000, fs_min=1.0)
    sta = rand() .* 10.0
    lta = sta * 2.0  
    U = classic_sta_lta(S, sta, lta)

    # should type, length, output  
    @test isa(U, SeisData)
    @test [eltype(U.x[ii]) for ii in 1:U.n] == [eltype(S.x[ii]) for ii in 1:S.n]
    @test [length(S[ii].x) for ii in 1:S.n] == [length(U[ii].x) for ii in 1:S.n]
    @test U.x != S.x 

    # test when sta > lta 
    @test_throws DomainError classic_sta_lta(S,lta,sta)

    # test in-place 
    classic_sta_lta!(S, sta, lta)
    @test U == S

    # test SeisChannel shorter than STA or LTA 
    S = randSeisData(s=1.0,c=0.0,nx=100,fs_min=1.)
    @test_throws DomainError classic_sta_lta(S,-0.5,lta)
    @test_throws DomainError classic_sta_lta(S,-0.5,-0.25)

end

@testset "Recursive STA/LTA" begin 
    # use random SeisChannel with fs > 1.0 
    S = randSeisData(s=1.0, c=0.0, nx=10000, fs_min=1.0)
    sta = rand() .* 10.0
    lta = sta * 2.0  
    U = recursive_sta_lta(S, sta, lta)

    # should type, length, output  
    @test isa(U, SeisData)
    @test [eltype(U.x[ii]) for ii in 1:U.n] == [eltype(S.x[ii]) for ii in 1:S.n]
    @test [length(S[ii].x) for ii in 1:S.n] == [length(U[ii].x) for ii in 1:S.n]
    @test U.x != S.x 

    # test when sta > lta 
    @test_throws DomainError recursive_sta_lta(S,lta,sta)

    # test in-place 
    recursive_sta_lta!(S, sta, lta)
    @test U == S

    # test SeisChannel shorter than STA or LTA 
    S = randSeisData(s=1.0,c=0.0,nx=100,fs_min=1.)
    @test_throws DomainError recursive_sta_lta(S,-0.5,lta)
    @test_throws DomainError recursive_sta_lta(S,-0.5,-0.25)
end

@testset "Delayed STA/LTA" begin 
        # use random SeisChannel with fs > 1.0 
        S = randSeisData(s=1.0, c=0.0, nx=10000, fs_min=1.0)
        sta = rand() .* 10.0
        lta = sta * 2.0  
        U = delayed_sta_lta(S, sta, lta)
    
        # should type, length, output  
        @test isa(U, SeisData)
        @test [eltype(U.x[ii]) for ii in 1:U.n] == [eltype(S.x[ii]) for ii in 1:S.n]
        @test [length(S[ii].x) for ii in 1:S.n] == [length(U[ii].x) for ii in 1:S.n]
        @test U.x != S.x 
    
        # test when sta > lta 
        @test_throws DomainError delayed_sta_lta(S,lta,sta)
    
        # test in-place 
        delayed_sta_lta!(S, sta, lta)
        @test U == S
    
        # test SeisChannel shorter than STA or LTA 
        S = randSeisData(s=1.0,c=0.0,nx=100,fs_min=1.)
        @test_throws DomainError delayed_sta_lta(S,-0.5,lta)
        @test_throws DomainError delayed_sta_lta(S,-0.5,-0.25)
end

@testset "Carl STA/LTA" begin 
        # use random SeisChannel with fs > 1.0 
        S = randSeisData(s=1.0, c=0.0, nx=10000, fs_min=1.0)
        sta = rand() .* 10.0
        lta = sta * 2.0  
        ratio = 0.8 
        quiet = 0.8
        U = carl_sta_trig(S, sta, lta, ratio, quiet)
    
        # should type, length, output  
        @test isa(U, SeisData)
        @test [eltype(U.x[ii]) for ii in 1:U.n] == [eltype(S.x[ii]) for ii in 1:S.n]
        @test [length(S[ii].x) for ii in 1:S.n] == [length(U[ii].x) for ii in 1:S.n]
        @test U.x != S.x 
    
        # test when sta > lta 
        @test_throws DomainError carl_sta_trig(S, lta, sta, ratio, quiet)
    
        # test in-place 
        carl_sta_trig!(S, sta, lta, ratio, quiet)
        @test U == S
    
        # test SeisChannel shorter than STA or LTA 
        S = randSeisData(s=1.0,c=0.0,nx=100,fs_min=1.)
        @test_throws DomainError carl_sta_trig(S,-0.5,lta,ratio,quiet)
        @test_throws DomainError carl_sta_trig(S,-0.5,-0.25,ratio,quiet)
end

@testset "Z-detect STA/LTA" begin 
        # use random SeisChannel with fs > 1.0 
        S = randSeisData(s=1.0, c=0.0, nx=10000, fs_min=1.0)
        sta = rand() .* 10.0
        U = z_detect(S, sta)
    
        # should type, length, output  
        @test isa(U, SeisData)
        @test [eltype(U.x[ii]) for ii in 1:U.n] == [eltype(S.x[ii]) for ii in 1:S.n]
        @test [length(S[ii].x) for ii in 1:S.n] == [length(U[ii].x) for ii in 1:S.n]
        @test U.x != S.x 
    
        # test in-place 
        z_detect!(S, sta)
        @test U == S
    
        # test SeisChannel shorter than STA or LTA 
        S = randSeisData(s=1.0,c=0.0,nx=100,fs_min=1.)
        @test_throws DomainError z_detect(S,-0.5)
end

@testset "Trigger on/off" begin 
    # use random SeisChannel with fs > 1.0 
    C = randSeisChannel(s=true,c=false,nx=10000,fs_min=1.)
    # make a fake earthquake 
    C.x[4000:6000] .*= exp.(-log.(0.001:0.0005:1.001)) .+ 1
    C.x[7000:9000] .*= exp.(-log.(0.001:0.0005:1.001)) .+ 1

    # choose random window 
    sta = rand() .* 10.0
    lta = sta * 5.0  
    classic_sta_lta!(C,sta,lta)
    cft = C.x

    # on/off thresholds 
    thresh1 = 4.0
    thresh2 = 2.0

    @test isa(trigger_onset(cft,thresh1,thresh2),Matrix{Int64})

    # test with incredibly high threshold 
    thresh1 = Inf
    @test isa(trigger_onset(cft,thresh1,thresh2),Matrix)
    @test isempty(trigger_onset(cft,thresh1,thresh2)) 

    # test max event delete 
    thresh1 = 3.5 
    @test isempty(trigger_onset(cft,thresh1,thresh2,max_len=2,max_len_delete=true)) 
end