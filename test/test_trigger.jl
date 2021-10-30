# tests for the triggering module 

@testset "Classic STA/LTA" begin 
    # use random SeisChannel with fs > 1.0 
    C = randSeisChannel(s=true,c=false,nx=10000,fs_min=1.)

    sta = rand() .* 10.0
    lta = sta * 2.0  

    # should type, length, output  
    @test isa(classic_sta_lta(C,sta,lta), Array)
    @test length(classic_sta_lta(C,sta,lta)) == 10000 
    @test eltype(classic_sta_lta(C,sta,lta)) == eltype(C.x)

    # test when sta > lta 
    @test_throws DomainError classic_sta_lta(C,lta,sta)

    # test SeisChannel shorter than STA or LTA 
    C = randSeisChannel(s=true,c=false,nx=100,fs_min=1.)
    sta = -0.5 # too short! 
    @test_throws DomainError classic_sta_lta(C,sta,lta)
    lta = -0.25 
    @test_throws DomainError classic_sta_lta(C,sta,lta)
end

@testset "Recursive STA/LTA" begin 
    # use random SeisChannel with fs > 1.0 
    C = randSeisChannel(s=true,c=false,nx=10000,fs_min=1.)

    sta = rand() .* 10.0
    lta = sta * 2.0  

    # should type, length, output  
    @test isa(recursive_sta_lta(C,sta,lta), Array)
    @test length(recursive_sta_lta(C,sta,lta)) == 10000 
    @test eltype(recursive_sta_lta(C,sta,lta)) == eltype(C.x)

    # test when sta > lta 
    @test_throws DomainError recursive_sta_lta(C,lta,sta)

    # test SeisChannel shorter than STA or LTA 
    C = randSeisChannel(s=true,c=false,nx=100,fs_min=1.)
    sta = -0.5 # too short! 
    lta = 10. 
    @test_throws DomainError recursive_sta_lta(C,sta,lta)
    lta = -0.25 
    @test_throws DomainError recursive_sta_lta(C,sta,lta)
end

@testset "Delayed STA/LTA" begin 
    # use random SeisChannel with fs > 1.0 
    C = randSeisChannel(s=true,c=false,nx=10000,fs_min=1.)

    sta = rand() .* 10.0
    lta = sta * 2.0  

    # should type, length, output  
    @test isa(delayed_sta_lta(C,sta,lta), Array)
    @test length(delayed_sta_lta(C,sta,lta)) == 10000 
    @test eltype(delayed_sta_lta(C,sta,lta)) == eltype(C.x)

    # test when sta > lta 
    @test_throws DomainError delayed_sta_lta(C,lta,sta)

    # test SeisChannel shorter than STA or LTA 
    C = randSeisChannel(s=true,c=false,nx=100,fs_min=1.)
    sta = -0.5 # too short! 
    lta = 10. 
    @test_throws DomainError delayed_sta_lta(C,sta,lta)
    lta = -0.25 
    @test_throws DomainError delayed_sta_lta(C,sta,lta)
end

@testset "Carl STA/LTA" begin 
    # use random SeisChannel with fs > 1.0 
    C = randSeisChannel(s=true,c=false,nx=10000,fs_min=1.)

    sta = rand() .* 10.0
    lta = sta * 2.0  
    ratio = 0.8 
    quiet = 0.8 

    # should type, length, output  
    @test isa(carl_sta_trig(C,sta,lta,ratio,quiet), Array)
    @test length(carl_sta_trig(C,sta,lta,ratio,quiet)) == 10000 
    @test eltype(carl_sta_trig(C,sta,lta,ratio,quiet)) == eltype(C.x)

    # test when sta > lta 
    @test_throws DomainError carl_sta_trig(C,lta,sta,ratio,quiet)

    # test SeisChannel shorter than STA or LTA 
    C = randSeisChannel(s=true,c=false,nx=100,fs_min=1.)
    sta = -0.5 # too short! 
    lta = 10. 
    @test_throws DomainError carl_sta_trig(C,sta,lta,ratio,quiet)
    lta = -0.25 
    @test_throws DomainError carl_sta_trig(C,sta,lta,ratio,quiet)
end

@testset "Z-detect STA/LTA" begin 
    # use random SeisChannel with fs > 1.0 
    C = randSeisChannel(s=true,c=false,nx=10000,fs_min=1.)

    sta = rand() .* 10.0
    lta = sta * 2.0  

    # should type, length, output  
    @test isa(z_detect(C,sta), Array)
    @test length(z_detect(C,sta)) == 10000 
    @test eltype(z_detect(C,sta)) == eltype(C.x)

    # test when sta > lta 
    @test_throws DomainError delayed_sta_lta(C,lta,sta)

    # test SeisChannel shorter than STA or LTA 
    C = randSeisChannel(s=true,c=false,nx=100,fs_min=1.)
    sta = -0.5 # too short! 
    lta = 10. 
    @test_throws DomainError z_detect(C,sta)
    lta = -0.25 
    @test_throws DomainError z_detect(C,sta)
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
    cft = classic_sta_lta(C,sta,lta)

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