export classic_sta_lta, classic_sta_lta!, delayed_sta_lta, delayed_sta_lta!
export recursive_sta_lta, recursive_sta_lta!, carl_sta_trig, carl_sta_trig!
export z_detect, z_detect!, trigger_onset

@doc """
  classic_sta_lta!(S, sta, lta)

Computes classic STA/LTA. 

The length of the STA and LTA are given by `sta` and `lta` in seconds, respectively.

!!! note

    The in-place version (classic_sta_lta!) overwrites the seismic trace with the STA/LTA 
    characteristic function. 

# Arguments
- `S::GphysData`: GphysData (SeisData, SeisChannel, SeisIO.Nodal.NodalChannel, SeisIO.Quake.EventChannel).
- `sta:Int`: Length of short time average window in seconds
- `lta:Int`: Length of long time average window in seconds
""" classic_sta_lta!
function classic_sta_lta!(S::GphysData, sta::Real, lta::Real)
    for ii in 1:S.n 
        classic_sta_lta!(S[ii], sta, lta)
    end
    return nothing
end
function classic_sta_lta!(C::GphysChannel, sta::Real, lta::Real)
    nsta = round(Int64, C.fs * sta)
    nlta = round(Int64, C.fs * lta)
    C.x .= classic_sta_lta(C.x, nsta, nlta)
    return nothing
end

@doc classic_sta_lta!
function classic_sta_lta(S::GphysData, sta::Real, lta::Real)
    U = deepcopy(S)
    classic_sta_lta!(U, sta, lta)
    return U 
end
function classic_sta_lta(C::GphysChannel, sta::Real, lta::Real)
    U = deepcopy(C)
    classic_sta_lta!(U, sta, lta)
    return U 
end
function classic_sta_lta(A::AbstractArray, nsta::Real, nlta::Real)
    nsta > nlta ? throw(DomainError((nsta, nlta), "Short-term average must be less than long-term average")) : nothing
    nsta < 0 ? throw(DomainError(nsta, "Short-term average must be less greater than zero")) : nothing
    nlta < 0 ? throw(DomainError(nlta, "Long-term average must be greater than zero")) : nothing
    N = size(A,1)
    T = eltype(A)
    nsta > N ? throw(DomainError(nsta, "Short-term average must be less than length of time series")) : nothing
    nlta > N ? throw(DomainError(nlta, "Long-term average must be less than length of time series")) : nothing

    # cumulative sum can be exploited to calculate a moving average 
    sta = cumsum(A .^ 2)
    lta = deepcopy(sta)

    # compute the STA and LTA 
    sta[nsta + 1:end] = sta[nsta + 1:end] .- sta[1:end-nsta]
    sta ./= nsta 
    lta[nlta + 1:end] = lta[nlta + 1:end] .- lta[1:end-nlta]
    lta ./= nlta 

    # pad zeros 
    sta[1:nlta] .= 0.0 

    # avoid division by zero by setting zero values to tiny float 
    tinyfloat = eps(T)
    ind = findall(lta .< tinyfloat)
    lta[ind] .= tinyfloat 
    return sta ./ lta 
end

@doc """
  delayed_sta_lta!(S, sta, lta)

Delayed STA/LTA from Withers, 1998. 

!!! note

    The in-place version (delayed_sta_lta!) overwrites the seismic trace with the STA/LTA 
    characteristic function.

# Arguments
- `S::GphysData`: GphysData (SeisData, SeisChannel, SeisIO.Nodal.NodalChannel, SeisIO.Quake.EventChannel).
- `sta::Real`: Length of short time average window in seconds
- `lta::Real`: Length of long time average window in seconds
""" delayed_sta_lta!
function delayed_sta_lta!(S::SeisData, sta::Real, lta::Real)
    for ii in 1:S.n
        delayed_sta_lta!(S[ii], sta, lta)
    end
    return nothing
end
function delayed_sta_lta!(C::GphysChannel, sta::Real, lta::Real)
    nsta = round(Int64, C.fs * sta)
    nlta = round(Int64, C.fs * lta)
    C.x .= delayed_sta_lta(C.x, nsta, nlta)
    return nothing
end

@doc delayed_sta_lta!
function delayed_sta_lta(S::GphysData, sta::Real, lta::Real)
    U = deepcopy(S)
    delayed_sta_lta!(U, sta, lta)
    return U 
end
function delayed_sta_lta(C::GphysChannel, sta::Real, lta::Real)
   U = deepcopy(C)
   delayed_sta_lta!(U, sta, lta)
   return U
end
function delayed_sta_lta(A::AbstractArray, nsta::Int, nlta::Int)
    nsta > nlta ? throw(DomainError((nsta, nlta), "Short-term average must be less than long-term average")) : nothing
    nsta < 0 ? throw(DomainError(nsta, "Short-term average must be less greater than zero")) : nothing
    nlta < 0 ? throw(DomainError(nlta, "Long-term average must be greater than zero")) : nothing
    N = size(A, 1)
    T = eltype(A)
    nsta > N ? throw(DomainError(nsta, "Short-term average must be less than length of time series")) : nothing
    nlta > N ? throw(DomainError(nlta, "Long-term average must be less than length of time series")) : nothing

    sta = zeros(T, N)
    lta = zeros(T, N)
    sta[1] = A[1] ^ 2 / nsta
    for ii in 2:N
        sta[ii] = (A[ii] ^ 2 + A[max(1, ii - nsta)] ^ 2) / nsta + sta[ii - 1]
        lta[ii] = (A[max(1, ii - nsta - 1)] ^ 2 + A[max(1, ii - nsta - nlta -1)] ^ 2) / nlta + lta[max(1, ii -1)]
    end
    sta[1:min(N,nlta + nsta + 50)] .= 0.0 
    lta[1:min(N,nlta + nsta + 50)] .= 1.0 # avoid division by zero 
    return sta ./ lta 
end

@doc """
  recursive_sta_lta!(S, sta, lta)

Computes recursive STA/LTA. 

The length of the STA and LTA are given by `sta` and `sta` in seconds, respectively. 

!!! note

    The in-place version (recursive_sta_lta!) overwrites the seismic trace with the STA/LTA 
    characteristic function. 

# Arguments
- `S::GphysData`: GphysData (SeisData, SeisChannel, SeisIO.Nodal.NodalChannel, SeisIO.Quake.EventChannel).
- `sta:Int`: Length of short time average window in seconds
- `lta:Int`: Length of long time average window in seconds

See also: [`classic_sta_lta`](@ref), [`carl_sta_trig`](@ref), [`delayed_sta_lta`](@ref) [`z_detect`](@ref)
""" recursive_sta_lta!
function recursive_sta_lta!(S::GphysData, sta::Real, lta::Real)
    for ii in 1:S.n
        recursive_sta_lta!(S[ii], sta, lta)
    end
    return nothing
end
function recursive_sta_lta!(C::GphysChannel, sta::Real, lta::Real)
    nsta = round(Int64, C.fs * sta)
    nlta = round(Int64, C.fs * lta)
    C.x .= recursive_sta_lta(C.x, nsta, nlta)
    return nothing
end

@doc recursive_sta_lta!
function recursive_sta_lta(S::GphysData, sta::Real, lta::Real) 
    U = deepcopy(S)
    recursive_sta_lta!(U, sta, lta)
    return U
end
function recursive_sta_lta(C::GphysChannel, sta::Real, lta::Real) 
    U = deepcopy(C)
    recursive_sta_lta!(U, sta, lta)
    return U
end
function recursive_sta_lta(A::AbstractArray, nsta::Int, nlta::Int)
    nsta > nlta ? throw(DomainError((nsta, nlta), "Short-term average must be less than long-term average")) : nothing
    nsta < 0 ? throw(DomainError(nsta, "Short-term average must be less greater than zero")) : nothing
    nlta < 0 ? throw(DomainError(nlta, "Long-term average must be greater than zero")) : nothing
    N = size(A,1)
    T = eltype(A)
    nsta > N ? throw(DomainError(nsta, "Short-term average must be less than length of time series")) : nothing
    nlta > N ? throw(DomainError(nlta, "Long-term average must be less than length of time series")) : nothing

    # compute short time average (STA) and long term average (LTA)
    # given by Evans and Allen 
    csta = 1.0 / nsta 
    clta = 1.0 / nlta 
    sta = 0.0 
    lta = 1e-99 # avoid zero division 
    charfct = zeros(eltype(A), size(A))
    icsta = 1.0 - csta 
    iclta = 1.0 - clta 

    for ii in 2:N
        sq = A[ii] ^ 2 
        sta = csta * sq + icsta * sta 
        lta = clta * sq + iclta * lta 
        charfct[ii] = sta / lta 

        if ii < nlta 
            charfct[ii] = 0.0 
        end
    end 
    return charfct 
end

@doc """
  carl_sta_trig!(S, sta, lta, ratio, quiet)

Computes the CarlSTAtrig characteristic function. 

eta = star - (ratio * ltar) - abs(sta - lta) - quiet

!!! note

    The in-place version (carl_sta_trig!) overwrites the seismic trace with the STA/LTA 
    characteristic function. 

# Arguments
- `S::GphysData`: GphysData (SeisData, SeisChannel, SeisIO.Nodal.NodalChannel, SeisIO.Quake.EventChannel).
- `sta::Real`: Length of short time average window in seconds
- `lta::Real`: Length of long time average window in seconds
- `ratio::AbstractFloat`: as `ratio` gets smaller, carl_sta_trig gets more sensitive
- `Quiet::AbstractFloat`: as `quiet` gets smaller, carl_sta_trig gets more sensitive
""" carl_sta_trig!
function carl_sta_trig!(S::GphysData, sta::Real, lta::Real, ratio::Real, quiet::Real)
    for ii in 1:S.n 
        carl_sta_trig!(S[ii],sta, lta, ratio, quiet)
    end
    return nothing
end
function carl_sta_trig!(C::GphysChannel, sta::Real, lta::Real, ratio::Real, quiet::Real)
    nsta = round(Int64, C.fs * sta)
    nlta = round(Int64, C.fs * lta)
    C.x .= carl_sta_trig(C.x, nsta, nlta, ratio, quiet)
    return nothing
end

@doc carl_sta_trig!
function carl_sta_trig(S::GphysData, sta::Real, lta::Real, ratio::Real, quiet::Real)
    U = deepcopy(S)
    carl_sta_trig!(U, sta, lta, ratio, quiet)
    return U
end
function carl_sta_trig(C::GphysChannel, sta::Real, lta::Real, ratio::Real, quiet::Real)
    U = deepcopy(C)
    carl_sta_trig!(U, sta, lta, ratio, quiet)
    return U
end
function carl_sta_trig(A::AbstractArray, nsta::Int, nlta::Int, ratio::Real, quiet::Real)
    nsta > nlta ? throw(DomainError((nsta,nlta), "Short-term average must be less than long-term average")) : nothing
    nsta < 0 ? throw(DomainError(nsta, "Short-term average must be less greater than zero")) : nothing
    nlta < 0 ? throw(DomainError(nlta, "Long-term average must be greater than zero")) : nothing
    N = size(A,1)
    T = eltype(A)
    nsta > N ? throw(DomainError(nsta, "Short-term average must be less than length of time series")) : nothing
    nlta > N ? throw(DomainError(nlta, "Long-term average must be less than length of time series")) : nothing
    sta = zeros(T,N)
    lta = deepcopy(sta)
    star = deepcopy(sta)
    ltar = deepcopy(sta)
    pad_sta = zeros(T, nsta)
    pad_lta = zeros(T, nlta)

    # compute the short time average (STA)
    for ii in 1:nsta
        sta .+= vcat(pad_sta, A[ii:N - nsta - 1 + ii])
    end
    sta ./= nsta

    # compute the long time average (LTA)
    for ii in 1:nlta
        lta .+= vcat(pad_lta, sta[ii:N - nlta - 1 + ii])
    end
    lta ./= nlta 
    lta = vcat(T(0), lta[1:end-1])
    
    # compute star, average of abs diff between trace and lta 
    for ii in 1:nsta
        star .+= vcat(pad_sta, abs.(A[ii:N - nsta - 1 + ii] .- lta[ii:N - nsta - 1 + ii]))
    end
    star ./= nsta 

    # compute ltar, long term average over star 
    for ii in 1:nlta 
        ltar .+= vcat(pad_lta, star[ii:N - nlta - 1 + ii])
    end 
    ltar ./= nlta 

    eta = star .- (T(ratio) .* ltar) .- abs.(sta .- lta) .- T(quiet)
    eta[1:nlta] .= T(-1.0) 
    return eta 
end

@doc """
  z_detect!(S, sta)

Z-detector from Withers, 1998. 

!!! note

    The in-place version (delayed_sta_lta!) overwrites the seismic trace with the STA/LTA 
    characteristic function.

# Arguments
- `S::GphysData`: GphysData (SeisData, SeisChannel, SeisIO.Nodal.NodalChannel, SeisIO.Quake.EventChannel).
- `sta::Real`: Length of short time average window in seconds
""" z_detect!
function z_detect!(S::GphysData, sta::Real)
    for ii in 1:S.n
        z_detect!(S[ii], sta)
    end
    return nothing 
end
function z_detect!(C::GphysChannel, sta::Real)
    nsta = round(Int64, C.fs * sta)
    C.x .= z_detect(C.x, nsta)
    return nothing
end

@doc z_detect!
function z_detect(S::GphysData, sta::Real)
    U = deepcopy(S)
    z_detect!(U, sta)
    return U
end
function z_detect(C::GphysChannel, sta::Real)
   U = deepcopy(C)
   z_detect!(U, sta)
   return U 
end
function z_detect(A::AbstractArray, nsta::Int)
    nsta < 0 ? throw(DomainError(nsta, "Short-term average must be less greater than zero")) : nothing
    N = size(A, 1)
    T = eltype(A)
    nsta > N ? throw(DomainError(nsta, "Short-term average must be less than length of time series")) : nothing

    # Z-detector given by Swindell and Snell, 1977 
    sta = zeros(T, N)
    A2 = A .^ 2

    # replace this with convolution
    for ii in 1:nsta 
        sta[ii:N - nsta - 1 + ii] .+= A2[ii:N - nsta - 1 + ii]
    end
    return (sta .- mean(sta)) ./ std(sta)
end

@doc """
  trigger_onset(charfct, thresh1, thresh2; max_len=1e99, max_len_delete=false)

Calculate trigger on and off times.

Given `thresh1` and `thresh2` calculate trigger on and off times from
characteristic function.
    
# Arguments
- `charfct::AbstractArray`: Characteristic function of e.g. STA/LTA trigger.
- `thresh1::Real`: Value above which trigger (of characteristic function)
    is activated (higher threshold).
- `thresh2::Real`: Value below which trigger (of characteristic function)
    is deactivated (lower threshold).
# Keywords
- `max_len::Int`: Maximum length of triggered event in samples. A new
    event will be triggered as soon as the signal reaches again above thresh1.
- `max_len_delete::Bool`: Do not write events longer than max_len into report file.
""" trigger_onset
function trigger_onset(charfct::AbstractArray, thresh1::Real, thresh2::Real;
    max_len::Int=10^10, max_len_delete::Bool=false,
)
    thresh2 > thresh1 ? throw(DomainError(thresh1, "thresh1 $thresh1 must be greater than thresh2 $thresh2")) : nothing
    ind1 = findall(charfct .> thresh1)
    if length(ind1) == 0
        return Array{Int64}(undef,0,0)
    end
    ind2 = findall(charfct .> thresh2)

    on = [ind1[1]]
    off = [-1]

    # determine the indices where charfct falls below off-threshold
    ind2_ = Array{Bool}(undef,length(ind2))
    ind2_[1:end-1] .= diff(ind2) .> 1
    # last occurence is missed by diff, add it manually
    ind2_[end] = true
    append!(off,ind2[ind2_])
    append!(on,ind1[findall(diff(ind1) .> 1) .+ 1])
    # include last pick if trigger is on or drop it
    if max_len_delete
        # drop it
        append!(off,max_len)
        append!(on,on[end])
    else
        # include it
        append!(off,ind2[end])
    end

    pick = []
    while on[end] > off[1]
        while on[1] <= off[1]
            deleteat!(on,1)
        end
        while off[1] < on[1]
            deleteat!(off,1)
        end
        if off[1] - on[1] > max_len
            if max_len_delete
                deleteat!(on,1)
                continue
            end
            prepend!(off,on[1] + max_len)
        end
        push!(pick,[on[1],off[1]])
    end
    return permutedims(hcat(pick...))
end