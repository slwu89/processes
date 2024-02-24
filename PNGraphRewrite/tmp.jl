using AlgebraicRewriting
using Catlab, AlgebraicPetri
using SpecialFunctions, Fleck
using Random
using Distributions
using Makie,CairoMakie

# --------------------------------------------------------------------------------
# given mean and sd, solve for the parameters of a Weibull distribution,
# so people believe us when we say it's actually non-Markov (we can cheat for Gamma with integer shape parameter)
# (from https://github.com/cran/mixdist/blob/master/R/weibullpar.R)
function weibullpar(mu, sigma)
    cv = sigma / mu
    if cv < 1e-06
        nu = cv / (sqrt(trigamma(1)) - cv * digamma(1))
        shape = 1 / nu
        scale = mu / (1 + nu * digamma(1))
    else
        aa = log(cv^2 + 1)
        nu = 2 * cv / (1 + cv)
        while true
            gb = (lgamma(1 + 2 * nu) - 2 * lgamma(1 + nu) - aa) / (2 * (digamma(1 + 2 * nu) - digamma(1 + nu)))
            nu -= gb
            if abs(gb) < 1e-12
                break
            end
        end
        shape = 1 / nu
        scale = exp(log(mu) - lgamma(1 + nu))
    end
    return shape, scale
end

# --------------------------------------------------------------------------------
# we want a PN with a marking
@present SchMarkedLabelledPetriNet <: SchLabelledPetriNet begin
    M::Ob
    m::Hom(M,S)
end

@acset_type MarkedLabelledPetriNet(SchMarkedLabelledPetriNet, index=[:it, :is, :ot, :os]) <: AbstractLabelledPetriNet

to_graphviz(SchMarkedLabelledPetriNet)

# a little function to see how many tokens are in each place
function marking(pn)
    names = Tuple(pn[:,:sname])
    vals = length.(incident(pn, parts(pn,:S), :m))
    return NamedTuple{names}(vals)
end

# --------------------------------------------------------------------------------
# make the set of rewrite rules

# function to make a DPO rule for a transition in a marked PN
function make_rule(pn::T, t) where {T<:MarkedLabelledPetriNet}
    input_m = inputs(pn, t)
    output_m = outputs(pn, t)

    # get L
    L = MarkedLabelledPetriNet{Symbol}()
    copy_parts!(L, pn, S=:)

    if length(input_m) > 0
        # add tokens
        add_parts!(
            L, :M, length(input_m),
            m = vcat(incident(L, pn[input_m, :sname], :sname)...)
        )
    end

    # get R
    R = MarkedLabelledPetriNet{Symbol}()
    copy_parts!(R, pn, S=:)

    if length(output_m) > 0
        # add tokens
        add_parts!(
            R, :M, length(output_m),
            m = vcat(incident(R, pn[output_m, :sname], :sname)...)
        )
    end

    _, span = only(maximum_common_subobject([L,R], abstract=false))

    return Rule(legs(first(span))[1], legs(first(span))[2])
end


# --------------------------------------------------------------------------------
# we want something to store the rules, clocks associated to each, and their type

@present SchClockSystem(FreeSchema) begin
    (Clock,Event)::Ob # clock is a single ID, event is an entire class of clocks (e.g. "typed" clocks)
    NameType::AttrType # each event has a name
    RuleType::AttrType # each event has a rewrite rule
    DistType::AttrType # each event has a function that returns a firing time when it's enabled
    MatchType::AttrType # each event has set of matches associated with it (in theory, this would instead be bijective with Clock)
    KeyType::AttrType
    event::Hom(Clock,Event)
    key::Attr(Clock,KeyType)
    name::Attr(Event,NameType)
    rule::Attr(Event,RuleType)
    dist::Attr(Event,DistType)
    match::Attr(Event,MatchType)
end

to_graphviz(SchClockSystem)

@acset_type ClockSystem(SchClockSystem, index=[:event], unique_index=[:name,:key])

# --------------------------------------------------------------------------------
# awful simulator function to avoid julia global weirdness

function run_sir(S,I,R,maxevent)

    # PN stores the model state
    # sirpn = @acset MarkedLabelledPetriNet{Symbol} begin
    #     S=3
    #     sname=[:S,:I,:R]
    #     T=6
    #     tname=[:inf,:rec,:birth,:deathS,:deathI,:deathR]
    #     I=6
    #     it=[1,1,2,4,5,6]
    #     is=[1,2,2,1,2,3]
    #     O=4
    #     ot=[1,1,2,3]
    #     os=[2,2,3,1]
    #     M=sum([S,I,R])
    #     m=[fill(1,S);fill(2,I);fill(3,R)]
    # end

    sirpn = @acset MarkedLabelledPetriNet{Symbol} begin
        S=3
        sname=[:S,:I,:R]
        T=7
        tname=[:inf,:rec,:birth,:deathS,:deathI,:deathR,:wane]
        I=7
        it=[1,1,2,4,5,6,7]
        is=[1,2,2,1,2,3,3]
        O=5
        ot=[1,1,2,3,7]
        os=[2,2,3,1,1]
        M=sum([S,I,R])
        m=[fill(1,S);fill(2,I);fill(3,R)]
    end

    # clock system stores the sampling method
    sirclock = @acset ClockSystem{Symbol,Rule,Function,Incremental.IncCCHomSet,Tuple{Int,Int,Int}} begin
        Event=nt(sirpn)
        name=sirpn[:,:tname]
    end

    # add rules
    for t in parts(sirclock,:Event)
        sirclock[t,:rule] = make_rule(sirpn, t)
    end

    # which rules are always enabled?
    always_enabled_t = findall(nparts.(codom.(getfield.(sirclock[:, :rule], :L)), :M) .== 0)

    # add clock distributions
    pop = nparts(sirpn, :M)
    lifespan = 65*365
    μ = 1/lifespan
    β = 0.001
    wane = 60

    # should somehow check that if this fn returns Exponential, its mapped to type Markov
    sirclock[only(incident(sirclock, :inf, :name)), :dist] = (t) -> Exponential(1 / β)
    sirclock[only(incident(sirclock, :birth, :name)), :dist] = (t) -> Exponential(1 / (μ*pop))
    sirclock[only(incident(sirclock, :deathS, :name)), :dist] = (t) -> Exponential(1 / μ)
    sirclock[only(incident(sirclock, :deathI, :name)), :dist] = (t) -> Exponential(1 / μ)
    sirclock[only(incident(sirclock, :deathR, :name)), :dist] = (t) -> Exponential(1 / μ)
    sirclock[only(incident(sirclock, :wane, :name)), :dist] = (t) -> Exponential(wane)

    α, θ = weibullpar(30, 5)
    # the non Markovian clock
    sirclock[only(incident(sirclock, :rec, :name)), :dist] = (t) -> Weibull(α,θ)

    # add incremental homsets
    for t in parts(sirclock,:Event)
        sirclock[t,:match] = IncHomSet(codom(sirclock[t,:rule].L), getfield.(sirclock[:,:rule], :R), sirpn, single=true)
    end

    # Fleck sampler
    sampler = FirstToFire{Tuple{Int,Int,Int}, Float64}()
    rng = Random.RandomDevice()
    tnow = 0.0

    for t in parts(sirclock,:Event)
        newkeys = collect(keys(sirclock[t,:match]))
        newkeys = map(newkeys) do k
            (t,k...)
        end
        add_parts!(
            sirclock, :Clock, length(newkeys),
            key = newkeys, event = fill(t, length(newkeys))
        )
        for c in newkeys
            enable!(sampler, c, sirclock[t, :dist](tnow), tnow, tnow, rng)
        end
    end

    # record initial state
    output = Vector{typeof((t=0.0,marking(sirpn)...))}()
    push!(output, (t=0.0,marking(sirpn)...))

    # when and what will happen next?
    (tnow, which) = next(sampler, tnow, rng)

    while length(output) < maxevent
        println("event $which fired at $tnow, total number of events: $(length(output))")
        event = first(which)
        update_maps = rewrite_match_maps(
            sirclock[event, :rule], 
            sirclock[event, :match][Pair(which[2:end]...)]
        )
        sirpn = codom(update_maps[:rh])
        push!(output, (t=tnow,marking(sirpn)...))

        # "always enabled" transitions need special treatment when they fire
        # because their hom-set won't change, the clocks won't be reset, so we do it manually
        if event ∈ always_enabled_t
            disable!(sampler, which, tnow)
            enable!(sampler, which, sirclock[event, :dist](tnow), tnow, tnow, rng)
        end        

        # update matches for all events
        for t in parts(sirclock, :Event)
            # update the hom-set
            del = Incremental.deletion!(sirclock[t,:match], update_maps[:kg])
            add = Incremental.addition!(sirclock[t,:match], event, update_maps[:rh], update_maps[:kh])
            del = map(del) do k
                (t,k...)
            end
            add = map(add) do k
                (t,k...)
            end
            # clocks that are disabled in the new marking (state):
            # 1: disable in the sampler
            # 2: remove from the clock system
            for c in del
                disable!(sampler, c, tnow)
            end        
            rem_parts!(sirclock, :Clock, sort(vcat(incident(sirclock, del, :key)...)))
            # clocks that are newly enabled in the new marking (state):
            # 1: enable in the sampler
            # 2: add to the clock system
            for c in add
                enable!(sampler, c, sirclock[t, :dist](tnow), tnow, tnow, rng)
            end
            add_parts!(sirclock, :Clock, length(add), key = add, event = t)             
            
        end
        (tnow, which) = next(sampler, tnow, rng)
    end

    return output
end

# non-markovian SIR with demography
sirout = run_sir(95,5,0,500)

f = Figure()
ax = Axis(f[1,1])
lines!(ax, first.(sirout), getindex.(sirout,2))
lines!(ax, first.(sirout), getindex.(sirout,3))
lines!(ax, first.(sirout), getindex.(sirout,4))
f

# can also run a simple birth-death model
sirout = run_sir(100,0,0,500)
lines(first.(sirout), getindex.(sirout,2))