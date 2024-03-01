using AlgebraicRewriting
using Catlab, AlgebraicPetri
using SpecialFunctions, Fleck
using Random
using Distributions
using Makie,CairoMakie

"""
    Get the shape and scale parameters for a Weibull distribution with a given mean and variance.
This code is a translation of the R "mixdist" package's function at https://github.com/cran/mixdist/blob/master/R/weibullpar.R
"""
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
# The PN schemas and types we are interested in

"""
    This is a presentation of a Petri net with a marking (a set of tokens, and for each token, an assignment to a place it is at).
"""
@present SchMarkedLabelledPetriNet <: SchLabelledPetriNet begin
    M::Ob
    m::Hom(M,S)
end

"""
    An abstract type for all marked labelled Petri nets.
"""
@abstract_acset_type AbstractMarkedLabelledPetriNet <: AbstractLabelledPetriNet

"""
    This is type for an acset to store instances of `SchMarkedLabelledPetriNet`, which is a subtype of `AbstractMarkedLabelledPetriNet`.
"""
@acset_type MarkedLabelledPetriNet(SchMarkedLabelledPetriNet, index=[:it, :is, :ot, :os]) <: AbstractMarkedLabelledPetriNet

# visualize the schema of the PNs we are considering here
# to_graphviz(SchMarkedLabelledPetriNet, graph_attrs=Dict(:dpi=>"72",:size=>"4",:ratio=>"expand"))

"""
    Given a petri net whose type `<:AbstractMarkedLabelledPetriNet`, return a `NamedTuple` giving the number
of tokens assigned to each place; this is called a marking.
"""
function marking(pn::T) where {T<:AbstractMarkedLabelledPetriNet}
    names = Tuple(pn[:,:sname])
    vals = length.(incident(pn, parts(pn,:S), :m))
    return NamedTuple{names}(vals)
end


# --------------------------------------------------------------------------------
# make the set of rewrite rules

"""
    Given an acset of type `<:AbstractMarkedLabelledPetriNet`, construct a DPO rewrite rule for transition `t`.
"""
function make_rule(pn::T, t) where {T<:AbstractMarkedLabelledPetriNet}
    input_m = inputs(pn, t)
    output_m = outputs(pn, t)

    # get L
    L = T()
    copy_parts!(L, pn, S=:)

    if length(input_m) > 0
        # add tokens
        add_parts!(
            L, :M, length(input_m),
            m = vcat(incident(L, pn[input_m, :sname], :sname)...)
        )
    end

    # get R
    R = T()
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

"""
    A presentation of a clock system, which is responsible for maintaining all the homsets for each
event possible in the system, and the set of enabled clocks for homset. It also stores the functions that
return a distribution over possible waiting times given a current simulation time.
"""
@present SchClockSystem(FreeSchema) begin
    (Clock,Event)::Ob # clock is a single ID, event is an entire class of clocks (e.g. "typed" clocks)
    NameType::AttrType # each event has a name
    RuleType::AttrType # each event has a rewrite rule
    DistType::AttrType # each event has a function that returns a firing time when it's enabled
    MatchType::AttrType # each event has a homset
    KeyType::AttrType
    event::Hom(Clock,Event)
    key::Attr(Clock,KeyType)
    name::Attr(Event,NameType)
    rule::Attr(Event,RuleType)
    dist::Attr(Event,DistType)
    match::Attr(Event,MatchType)
end

# let's look at it
# to_graphviz(SchClockSystem)

"""
    The acset type that stores instances of `SchClockSystem`
"""
@acset_type ClockSystem(SchClockSystem, index=[:event], unique_index=[:name,:key])

# --------------------------------------------------------------------------------
# set up and run a SPN model

"""
    Given a `spn<:AbstractMarkedLabelledPetriNet` and a dictionary mapping each transition name to
a function that takes in parameter `t` and returns a distribution object, sample a single trajectory
of the stochastic dynamics on the petri net until `maxevent` events occur. Print stuff (or not) for debugging
with `verbose`.
"""
function run_spn(spn::T,clockdists,maxevent,verbose=false) where {T<:AbstractMarkedLabelledPetriNet}

    # clock system stores the sampling method
    sirclock = @acset ClockSystem{Symbol,Rule,Function,Incremental.IncCCHomSet,Tuple{Int,Int,Int}} begin
        Event=nt(spn)
        name=spn[:,:tname]
    end

    # add rules, distributions
    for t in parts(sirclock,:Event)
        sirclock[t,:rule] = make_rule(spn, t)
        sirclock[t,:dist] = clockdists[sirclock[t,:name]]        
    end

    # make incremental homsets after all rules are made
    for t in parts(sirclock,:Event)
        sirclock[t,:match] = IncHomSet(codom(sirclock[t,:rule].L), getfield.(sirclock[:,:rule], :R), statepn, single=true)
    end

    # which rules are always enabled?
    always_enabled_t = findall(nparts.(codom.(getfield.(sirclock[:, :rule], :L)), :M) .== 0)

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
    output = Vector{typeof((t=0.0,marking(spn)...))}()
    push!(output, (t=0.0,marking(spn)...))

    # when and what will happen next?
    (tnow, which) = next(sampler, tnow, rng)

    while length(output) < maxevent
        !verbose || println("event $which fired at $tnow, total number of events: $(length(output))")
        event = first(which)
        update_maps = rewrite_match_maps(
            sirclock[event, :rule], 
            sirclock[event, :match][Pair(which[2:end]...)]
        )
        spn = codom(update_maps[:rh])
        push!(output, (t=tnow,marking(spn)...))

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

"""
    Given counts of `S`, `I`, `R` persons, build a `MarkedLabelledPetriNet`
that represents the SIRS model with demography.
"""
function build_sirs_pn(S,I,R)
    pn = @acset MarkedLabelledPetriNet{Symbol} begin
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
    return pn
end

# --------------------------------------------------------------------------------
# simulate a few models

# run the SIRS model with demography
statepn = build_sirs_pn(95,5,0)

# parameters to specify the random waiting times
pop = nparts(statepn, :M)
lifespan = 65*365
μ = 1/lifespan
β = 0.001
wane = 60

# functions which take a time point and return a distribution of waiting times
clockdists = Dict{Symbol,Function}()

# the Exponential clocks (Markov)
clockdists[:inf] = (t) -> Exponential(1 / β)
clockdists[:birth] = (t) -> Exponential(1 / (μ*pop))
clockdists[:deathS] = (t) -> Exponential(1 / μ)
clockdists[:deathI] = (t) -> Exponential(1 / μ)
clockdists[:deathR] = (t) -> Exponential(1 / μ)
clockdists[:wane] = (t) -> Exponential(wane)

# the Weibull clock (non-Markov)
α, θ = weibullpar(30, 5)
clockdists[:rec] = (t) -> Weibull(α,θ)

sirout = run_spn(statepn, clockdists, 2000, false)

f = Figure()
ax = Axis(f[1,1])
ln1 = lines!(ax, first.(sirout), getindex.(sirout,2))
ln2 = lines!(ax, first.(sirout), getindex.(sirout,3))
ln3 = lines!(ax, first.(sirout), getindex.(sirout,4))
Legend(f[1, 2],
    [ln1,ln2,ln3],
    ["S", "I","R"]
)
f
Makie.save("figures/SIRStrajectory.png", f, px_per_unit=1, size=(800,600))

# can also run a simple birth-death model
statepn = build_sirs_pn(100,0,0)
sirout = run_spn(statepn, clockdists, 2000, false)
f = lines(first.(sirout), getindex.(sirout,2))
Makie.save("figures/BDtrajectory.png", f, px_per_unit=1, size=(800,600))