using Catlab, AlgebraicPetri, AlgebraicRewriting, AlgebraicRewriting.Incremental
using Distributions, Fleck
using Random
using SpecialFunctions
using Plots

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
# an SIR model with demography (birth and death)
sirpn = @acset MarkedLabelledPetriNet{Symbol} begin
    S=3
    sname=[:S,:I,:R]
    T=6
    tname=[:inf,:rec,:birth,:deathS,:deathI,:deathR]
    I=6
    it=[1,1,2,4,5,6]
    is=[1,2,2,1,2,3]
    O=4
    ot=[1,1,2,3]
    os=[2,2,3,1]
end

to_graphviz(sirpn)

# add a marking
add_parts!(sirpn, :M, 100, m=[fill(1,95); fill(2,5)])

# check the marking
marking(sirpn)

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

# -----------------
# test manual updates

# rewrite rules for our model
sirpn_rules = Dict([
    sirpn[t,:tname] => make_rule(sirpn, t)
    for t in parts(sirpn, :T)
])

# ---------------
# "manual" update

hset_inf = IncHomSet(codom(sirpn_rules[:inf].L), [sirpn_rules[t].R for t in sirpn[:,:tname]], sirpn)
@assert prod(size(matches(hset_inf))) == marking(sirpn).S * marking(sirpn).I

# make an infection happen
event_maps = rewrite_match_maps(sirpn_rules[:inf], sample(matches(hset_inf)))

# the new state
marking(sirpn)
marking(codom(event_maps[:kh]))

Incremental.deletion!(hset_inf, event_maps[:kg])
Incremental.addition!(hset_inf, 1, event_maps[:rh], event_maps[:kh])

prod(size(matches(hset_inf)))

# ---------------
# rewrite! update
rewrite!(hset_inf, sirpn_rules[:inf], sample(matches(hset_inf)))





# --------------------------------------------------------------------------------
# we want something to store the rules, clocks associated to each, and their type

@present SchClockSystem(FreeSchema) begin
    (Clock,Event)::Ob # clock is a single ID, event is an entire class of clocks (e.g. "typed" clocks)
    NameType::AttrType # each event has a name
    RuleType::AttrType # each event has a rewrite rule
    DistType::AttrType # each event has a function that returns a firing time when it's enabled
    MatchType::AttrType # each event has set of matches associated with it (in theory, this would instead be bijective with Clock)
    event::Hom(Clock,Event)
    name::Attr(Event,NameType)
    type::Attr(Event,NameType) # markov or not; get rid of this on the next go?
    rule::Attr(Event,RuleType)
    dist::Attr(Event,DistType)
    match::Attr(Event,MatchType)
end

# to_graphviz(SchClockSystem)

@acset_type ClockSystem(SchClockSystem, index=[:event,:type], unique_index=[:name])

sirclock = @acset ClockSystem{Symbol,Rule,Function,Incremental.IncSumHomSet} begin
    Event=nt(sirpn)
    name=sirpn[:,:tname]
    type=[:Markov,:NonMarkov,:Markov,:Markov,:Markov,:Markov]
end

for t in parts(sirclock,:Event)
    sirclock[t,:rule] = make_rule(sirpn, t)
end

# --------------------------------------------------------------------------------
# add the distribution functions

pop = nparts(sirpn, :M)
lifespan = 65*365
μ = 1/lifespan

β = 0.05*5

# should somehow check that if this fn returns Exponential, its mapped to type Markov
sirclock[only(incident(sirclock, :inf, :name)), :dist] = (t) -> Exponential(1 / β)
sirclock[only(incident(sirclock, :birth, :name)), :dist] = (t) -> Exponential(1 / (μ*pop))
sirclock[only(incident(sirclock, :deathS, :name)), :dist] = (t) -> Exponential(1 / μ)
sirclock[only(incident(sirclock, :deathI, :name)), :dist] = (t) -> Exponential(1 / μ)
sirclock[only(incident(sirclock, :deathR, :name)), :dist] = (t) -> Exponential(1 / μ)

# get parameters of a Weibull (from https://github.com/cran/mixdist/blob/master/R/weibullpar.R)
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

α, θ = weibullpar(7, 2.6)

# the non Markovian clock
sirclock[only(incident(sirclock, :rec, :name)), :dist] = (t) -> Weibull(α,θ)


# --------------------------------------------------------------------------------
# add the match sets

for t in parts(sirclock,:Event)
    sirclock[t,:match] = IncHomSet(codom(sirclock[t,:rule].L), getfield.(sirclock[:,:rule], :R), sirpn)
end


# --------------------------------------------------------------------------------
# sample and add the clocks

# Fleck sampler
sampler = FirstToFire{Int, Float64}()
rng = Random.RandomDevice()
tnow::Float64 = 0.0

# add clocks
# we should actually start the recovery clocks some random time in the past, unless we assume they all got sick right at tnow, which might
# be fine for some scenarios.
for t in parts(sirclock,:Event)
    nmatch = prod(size(matches(sirclock[t,:match])))
    newclocks = add_parts!(
        sirclock, :Clock, nmatch,
        event = fill(t, nmatch)
    )
    for c in newclocks
        enable!(sampler, c, sirclock[t, :dist](tnow), tnow, tnow, rng)
    end
end

# when and what will happen next?
(tnow, which) = next(sampler, tnow, rng)

# step! should be here and act on the clock system
# sirclock[which, [:event,:rule]]

update_maps = rewrite_match_maps(sirclock[which, [:event,:rule]], sample(matches(sirclock[which, [:event,:match]])))

sirpn = codom(update_maps[:rh])

# update matches for all events
for t in parts(sirclock, :Event)
    Incremental.deletion!(sirclocl[t,:match], update_maps[:kg])
    Incremental.addition!(sirclocl[t,:match], sirclock[which, :event], update_maps[:rh], update_maps[:kh])
end



# step! ends

# -----------------------------------------
# check how matches and rewriting will work

# inf_match = get_matches(sirclock[1,:rule], sirpn)
# rec_match = get_matches(sirclock[2,:rule], sirpn)

hset_inf = IncHomSet(codom(sirclock[1,:rule].L), getfield.(sirclock[:,:rule],:R), sirpn)
hset_rec = IncHomSet(codom(sirclock[2,:rule].L), getfield.(sirclock[:,:rule],:R), sirpn)

# make an infection happen
event_maps = rewrite_match_maps(sirclock[1,:rule], sample(matches(hset_inf)))

# the new state
marking(codom(event_maps[:kh]))
sirpn = codom(event_maps[:kh])

Incremental.deletion!(hset_inf, event_maps[:kg])
Incremental.deletion!(hset_rec, event_maps[:kg])
Incremental.addition!(hset_inf, 1, event_maps[:rh], event_maps[:kh])
Incremental.addition!(hset_rec, 1, event_maps[:rh], event_maps[:kh])

rewrite!(hset_inf, sirclock[1,:rule], sample(matches(hset_inf)))e