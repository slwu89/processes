using Catlab, AlgebraicPetri, AlgebraicRewriting
using Distributions, Fleck
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
# we want something to store the rules, clocks associated to each, and their type

@present SchClockSystem(FreeSchema) begin
    (Clock,Event)::Ob
    NameType::AttrType
    RuleType::AttrType
    ClockType::AttrType
    event::Hom(Clock,Event)
    name::Attr(Event,NameType)
    type::Attr(Event,NameType)
    rule::Attr(Event,RuleType)
    clock::Attr(Event,ClockType)
end

# to_graphviz(SchClockSystem)

@acset_type ClockSystem(SchClockSystem, index=[:event,:type,:name])

sirclock = @acset ClockSystem{Symbol,Rule,Function} begin
    Event=nt(sirpn)
    name=sirpn[:,:tname]
    type=[:Markov,:NonMarkov,:Markov,:Markov,:Markov,:Markov]
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

    I, span = only(maximum_common_subobject([L,R], abstract=false))

    # retval = (
    #     rule = Rule(legs(first(span))[1], legs(first(span))[2]),
    #     L = L,
    #     R = R,
    #     I = I
    # )

    return Rule(legs(first(span))[1], legs(first(span))[2])
end

# # rewrite rules for our model
# sirpn_rules = Dict([
#     sirpn[t,:tname] => make_rule(sirpn, t)
#     for t in parts(sirpn, :T)
# ])

for t in parts(sirclock,:Event)
    sirclock[t,:rule] = make_rule(sirpn, t)
end

# --------------------------------------------------------------------------------
# continuous time discrete event simulation

pop = nparts(sirpn, :M)
lifespan = 65*365
μ = 1/lifespan

β = 0.05*5

# should somehow check that if this fn returns Exponential, its mapped to type Markov
sirclock[only(incident(sirclock, :inf, :name)), :clock] = (t) -> Exponential(1 / β)
sirclock[only(incident(sirclock, :birth, :name)), :clock] = (t) -> Exponential(1 / (μ*pop))
sirclock[only(incident(sirclock, :deathS, :name)), :clock] = (t) -> Exponential(1 / μ)
sirclock[only(incident(sirclock, :deathI, :name)), :clock] = (t) -> Exponential(1 / μ)
sirclock[only(incident(sirclock, :deathR, :name)), :clock] = (t) -> Exponential(1 / μ)

# the non Markovian clock
sirclock[only(incident(sirclock, :rec, :name)), :clock] = (t) -> Exponential(1 / μ)


# # test rules
# epiPN1 = deepcopy(epiPN);
# poptable(epiPN1)
# mset = get_matches(epiPN_rules[:birth].rule, epiPN1);
# epiPN1 = rewrite_match(epiPN_rules[:birth].rule, sample(mset));
# poptable(epiPN1)

# @inline function rate_to_proportion(r::T, t::T) where {T<:Float64}
#     1-exp(-r*t)
# end

# δt = 0.1
# nsteps = 400
# tmax = nsteps*δt
# tspan = (0.0,nsteps)
# t = 0.0:δt:tmax;
# p = [0.05,10.0,0.25,δt]; # β,c,γ,δt


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
mean(Weibull(α,θ))
std(Weibull(α,θ))
histogram(rand(Weibull(α,θ),1000))