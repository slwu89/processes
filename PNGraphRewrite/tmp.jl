using Catlab, AlgebraicPetri, AlgebraicRewriting
using Distributions, Fleck

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

    I, span = only(maximum_common_subobject([L,R], abstract=false))

    retval = (
        rule = Rule(legs(first(span))[1], legs(first(span))[2]),
        L = L,
        R = R,
        I = I
    )

    return retval
end

# rewrite rules for our model
sirpn_rules = Dict([
    sirpn[t,:tname] => make_rule(sirpn, t)
    for t in parts(sirpn, :T)
])

# --------------------------------------------------------------------------------
# continuous time discrete event simulation

# how to ensure `get_dist` returns a value of type `T`?
mutable struct process{T<:ContinuousUnivariateDistribution,C,F<:Function}
    dist::T
    clocks::C
    get_dist::F
end

pp = process(Exponential(), 1, x->x)


# clocks: each event has a set of clock IDs
sirpn_clocks = Dict([
    :inf => process{Exponential,Int,}
])

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