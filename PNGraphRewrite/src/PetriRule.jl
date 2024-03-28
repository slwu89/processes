using AlgebraicRewriting, Catlab, AlgebraicPetri
using Fleck
using Random
import AlgebraicRewriting.Incremental: IncSumHomSet
using AlgebraicRewriting.Incremental: pattern, key_vect, key_dict, IncCCHomSet

# This has been upstreamed and will be in next AlgRewriting release
IncSumHomSet(x::IncSumHomSet) = x
function IncSumHomSet(hs::IncCCHomSet) 
  pat = pattern(hs)
  kv = [[x] for x in key_vect(hs)]
  kd = Dict([k]=>v for (k,v) in key_dict(hs))
  IncSumHomSet(pat, coproduct([pat]), id(pat), [hs], kv, kd)
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


ob_name(pn::LabelledPetriNet, s::Int)::Symbol = pn[s, :sname]
ob_name(::PetriNet, s::Int)::Symbol = Symbol("S$s")

function PetriNetCSet(pn::AbstractPetriNet, args=[]; kw...)
  res = AnonACSet(BasicSchema(ob_name.(Ref(pn), parts(pn, :S)),[]))
  for (arg, s) in zip(args, parts(pn, :S))
    add_parts!(res, ob_name(pn, s), arg)  # Add tokens from Vector{Int}
  end
  for (s, arg) in pairs(kw)
    add_parts!(res, s, arg)  # Add tokens by name
  end
  res
end

"""Assumes that tokens are deleted and recreated, rather than preserved"""
function make_rule(pn::Union{PetriNet, LabelledPetriNet}, t::Int)
  L, R = LR = [PetriNetCSet(pn) for _ in 1:2]
  add_part!.(Ref(L), ob_name.(Ref(pn), pn[incident(pn, t, :it), :is]))
  add_part!.(Ref(R), ob_name.(Ref(pn), pn[incident(pn, t, :ot), :os]))
  Rule(create.(LR)...)
end

make_rules(pn) = make_rule.(Ref(pn), parts(pn, :T))


to_clocksys(spn::AbstractMarkedLabelledPetriNet, clockdists) = 
  to_clocksys(spn, spn, clockdists)

function to_clocksys(spn::T, init::ACSet, clockdists) where {T<:Union{AbstractPetriNet,AbstractMarkedLabelledPetriNet}}
  clock = @acset ClockSystem begin
    Global=1
    rng=[Random.RandomDevice()]
    sampler=[FirstToFire{ClockKeyType, Float64}()] # Fleck sampler
    Event=nt(spn)
    name=spn[:,:tname]
  end

  # add rules, distributions
  for t in parts(clock,:Event)
    clock[t,:rule] = make_rule(spn, t)
    clock[t,:dist] = clockdists[clock[t,:name]]    
  end

  # make incremental homsets after all rules are made
  for t in parts(clock, :Event)
    clock[t,:match] = IncSumHomSet(IncHomSet(codom(left(clock[t,:rule])),  
                                             right.(clock[:,:rule]), init))
  end

  # which rules are always enabled?
  aa = [t for t in parts(spn, :T) if isempty(incident(spn, t, :it))]
  add_parts!(clock, :AlwaysEnabled, length(aa); always_enabled=aa)

  for t in parts(clock,:Event)
    newkeys = [(t,k) for k in keys(clock[t,:match])]
    add_parts!(
      clock, :Clock, length(newkeys);
      key = newkeys, event = fill(t, length(newkeys))
    )
    for c in newkeys
      enable!(sampl(clock), c, clock[t, :dist](0.), 0., 0., rng(clock))
    end
  end
  
  return clock
end 
