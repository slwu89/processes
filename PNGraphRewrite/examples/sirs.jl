using PNGraphRewrite
using Catlab, AlgebraicPetri
using Distributions

using Makie, CairoMakie

"""
Given counts of `S`, `I`, `R` persons, build a `MarkedLabelledPetriNet`
that represents the SIRS model with demography.
"""
build_sirs_pn(S=0,I=0,R=0) = @acset MarkedLabelledPetriNet{Symbol} begin
  S=3; sname=[:S,:I,:R]
  T=7; tname=[:inf,:rec,:birth,:deathS,:deathI,:deathR,:wane]
  I=7
  it=[1,1,2,4,5,6,7]
  is=[1,2,2,1,2,3,3]
  O=5
  ot=[1,1,2,3,7]
  os=[2,2,3,1,1]
  M=sum([S,I,R])
  m=[fill(1,S); fill(2,I); fill(3,R)]
end

sirs_unlabeled = PetriNet()
copy_parts!(sirs_unlabeled, build_sirs_pn())
t1,t2,t3,t4,t5,t6,t7=PNGraphRewrite.make_rules(sirs_unlabeled)

sir_pn= LabelledPetriNet()
copy_parts!(sir_pn, build_sirs_pn())
t1,t2,t3,t4,t5,t6,t7=PNGraphRewrite.make_rules(sir_pn)

@assert PetriNetCSet(sir_pn, [1,2,3]) == PetriNetCSet(sir_pn; I=2,R=3,S=1)

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


# ------------------------------------------------------------------------------
# simulate a few models
sirout = run_spn!(statepn, clockdists; save=marking, maxevent=2000)

X = first.(sirout)
SIR = [getindex.(last.(sirout), x) for x in [:S, :I, :R]]
f = Figure();
ax = Axis(f[1,1])
ln1, ln2, ln3 = lines!.(Ref(ax), Ref(X), SIR)
Legend(f[1, 2], [ln1,ln2,ln3], ["S", "I","R"])
f
Makie.save("figures/SIRStrajectory.png", f, px_per_unit=1, size=(800,600))

# can also run a simple birth-death model
statepn = build_sirs_pn(100,0,0)
BDout = run_spn!(statepn, clockdists; save=marking)
f = lines(first.(BDout), getindex.(last.(BDout), :S))
Makie.save("figures/BDtrajectory.png", f, px_per_unit=1, size=(800,600))

# Do computation on CSet schema with just discrete ob, rather than carrying 
# around Petri Net schema for every step

sirout = run_spn!(sir_pn, clockdists, PetriNetCSet(sir_pn, [100,5,0]); 
                 save=(X->nparts.(Ref(X),[:S,:I,:R])), maxevent=2000)
X = first.(sirout)
SIR = [getindex.(last.(sirout), x) for x in 1:3]
f = Figure();
ax = Axis(f[1,1])
ln1, ln2, ln3 = lines!.(Ref(ax), Ref(X), SIR)
Legend(f[1, 2], [ln1,ln2,ln3], ["S", "I","R"])
f
Makie.save("figures/SIRSDiscretetrajectory.png", f, px_per_unit=1, size=(800,600))
                 