using Catlab, AlgebraicPetri, AlgebraicRewriting

# we want a PN with a marking
@present SchMarkedLabelledPetriNet <: SchLabelledPetriNet begin
    M::Ob
    m::Hom(M,S)
end

@acset_type MarkedLabelledPetriNet(SchMarkedLabelledPetriNet, index=[:it, :is, :ot, :os]) <: AbstractLabelledPetriNet

to_graphviz(SchMarkedLabelledPetriNet)

# an SIR model with birth and death
epiPN = @acset MarkedLabelledPetriNet{Symbol} begin
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

to_graphviz(epiPN)

# add a marking
add_parts!(epiPN, :M, 500, m=[fill(1,450); fill(2,50)])

# how many S people?
length(incident(epiPN, :S, [:m,:sname]))

# make the set of rewrite rules
for t in parts(epiPN, :T)

end

t=1

# get L
input_m = epiPN[incident(epiPN, t, :it), [:is,:sname]]

L = MarkedLabelledPetriNet{Symbol}()

if length(input_m) > 0
    add_parts!(L, :S, length(unique(input_m)), sname=unique(input_m))
    add_parts!(L, :M, length(input_m), m=vcat(incident(L, input_m, :sname)...))
end

# get R
output_m = epiPN[incident(epiPN, t, :ot), [:os,:sname]]

R = MarkedLabelledPetriNet{Symbol}()

if length(output_m) > 0
    add_parts!(R, :S, length(unique(output_m)), sname=unique(output_m))
    add_parts!(R, :M, length(output_m), m=vcat(incident(R, output_m, :sname)...))
end

# get I
l, r = first(maximum_common_subobject([L,R]))