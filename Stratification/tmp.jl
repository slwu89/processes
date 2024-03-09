using AlgebraicPetri, AlgebraicPetri.TypedPetri
using Catlab


P_infectious = LabelledPetriNet(
  [:Pop],
  :infect=>((:Pop, :Pop)=>(:Pop, :Pop)),
  :disease=>(:Pop=>:Pop),
  :strata=>(:Pop=>:Pop)
)

# the SIS PN
sis_uwd = @relation (S, I) where (S::Pop, I::Pop) begin
  infect(S, I, I, I)
  disease(I, S)
end

to_graphviz(sis_uwd, box_labels=:name, junction_labels=:variable)

# 3rd argument is what the transitions will be called in P
sis_typed = oapply_typed(P_infectious, sis_uwd, [:infection, :recovery])
to_graphviz(sis_typed)

# the quarantine PN
quarantine_uwd = @relation (Q,NQ) where (Q::Pop, NQ::Pop) begin
  strata(Q,NQ) # enter quarantine
  strata(NQ,Q) # exit quarantine
end

to_graphviz(quarantine_uwd, box_labels=:name, junction_labels=:variable)

quarantine_typed = oapply_typed(P_infectious, quarantine_uwd, [:exit_Q, :enter_Q])
to_graphviz(quarantine_typed)

# add reflexive transitions to both
quarantine_typed = add_reflexives(
    quarantine_typed,
    [[:disease], [:disease, :infect]],
    P_infectious
)
to_graphviz(quarantine_typed)

sis_typed = add_reflexives(
    sis_typed,
    [[:strata], [:strata]],
    P_infectious
)
to_graphviz(sis_typed)

# take the pullback, or the product in the slice category Petri/P_infectious
sis_quarantine = typed_product([sis_typed, quarantine_typed])

to_graphviz(sis_quarantine)

# what if we now want to additionally stratify by age?

# the age PN
# age_uwd = @relation (Yng,Old) where (Yng::Pop,Old::Pop) begin
#     infect(Yng,Yng,Yng,Yng)
#     infect(Yng,Old,Yng,Old)
#     infect(Old,Yng,Old,Yng)
#     infect(Old,Old,Old,Old)
# end
# age_typed = oapply_typed(P_infectious, age_uwd, [:Yng_Yng, :Yng_Old, :Old_Yng, :Old_Old])

# easier way to make the above
age_typed = pairwise_id_typed_petri(P_infectious, :Pop, :infect, [:Yng, :Old])
age_typed = add_reflexives(
    age_typed,
    [[:strata,:disease], [:strata,:disease]],
    P_infectious
)
to_graphviz(age_typed)

sis_quarantine_age = typed_product([sis_typed, quarantine_typed, age_typed])
to_graphviz(sis_typed)
to_graphviz(quarantine_typed)
to_graphviz(age_typed)

to_graphviz(dom(sis_quarantine_age))