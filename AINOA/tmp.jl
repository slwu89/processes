using Catlab

""" Abstract type for C-sets that contain an AINOA graph.

This type encompasses C-sets where the schema for AINOA graphs is a subcategory of C.
"""
@abstract_acset_type HasAinoa

@present SchAinoa(FreeSchema) begin
    (A,I,N,O)::Ob
    in::Hom(I,N)
    ia::Hom(I,A)
    on::Hom(O,N)
    oa::Hom(O,A)
end

@abstract_acset_type AbstractAinoa <: HasAinoa

@acset_type Ainoa(SchAinoa, index=[:in,:on], unique_index=[:ia,:oa]) <: AbstractAinoa

# to_graphviz(SchAinoa)