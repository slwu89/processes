using Catlab
import Catlab: vertices, edges

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

# some simple accessors
vertices(g::HasAinoa) = begin
    parts(g, :N)
end

edges(g::HasAinoa) = begin
    parts(g, :A)
end

@present SchLabeledAinoa <: SchAinoa begin
    Label::AttrType
    nlabel::Attr(N,Label)
    alabel::Attr(A,Label)
end

# to_graphviz(SchLabeledAinoa)

@abstract_acset_type AbstractLabeledAinoa <: AbstractAinoa

@acset_type LabeledAinoa(SchLabeledAinoa, index=[:in,:on], unique_index=[:ia,:oa]) <: AbstractLabeledAinoa

# from ex 1.2 in Kock's paper
gr_12 = @acset LabeledAinoa{Symbol} begin
    A=5
    alabel=[:a,:b,:c,:d,:e]
    N=3
    nlabel=[:x,:y,:z]
    I=3
    ia=[2,3,4]
    in=[2,3,3]
    O=3
    oa=[3,4,5]
    on=[2,2,2]
end

# in(G) is complement of O>->A
inboundary(g::HasAinoa) = begin
    setdiff(parts(g,:A), g[:,:oa])
end

# Out(G) is completement of I>->A
outboundary(g::HasAinoa) = begin
    setdiff(parts(g,:A), g[:,:ia])
end

# edges that belong to in and out boundaries are isolated
isolated(g::HasAinoa) = begin
    intersect(inboundary(g), outboundary(g))
end

# inner edges are incoming and outgoing from some node
# so the intersection of the two maps O>->A<-<I
inner(g::HasAinoa) = begin
    intersect(g[:,:oa], g[:,:ia])
end

inboundary(gr_12)
outboundary(gr_12)
isolated(gr_12)
inner(gr_12)

# how to plot this thing?
pg = PropertyGraph{Any}()

for v in vertices(gr_12)
    add_vertex!(pg, label=string(gr_12[v,:nlabel]))
end

inner_edges = inner(gr_12)
isolated_edges = isolated(gr_12)

out_edges = outboundary(gr_12) # edges that are not incoming to anything
setdiff!(out_edges, isolated_edges)
setdiff!(out_edges, inner_edges)

in_edges = inboundary(gr_12) # edges that are not outgoing from anything
setdiff!(in_edges, isolated_edges)
setdiff!(in_edges, inner_edges)

for a in inner_edges
    add_edge!(
        pg,
        only(gr_12[incident(gr_12, a, :oa), :on]),
        only(gr_12[incident(gr_12, a, :ia), :in]),
        label = string(gr_12[a, :alabel])
    )
end

for a in out_edges
    tgt_invis = add_vertex!(pg, style="invis", shape="none", label="")
    add_edge!(
        pg,
        only(gr_12[incident(gr_12, a, :oa), :on]),
        tgt_invis,
        label = string(gr_12[a, :alabel])
    )
end

for a in in_edges
    src_invis = add_vertex!(pg, style="invis", shape="none", label="")
    add_edge!(
        pg,
        src_invis,
        only(gr_12[incident(gr_12, a, :ia), :in]),
        label = string(gr_12[a, :alabel])
    )
end

for a in isolated_edges
    src_invis = add_vertex!(pg, style="invis", shape="none", label="")
    tgt_invis = add_vertex!(pg, style="invis", shape="none", label="")
    add_edge!(
        pg,
        src_invis,
        tgt_invis,
        label = string(gr_12[a, :alabel])
    )
end

to_graphviz(pg)