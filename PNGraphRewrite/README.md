# Petri Nets and Graph Rewriting

It's a good idea to be able to formalize the "token game" that occurs on a marked Petri Net (henceforth, PN). This was first considered by [Kreowski (1981)](https://doi.org/10.1007/3-540-10291-4_22), but his graphs had the relations in an inconvenient direction (edges from places to tokens). [Genrich et al. (1983)](https://doi.org/10.1007/BFb0000102) fixed this. The reason is that the first way is not a function, we need to introduce another set in order to handle the one-to-many relationship (assuming we are working in **Set**), yet the second is a functional relation (a token tells you which place it belongs to).

# Notes

  * after Fleck is registered, move the dep to the one in the general repo

# References

  1. Kreowski, HJ. (1981). A comparison between petri-nets and graph grammars. In: Noltemeier, H. (eds) Graphtheoretic Concepts in Computer Science. WG 1980. Lecture Notes in Computer Science, vol 100. Springer, Berlin, Heidelberg. https://doi.org/10.1007/3-540-10291-4_22
  2. Genrich, H.J., Janssens, D., Rozenberg, G., Thiagarajan, P.S. (1983). Petri nets and their relation to graph grammars. In: Ehrig, H., Nagl, M., Rozenberg, G. (eds) Graph-Grammars and Their Application to Computer Science. Graph Grammars 1982. Lecture Notes in Computer Science, vol 153. Springer, Berlin, Heidelberg. https://doi.org/10.1007/BFb0000102