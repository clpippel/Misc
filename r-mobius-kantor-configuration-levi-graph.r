## r-mobius-kantor-configuration-levi-graph.r
## # Unlicense (Ø) 2022 CLP, MR-Amsterdam.
## ----------------------------------------------------------------------------
## Levi graph of the Möbius-Kantor configuration.
## SELF-DUAL CONFIGURATIONS AND REGULAR GRAPHS.
## Author: H. S. M. Coxeter.
## Journal: Bull. Amer. Math. Soc. 56 (1950), 413-455.
## DOI: https://doi.org/10.1090/S0002-9904-1950-09407-5 (free access) .
##
## Galois field GF(3^2) with modulus the irreducible quadratic p(x) = x^2 + x + 2.
## k.x^2 == -k.x - k.2, k in Z3 
## Successive powers of the primitive root λ = 10 in GF(3, 2):
##   (λ1) primitive root            λ = 10, shift left, overflow = 1.
##   (λ2) λ*10 = 1 | (0,0) + (-1, -2) = 21, shift left, overflow = 2.
##   (λ3) λ*21 = 2 | (1,0) + (-2, -4) = 22, shift left, overflow = 2.
##   (λ4) λ*22 = 2 | (2,0) + (-2, -4) = 02, shift left, overflow = 0.
##   (λ5) λ*02 = 0 | (2,0) + ( 0,  0) = 20, shift left, overflow = 2.
##   (λ6) λ*20 = 2 | (0,0) + (-2, -4) = 12, shift left, overflow = 1.
##   (λ7) λ*12 = 1 | (2,0) + (-1, -2) = 11, shift left, overflow = 1.
##   (λ8) λ*11 = 1 | (1,0) + (-1, -2) = 01
##
## λ1 .. λ8  = 10, 21, 22, 02, 20, 12, 11, 01 (coordinates).
## Note that for example 21 = 2x + 1 in GF(3, 2).
## Three points are collinear when their sum equals 0x + 0 in GF(3,2).
## This is a consequence of the identity λr(1 + λ + λ3) = 0, which holds since
##   λ3 + λ + 1 = (λ + 2)(λ2 + λ + 2), (p. 429).
## For example: (1, 2, 4) = 10 + 21 + 02 = 00 (mod 3).
##
## We have a system of eight triples 124, 235, 346, 457, 568, 671, 782, 813,
## derived from one another by cyclic permutation of the digits.
## Möbius-Kantor configuration in RxR:
##        6
##     7     3
##      8   5
## 1      2      4

###############################################################################
require(igraph)
## build Möbius-Kantor configuration. Euler path + two edges
mkc <- make_graph(~ 1-2-5-3-4-2-8-3-6-5-7-6-8-7-1, 1-8, 4-5)
mkc$main <- "Möbius-Kantor configuration\n Vertex ids from K[x] / p(x)\n p(x) = x^2 + x + 2, K = Z3"

## set edge attributes to curved, solid and dotted (line type lty)
E(mkc)$curved <- 0
mkc <- set_edge_attr(mkc, name="curved", E(mkc)[9,14], c(-1.5, 1.5))
mkc <- set_edge_attr(mkc, name="lty"   , E(mkc), 1L )
mkc <- set_edge_attr(mkc, name="lty"   , E(mkc)[9,14], c(3L,3L) )

q3 <- sqrt(3)
ltmkc <- c( 0, 1, 1.25, 1.50 , 2, .75 , 1, .50
          , 0, 0, q3/4, 1    , 0, q3/4, 2, 1  
          ) 
dim(ltmkc) <- c(8, 2)
plot(mkc, layout=ltmkc, edge.width=4, vertex.size=10, sub="Fig.1", vertex.color="green")


###############################################################################
mkconf <- list( c(1,2,4) , c(2,3,5) , c(3,4,6) , c(4,5,7)
              , c(5,6,8) , c(6,7,1) , c(7,8,2) , c(8,1,3)
              )

## compute edges of the Möbius-Kantor graph as sequence of vertex pairs.
## levi graph construction.
## append vertex pair (i, ll) if point i is incident with line ll.
vp <- c()
for ( i in seq(length(mkconf))){
  for (ll in mkconf) {
    if (!all(is.na(match(ll, i))) ) {
      vp <- append(vp, c(i, paste(ll, collapse="")))
    }
  }
}

## build Levi graph.
## and reorder to match Galois sequence.
mkl <- make_graph(vp, directed=FALSE)
mkl$main <- "Möbius-Kantor graph\nLevi graph construction"
degree(mkl)
## vertex <-> index
## 1  124  671  813  2 235  782   3 346   4 457   5 568   6   7   8 
## 2  1    12   16   3   4   14   5 6     7   8   9  10  11  13  15
vn <- V(mkl)$name                      # vertex names.
pp <- match(vn, sort(vn))
mkg <- permute(mkl, pp)

## Levi / incidence graph of the Möbius–Kantor configuration.
dev.new(); plot(mkg, layout=layout_in_circle, edge.width=2, sub="Fig. 2")

degree(mkg)
## plot coordinates
##        001 124 002 235 003 346 004 457 005 568 006 671 007 782 008 813
ltmkg <- c(10,  8,  3,1.5,  5,  5,  3,1.5,  0,  2,  7,8.5,  5,  5,  7,8.5,
            5,  5,  7,8.5, 10,  8,  3,1.5,  4,  5,  3,1.5,  0,  2,  7,8.5
          )
dim(ltmkg) <- c(16, 2)
dev.new(); plot(mkg, layout=ltmkg, edge.width=2, main = "Möbius-Kantor graph\n{8} + {8/3} ", sub="Fig. 3")


###############################################################################
# manual construction Möbius–Kantor graph ({8} + {8/3}).
mkd <-
  make_graph(~ 1, 124, 2, 235, 3, 346, 4, 457
             , 5, 568, 6, 671, 7, 782, 8, 813         # vertex ids.
             , 1-813-3-235-5-457-7-671-1              # outer octagon.
             , 124-2-782-8-568-6-346-4-124            # inner octagon.
             , 1-124, 813-8, 3-346, 235-2         
             , 5-568, 457-4, 7-782, 6-671             # connections.
             )
mkd$main <- "Möbius-Kantor graph\nManual construction"
degree(mkd)
dev.new(); plot(mkd, layout=layout_in_circle, edge.width =2, sub="Fig. 4")


###############################################################################
# graph as a Swiss tournament
vps <- c( 1,124, 2,235, 3,346, 4,457, 5,568 ,6,671, 7,782, 8,813 # round 1.
        , 124,2, 235,3, 346,4, 457,5, 568,6, 671,7, 782,8, 813,1 # round 2.
        , 1,671, 2,782, 3,813, 4,124, 5,235, 6,346, 7, 457,8,568 # round 3.
        )
mks <- make_graph(as.character(vps), directed=FALSE)
mks$main <- "Möbius-Kantor graph\nConstructed as Swiss tournament"
E(mks)[ 1: 8]$color <- 'red'
E(mks)[ 9:16]$color <- 'blue'
E(mks)[17:24]$color <- 'green'
dev.new(); plot(mks, layout=layout_in_circle, edge.width=2, sub="Fig. 5")
dev.new(); plot(mks, layout=ltmkg, edge.width=2, sub="Fig. 6")
