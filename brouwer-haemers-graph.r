## brouwer-haemers-graph.r
## Unlicense (Ø) 2023 CLP, MR-Amsterdam.
## -------------------------------------
## Strongly regular graphs SRG(v, k, λ, μ).
## Wikipedia (en), Strongly regular graph.
## Paley graphs,          SRG(v = 4t + 1, k = 2t, λ = t − 1, μ = t). 
## Brouwer Haemers graph. SRG(v = 81, k = 22, λ = 1, μ = 6).  
##
## References.
## Andries E. Brouwer & Hendrik Van Maldeghem, Strongly regular graphs, Amsterdam and Ghent, March 2021.
## Ch. 7.4.4 Paley graphs, Ch.10.28 The Brouwer-Haemers graph.
## [🔗](https://homepages.cwi.nl/~aeb/math/srg/rk3/srgw.pdf)
## 
## Brouwer, A.E.; Haemers, W.H, Structure and uniqueness of the (81, 20, 1, 6) strongly regular graph, Tilburg University.
## [🔗](https://pure.uvt.nl/ws/portalfiles/portal/996688/structur.pdf)
##
##  Brouwer, A.E., Paley graphs.
## [🔗](https://www.win.tue.nl/~aeb/graphs/Paley.html)
##
## Brouwer, A.E., Brouwer-Haemers graph
## [🔗](https://www.win.tue.nl/~aeb/graphs/Brouwer-Haemers.html)
##
##  Code layout.
## (i)   Theory, modern algebra
## (ii)  Finite field operations and structures
## (iii) Creating finite fields
## (iv)  Calculating residuals
## (v)   Creating and plotting the graph
## (vi)  Statistics
## (vii) Examples   
## 
## ----------------------------------------------------------------- (i)
## Some modern algebra (Artin, Galois Theory, p28).
## q is prime. Zq, or Z/qZ, integers modulo q, %% q in R, Zq is a field (K) iff q is prime.
## https://math.stackexchange.com/questions/252046/the-paley-graph-with-9-vertices
##
## q is a prime power p^d. GF(m)[X]/P is a field iff P is a prime ideal.
## Let K[x] be the polynomials with coefficients in K.
## Let p(x) be polynomial anx^n + an-1x^n-1 ... a0.
## Let p(x)K[x] be ideal P, consisting of all polynomials with factor p(x).
## Let F = K[x]/P be the ring of cosets. A field iff p is irreducible in K.
## Let q be the order of F. Note that P is the zero-element in F.
## F is homeomorphic to K(ξ) the field of formal expressions an-1ξ^n-1 ... a0.ξ^0.
## Powers of ξ less then n are linearly independent and
## ξ^n ≡ -(an-1ξ^n-1 ... a0.ξ^0).
##  
## An element is represented in R as vector c(an-1, ... a0).
## ξ^k corresponds with position q-k in vector c(...).
## Multiplication is implemented as shift register (to the left).
## Let p(x) = x^2 + x + 2, coefficients in K = Z3, -1 ≡ 2.
## ξ^2 ≡ -1ξ - 2, or  ξ^2 ≡ 2ξ + 1 this generates the cyclic field F9:
## ξ^0 = 1, ξ1 = ξ, ξ2 = 2ξ + 1, ξ3 = 2ξ + 2, ξ4 = 2, ξ5 = 2 ξ6 = ξ, ξ7 = ξ + 2, ξ8 = ξ + 1.
##
## Irreducible, primitive polynomials over GF(3)
## https://baylor-ir.tdl.org/bitstream/handle/2104/8793/GF3%20Polynomials.pdf?sequence=1
##
## ------------------------------------------------------------------- (ii)
## Primitive polynomials over finite fields GF(q).
## p(x) = m               GF(m), m is prime, d=1.
## px  <- c()
##
## p(x) = x^2 + x + 2     GF(3^2). m=3, d=2.
## px  <- -c(1, 2)        # primitive polynomial in K.
##
## p(x) = x^2 + x + 2     GF(5^2). m=5.
## px  <- -c(1, 2)        # primitive polynomial in K.
##
## p(x) = x^2 + x + 2     GF(7^2). m=7.
## px  <- -c(1, 3)        # primitive polynomial in K.
##
## p(x) = x^2 + x + 7     GF(11^2). m=11.
## px  <- -c(1,7)         # primitive polynomial in K.
##
## p(x) = x^2 + x + 2     GF(13^2). m=13.
## px  <- -c(1,2)         # primitive polynomial in K.
##
## p(x) = x^3 + x + 2     GF(5^3). m=5, d = 3.
## px  <- -c(1,1,3)       # primitive polynomial in K.
##
## BrHa                   Brouwers Haemers graph when TRUE
## p(x) = x^4 + x + 2     GF(3^4). m=3, d=4.
## px  <- -c(0, 0, 1, 2)  # primitive polynomial in K.
## -----------------------------------------------------------------------
BrHa  <- TRUE             # Adjacent vertices differ nonzero fourth power.
m     <- 3                # K = GF(m), E = GF(m^d).
px    <- -c(0,0,1,2)      # Primitive polynomial in K. 
d     <- length(px)       # d is degree of extension.

## calculate powers of ξ^n, ξ^n+1
pwx <- list()
nxt <- rep(0, d); nxt[1] <- 1     
for (i in seq_len(d-1) ) {
   nxt      <- (c(tail(nxt, d-1), 0) + head(nxt, 1) * px) %% m
   pwx[[i]] <- nxt
}

## multiply by ξ and reduce
gf_shift <- function(x, px, q=m) { n <- length(px); (c(tail(x, n-1), 0) + head(x, 1) * px) %% q }

## multiply in quotient ring K[x]/P)
## polynom package
## https://stackoverflow.com/questions/39979884/calculate-the-product-of-two-polynomial-in-r
gf_mult <- function(a, b, px=px, q=m) {
  n  <- length(a)
  ab <- outer(a, b)
  p  <- as.vector(tapply(ab, row(ab) + col(ab), sum))

  ## reduce powers x^q, x^q+1
  for (i in seq_len(n-1)) {
    p[- (1:n-1)] <-  tail(p, n) + pwx[[i]] * p[n - i]
  }
  return(tail(p, n) %% q)
}

## ---------------------------------------- (iii)
## create Galois field of order q (GFq).
## points in GF(q) == graph vertices.
if (d > 1) {
  start <- rep(0, d); start[d] <- 1
  alpha <- c(rep(0, d-2), c(1,0))
  nxt   <- start

  GFq <- list(rep(0, d))                # Galois field including zero.
  for (i in seq_len(m^d)+1) {
    GFq[i] <- list(nxt)
    nxt <- gf_mult(nxt, alpha)          # Shift left and reduce, multiply by ξ.
    if (all(nxt == start)) break        # Full cycle iff polynomial is primitive.
  }
} else {
  GFq <- seq_len(m)-1                   # Short cut as we do not know the primitive element.
  split(GFq, seq_along(GFq))            # transpose vector
}

q <- length(GFq)                        # Order of field.
stopifnot(q > 1)                        # Non trivial.
if (q != m^d) warning("Polynomial px is not primitive.") # Warning.
if (q %% 4 != 1) warning("Not a Paley graph")            # Warning.

## ----------------------------------------------------------- (iv)
## calculate all quadratic residuals in (QR).
QR <- list()
for ( i in seq_along(GFq[]) ){
    vtx <- GFq[[i]]
    if (BrHa && q == 81) vtx <- gf_mult(vtx, vtx)  # ^2
    QR[[i]] <- gf_mult(vtx, vtx)                   # ^2
}
QR <- unique(QR)

## ------------------------------------------------------------- (v)
## Calculate edges for all vertex pairs in upper half (el)
## Determmine difference k, i.
## If difference is a nonzero quadratic residual then create edge
idx <- 0
df <- data.frame(matrix(ncol = 3, nrow = 0, dimnames=list(NULL, c("to", "from", "diff"))))
for (i in seq_len(q-1)){
  for (k in seq(i+1, q) ){
    idx      <- idx + 1
    diff     <- GFq[[k]] - GFq[[i]]                 # difference between point k, i.
    diff <- (diff + m) %% m                         # all positive.
    df[idx,] <- c(i, k, list(list(diff)) )
  }
}

qresi <- match(df[,3], QR)                          # index in residuals QR.
df[, "qr-idx"] <- qresi
el <- df[!is.na(df[,"qr-idx"]),]                    # Remove non residuals.

## plot and statistics
library(igraph)
g <- graph_from_data_frame(el, directed = FALSE)
k   <- mean(min(degree(g)), max(degree(g)))
g$layout  <- layout_in_circle(g)
g$main <- paste0(ifelse(!exists("BrHa") || BrHa, "Brouwer, Haemers", "Payley"), " graph (v = ", q, ", k = ", k, ")")
V(g)$size <- 5
V(g)$label <- NA
c(q * (q-1) / 4, ecount(g))
plot(g)

## ----------------------------------------------------------- (vi)
ev  <- round(eigen(g[])$values, digits = 6)
nev <- length(unique(ev))
t   <- (q-1) / 4
{
message(sprintf(g$main))
message(sprintf("Vertices               : %i"      , vcount(g)))
message(sprintf("Degree (k)             : %i, %i"  , 2*t/(BrHa+1), k ))
message(sprintf("Eigenvalues (r, s, #)  : %.1f, %.1f, %i", max(ev), min(ev), nev) )
message(sprintf("Diameter               : %i"      , diameter(g)) )
message(sprintf("Radius                 : %i"      , radius(g)) )
message(sprintf("Girth (shortes cycle)  : %i"      , girth(g)$girth) )
message(sprintf("Independent set (α)    : %i"      , length(largest_ivs(g)[[1]])) )
message(sprintf("Edge connectivity      : %.1f"    , edge_connectivity(g)) )
message(sprintf("Edges                  : %.1f, %i", q * t/ (BrHa+1), ecount(g) ))
message(sprintf("Automorphisms          : %.1f, %s", d*q*(q-1) / 2, automorphisms(g)$group_size) )
}
max(greedy_vertex_coloring(g))
stopifnot(min(degree(g)) == max(degree(g)))
stopifnot(max(ev) == max(degree(g)))
stopifnot(nev == 3)

## ------------------------------------------------------------------------------------------------
## https://www.win.tue.nl/~aeb/graphs/srg/srgtab.html
## SRG parameters.
##  v   - number of vertices
##  k   - valency, degree
##  λ   - number of common neighbours of two adjacent vertices (lambda)
##  μ   - number of common neighbours of two nonadjacent vertices (mu)
##  r^f - positive eigenvalue with multiplicity
##  s^f - negative eigenvalue with multiplicity
##
## Paley graph
## For q = 4t + 1, the parameters are v = 4t + 1, k = 2t, λ = t − 1, μ = t.
## q 5 9 13 17 25 29 37 41 49 53 61 73 81 89 97 101 109 113 121 125 137 149 157 169 173 181 193 197
## α 2 3  3  3  5  4  4  5  7  5  5  5  9  5  6   5   6   7  11   7   7   7   7  13   8   7   7   8
## χ 3 3  5  6  5  8 10  9  7 11 13 15  9 18 17  21  19  17  11  18  20  22  23  13  22  26  28  25
## https://www.win.tue.nl/~aeb/graphs/Paley.html
##
## Brouwer Haemers graph, v = 81, k = 20,  λ = 1, μ = 6, α = 15, χ = 6 of 7(E van Dam).
## Radius = 2, Diameter = 2, Girth = 3.
## Automorphisms    : 233,280 (2^6 * 3^6 * 5)
## https://www.win.tue.nl/~aeb/graphs/Brouwer-Haemers.html
##
## ------------------------------------------------------------- (vii)
## Manual construction of the Paley graph for primes, q = 13.
q <- 13
g13 <- graph_from_literal( 0-1-2-3-4-5-6-7-8-9-10-11-12-0    # (1)
                       , 0-3-6-9-12-2-5-8-11-1-4-7-10-0    # (3)
                       , 0-4-8-12-3-7-11-2-6-10-1-5-9-0    # (4)
     )
dev.new(); plot(g13, layout=layout_in_circle)

# Or more general for q is prime.
q <- 13; QRg <- c(1,3,4)
q <- 17; QRg <- c(1,2,3,5)
q <- 5;  QRg <- c(1)
BrHA <- FALSE
g2 <- graph_from_data_frame(data.frame(from=numeric(), to=numeric()),
                           directed = FALSE,
                           vertices = seq(q) - 1)
for (h in QRg) {
     g2 <- g2 + path(as.character(seq(from=0, to=q*h, by=h) %% q) )
}
dev.new(); plot(g2, layout=layout_in_circle)

## ----------------------------------------------------------------
## compute generators of quadratic residues mod q, q is prime (QR).
q <- 29
q1 <- (seq(q/2) * seq(q/2)) %% q
q2 <- c(q1, -q1 + q) 
QRg <- sort(unique(q2[q2 < q/2]))
