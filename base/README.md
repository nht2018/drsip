# Framework for conic programming

## basic data structure

### MatCell

MatCell class is used store data in blocks. This data structure is used for storing data matrix
(e.g., data matrix in SDP) or function handles (e.g. Hessian in SDP)
in conic programming

MatCell can be used as cell equiped with blockwise operators, e.g.,+ - .* and so on.

### BasicCone

BasicCone class is used to store one single cone.

a cell of BasicCone can be used to describe the cartesian product of numbers of cones

### Cone

Cone is used to represent concatenation of several BasisCones

## Expample

### standard conic linear programming

![](figs/eqn_conicLP.svg)

where
![](figs/eqn_K.svg)
Here the superscript $s, q, l, u$ are the leading letters of "semidefinite", "quadratic", "linear", "unbounded" respectively.

For this problem, we descibe the cone $K$ using a cell of BasicCone objects. Example code is as follows:

```matlab
Ks = BasicCone( 's', [1, 2, 3] )  % three semidefinite cones of size 1, 2, 3 respectively
Kq = BasicCone( 'q', [4, 5, 6, 7] )  % four quadratic cones of size 4, 5, 6, 7 respectively
Kl = BasicCone( 'l', [8] )  % nonnegative orthant of size 8
Ku = BasicCone( 'u', [9] )   % unbounded variables of size 9
K = Cone({Ks, Kq, Kl, Ku} )  % concatenate Ks, Kq, Kl, Ku
```

The data matrix $\mathcal{A}$ is stored in a MatCell object, i.e.,

```matlab
A = MatCell( {As, Aq, Al, Au } )  
```

where

$A^s$ is of size $m \times \sum_{j=1}^{n_s} s_j (s_j+1)/2$,
$A^q$ is of size $m \times \sum_{i=1}^{n_q} q_i$,
$A^l$ is of size $m \times n_l$,
$A^u$ is of size $m \times n_u$.

The data vector $c$ is stored in a MatCell object, i.e.,

```matlab
c = MatCell( {cs, cq, cl, cu } )  
```

where

$c^s$ is of size $\sum_{j=1}^{n_s} s_j (s_j+1)/2$,
$c^q$ is of size $\sum_{i=1}^{n_q} q_i$,
$c^l$ is of size $n_l$,
$c^u$ is of size $n_u$.

The vector $b$ is stored in a doube array of size $m$.

The variable $x$ is stored in a MatCell object, i.e.,

```matlab
x = MatCell( { xs, xq, xl, xu } )  
```
