<!--
 * Author: Hantao Nie (nht@pku.edu.cn)
 * Date: 2023-09-11 21:33:06
 * LastEditors: Hantao Nie (nht@pku.edu.cn)
 * LastEditTime: 2023-09-18 12:26:32
 * Description: 
 * 
 * Copyright (c) 2023, Hantao Nie, Peking University. 
-->
This folder contains files for linear algebra operations, including constructing matrices, factorizing matrices, and solving linear systems. The codes in this folder are just examples. You can copy and modify them to fit your custumed algorithm.

## Construct linear system
In many algorithms, we construct the coefficient matrix of the Newton system as 
$$ A  H  A^\top$$
For different type of cones, $H$ has different forms. To handle this, we provide different strategies to construct $H$, considering the sparsity of $A$.

we split the lhs coefficien matrix into 3 * 3 subblocks as
```matlab
lhs.mat = [lhs.mat11, lhs.mat12, lhs.mat13;
           lhs.mat21, lhs.mat22, lhs.mat23;
           lhs.mat31, lhs.mat32, lhs.mat33];
```
lhs.mat11 is constructed by the main part of $A  H  A^\top$ .
lhs.mat12, lhs.mat12 is used for storing the part that are generated from Schur complement in order to handling handling free variables .
lhs.mat31, lhs.mat32, lhs.mat33 is used for storing the part that are generated from Schur complement in order to dense columns/rows. 

Ususally lhs.mat is sysmmetric, i.e., lhs,mat21 = lhs.mat12', lhs.mat31 = lhs.mat13', lhs.mat32 = lhs.mat23'. However, lhs.mat idoes not have to be symmetric. 

## factorize linear system & solve linear system
We provide a interface for factorizing and solving linear system, calling matlab's build-in functions and also external packages like pardiso. You can modify the codes to fit your custumed algorithm.