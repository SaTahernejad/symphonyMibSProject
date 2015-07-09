# symphony

MibS README
===========================
MibS solves the general integer bilevel linear problems and interdiction problems (which is a special case of integer bilevel problems).

Solving General integer bilevel linear problems
===========================
For solving general integer bilevel linear problems, MibS requires two files. The first one, (which is a mps file) includes the data of the relaxed problem. The second one, includes the data of lower level problem. The format of this file is:

N : The number of lower level variables

M : The number of lower level variables

LC : The index of lower level variables

LR : The index of lower level constraints

LO : The coefficients of lower level objective function

OS : Lower level objective sense

Solving Interdiction problems
===========================
For solving the interdiction problems, MibS requires two files. The first one, (which is a mps file), includes the data of lower level problem. The format of the second file is:

N : The number of lower level variables

M : The number of lower level constraints

LC : The index of lower level variables

LR : The index of lower level constraints

LO : The coefficients of lower level objective function

OS : Lower level objective sense

IC : The coefficients of budget constraint

IB : Interdiction budget

the following command would solve the instance "sample.mps"
and “sample.txt”:

./mibs -F sample.mps -Z sample.txt
