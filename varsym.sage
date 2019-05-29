#!/usr/bin/env sage

### Clear "varnames" if set
try:
	varnames
except NameError:
	0
else:
	del varnames
	
### Load input
# load('varsym_input.sage')

secondorder = 0
double = False
order = 2*max(weights)

### Initizalize variables and load methods
load('tools/2d_variational_calculus.sage')
load('tools/output.sage')

init_output("varsym-" + str(switch), even)
		
### start timing
w = walltime()

### Print mission statement
textadd( switch )

### 

warning = False

[laglist,pde,constraints] = get_input(switch)

textadd("Equations:")
for eqns in pde:
	for i in eqns:
		latexadd(eqns[i])

replace = get_replace(pde)

if checkcomm:
	check_commutativity(pde)
	for i in pde[0]:
		if not( replace(varder([1,i],laglist[i-1], field([0]*numtimes),j )) == 0):
			warning = True
			textadd("PDEs do not satisfy EL equations for L_1" + str(i))

def crossdiff(laglist):
	var('a')
	dL = matrix([[0*a for j in [0..dim-1]] for i in [0..dim-1]])
	for i in [1..numtimes-1]:
		for j in [i+1..numtimes-1]:
			dL[i,j] = replace( vdiff(laglist[j],t[i]) - vdiff(laglist[i],t[j]) )
	return dL

def matrix_integrate(m):
	F = 0*copy(m)
	for i in [0..numtimes-1]:
		for j in [0..numtimes-1]:
			F[i,j] = vintegrate(m[i,j])
	return F

def plurilag(laglist):
	dL = crossdiff(laglist)
	textadd("dL")
	latexadd(dL)

	F = matrix_integrate(dL)
	textadd("F")
	latexadd(F)
	
	L = copy(F)
	L[0] = laglist
	for i in [1..numtimes-1]:
		for j in [i+1..numtimes-1]:
			for k in [0..order-1]:
				index = [k] + [0]*(dim-1)
				index1 = [k+1] + [0]*(dim-1)
				for l in [1..components]:
					if i+1 in pde[l-1]:
						L[i,j] += varder([1,j+1],L[0][j],index1,l) * vdiff( pde[l-1][i+1].lhs() - pde[l-1][i+1].rhs(), deri(index) )
					elif switch=="SG" and i==2 and k > 0: #workaround for non-evolutionary SG equation
						L[i,j] += varder([1,j+1],L[0][j],index1,l) * vdiff( constraints[0].lhs() - constraints[0].rhs(), [k-1] + [0]*(dim-1) )
					if j+1 in pde[l-1]:
						L[i,j] += -varder([i+1,1],L[0][i],index1,l) * vdiff( pde[l-1][j+1].lhs() - pde[l-1][j+1].rhs(), deri(index) )
	return expand(L - L.transpose())

L = plurilag(laglist)

textadd("Lagrangian 2-form")
latexadd(utriang(L))

if len(constraints) > 0:
	replace = get_replace(pde,constraints)

elcheck(L)

### 

end_output(viewpdf)
if save:
	save_lagrangian("varsym-" + str(switch), [pde,utriang(L)])
