#!/usr/bin/env sage

load("tools/init_continuous.sage")

def get_allpde(pde,constraints=[]):
	"""Returns a list of differential consequences of pde.
	pde = [{i: v1_i == ...}, {i: v2_i == ...}, ...]"""
	allpde = []
	
	def replace(func,iterations=components*dt):
		out = func
		for i in [1..iterations]:
			out = out.subs(allpde)
		return expand(out)
	
	for e in constraints:
		for index in indices(order - dt,dim):
			ediff = vdiff(e,deri(index))
			if not( ediff.lhs() in [eqn.lhs() for eqn in allpde] ):  #avoid multiple substitution
				allpde += [replace(ediff)]
	
	for component in [1..components]:
		for index in indices(order,dim):
			add = True
			### Determine lowest equation we can use to replace
			lowest = 1
			for i in reversed([3..numtimes] + [numtimes+2..dim]): 
				if index[i-1] > 0 and i in pde[component-1]:
					lowest = i
			if component == secondorder:
				if index[1] > 1 and 2 in pde[component-1]:
					lowest = 2
			else:
				if index[1] > 0 and 2 in pde[component-1]:
					lowest = 2
			### Do replacement
			if not(lowest == 1 or lowest == numtimes + 1):
				lowindex = copy(index)
				if component == secondorder and lowest == 2:
					lowindex[lowest-1] += -2
				else:
					lowindex[lowest-1] += -1
				lhs = vdiff(pde[component-1][lowest].lhs(),deri(lowindex))
				if not( lhs in [eqn.lhs() for eqn in allpde] ):  #avoid multiple substitution
					rhs = replace(vdiff(pde[component-1][lowest].rhs(),deri(lowindex)))
					allpde += [lhs == expand(cleantrig(rhs))]
		allpde = [ e.lhs() == replace(e.rhs()) for e in allpde] # catch incomplete replacements
	return allpde
	
def get_replace(pde,constraints=[]):
	"""Returns a function that eliminates time derivatives using pde and its consequences.
	pde = [{i: v1_i == ...}, {i: v2_i == ...}, ...]"""
	allpde = get_allpde(pde,constraints)
	def replace(func,iterations=components*dt):
		out = func
		for i in [1..iterations]:
			out = out.subs(allpde)
		return expand(out)
	return replace

### Check commutativity
def check_commutativity(pde):
	"""Checks that the flows of given equations commute
	pde = [{i: v1_i == ...}, {i: v2_i == ...}, ...]"""
	global warning
	for component in [1..components]:
		for i in pde[component-1]:
			for j in pde[component-1]:
				if i < j:
					if not( component == secondorder and i == 2):
						test = cleantrig(replace( vdiff(pde[component-1][i].rhs(), eval("t"+str(j)) ) - vdiff(pde[component-1][j].rhs(), eval("t"+str(i)) )))
					else:
						test = cleantrig(replace( vdiff(pde[component-1][2].rhs(), eval("t"+str(j)) ) - vdiff(pde[component-1][j].rhs(), t2,2 )))
					if not( test == 0 ):
						warning = True
						textadd("Nonzero commutator! " + str(i) + str(j))
						latexadd(expand( test ))
					else:
						textadd("PDEs for t" + str(i) + " and t" + str(j) + " commute.")
						
##################

def check(matrix,message="All entries = 0 mod eqns",verbose=True):
	"""Checks if a matrix is 0 modulo pde. If yes, prints message. If no, prints matrix modulo pde."""
	size = matrix.dimensions()[1]

	remainder = 0*copy(matrix)
	for i in [1..size]:
		for j in [1..size]:
			out = expand(replace(matrix[i-1,j-1]))
			if out == 0:
				remainder[i-1,j-1] = 0
			else:
				remainder[i-1,j-1] = cleantrig(out)
	if verbose:
		if remainder == 0*remainder:
			textadd(message)
		else:
			textadd("Nonvanishing terms!")
			latexadd(remainder)
	return remainder == 0*remainder
	

def diff_in_columns(matrix, lagnumtimes=numtimes):
	"""Takes difference of rows of matrix."""
	out = copy(matrix)
	for i in range(lagnumtimes):
		for j in range(lagnumtimes):
			if not(i==j) or double:
				if j==1: #first column
					out[i-1,j-1] = matrix[range(lagnumtimes)[1]-1,j-1] - matrix[i-1,j-1]
				else:
					out[i-1,j-1] = matrix[0,j-1] - matrix[i-1,j-1]
	return out
	
### 
@parallel
def elplane(matrix,index,component=1, lagnumtimes=numtimes):
	"""Calculates multi-time EL expressions of 1st kind and checks if they are satisfied."""
	global warning
	relevantmatrix = copy(matrix)
	for i in [1..lagnumtimes]: # eliminate irrelevant terms
		if index[i-1] > 0:
			for j in [1..lagnumtimes]:
				relevantmatrix[j-1,i-1] = 0
				if not(double):
					relevantmatrix[i-1,j-1] = 0
		if double and index[numtimes+i-1] > 0:
			for j in [1..lagnumtimes]:
				relevantmatrix[i-1,j-1] = 0
	vararray = varders(relevantmatrix,index,False,component)
	if viewELeqs and vararray != 0*vararray:
		textadd("Multi-time EL equations; - " + str(component) + " - " + str(index) + " - plane:")
		latexadd(vararray)
	if not(check(vararray,verbose=False)):
		warning = True
		textadd(str(index) + " - " + str(component) + " - plane")
		check(vararray,verbose=True)


@parallel
def eledge(matrix,index,component=1, lagnumtimes=numtimes):
	"""Calculates multi-time EL expressions of 2nd kind and checks if they are satisfied."""
	global warning
	relevantmatrix = copy(matrix)
	for i in [1..lagnumtimes]: # eliminate irrelevant terms
		if index[i-1] > 0:
			for j in [1..lagnumtimes]:
				if not(i==j) or double:
					relevantmatrix[j-1,i-1] = 0
	vararray = varders(relevantmatrix,index,True,component)
	diffvararray = replace(diff_in_columns(vararray))
	if viewELeqs and vararray != 0*vararray:
		textadd("Multi-time EL equations; - " + str(component) + " - " + str(index) + " - edge:")
		latexadd(vararray)
	if not(check(diffvararray,verbose=False)):
		warning = True
		textadd(str(index) + " - " + str(component) + " - edge")
		latexadd(vararray)
		latexadd(diffvararray)

### returns upper triangular part of a matrix
def utriang(matrix):
	"""Returns the upper-triangular part of the matrix."""
	triang = 0*matrix
	for i in [0..matrix.dimensions()[1]-1]:
		triang[i,i:] = matrix[i,i:]
	return triang
	
### check multi-time EL equations for given matrix of coefficients
def elcheck(matrix):
	"""Checks if the equations in pde satisfy the pluri-Lagrangian principle for the 2-form represented by matrix."""
	triang = utriang(matrix)
	global warning
	if elcheckdepth > -1:
		textadd("EL check...")
		for index in indices(elcheckdepth,dim):
			for component in [1..components]:
				if double:
					elplane(matrix,index,component)
				else:
					elplane(triang,index,component)
		for index in indices(elcheckdepth-1,dim):
			for component in [1..components]:
				eledge(matrix,index,component)
		if not(warning):
			textadd("EL check successful.")
