#!/usr/bin/env sage

### the discrete field u_lattice
def ul(n,m,component=1):
	"""Interpolating field of the lattice variable U at coordinates (n,m)."""
	times = [eval('t' + str(i)) for i in [1..dim]]
	for i in range(dim):
		if double:
			if i <= numtimes:
				times[i-1] += miwaconst * (-1)^(weights[i-1]+1)* (n/weights[i-1] * a^weights[i-1])
			else:
				times[i-1] += miwaconst * (-1)^(weights[i-1]+1)* (m/weights[i-1] * b^weights[i-1])
		else:
			times[i-1] += miwaconst * (-1)^(i+1)* (n/i * a^i + m/i * b^i)
	out = 'u('
	for i in [1..dim]:
		out += 'times[' + str(i-1) + '],'
	out = eval(out + str(component) + ')')
	return out
	
	
var('a', latex_name='\\alpha')
var('b', latex_name='\\beta')

def killpowers(f,truncorder=numtimes):
	"""Removes all terms with a^k or b^k for k > truncorder."""
	return expand(f).subs([a^(truncorder+i)==0 for i in [1..2*numtimes]]).subs([b^(truncorder+i)==0 for i in [1..2*numtimes]])
