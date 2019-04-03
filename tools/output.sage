#!/usr/bin/env sage

import os 

### Open output files
def init_output(filename,even):
	global output
	global warning
	global plaindoc
	global latexdoc
	if not os.path.exists('logs/'):
		os.makedirs('logs/')
	filename += '-' + str(numtimes)
	if even:
		filename += 'full'
	plaindoc = open('logs/' + filename + '-plain','w')
	latexdoc = open('logs/' + filename + '-tex','w')
	
	output = []	
	warning = False

### Customized latex function to use to export tex file
### For immediate display, the built in latex() should be used
def cleanlatex(eqn):
	string = latex(eqn)
	for power in [1..9]:
		string = string.replace('^{'+str(power)+'}','^'+str(power))
	string = string.replace('{v_{','v_{') 
	string = string.replace('{\\bar v_{','\\bar v_{') 
	string = string.replace('}}','}') 
	string = string.replace('\\, ','') 
	string = string.replace('\\left(v_{}\\right)','(v)') 
	return string

### print to latex and plaintext files
### if all==True, print to console and pdf as well

def latexadd(eqn,all=True):
	global output
	print ""
	print eqn
	print ""
	if all:
		output += ['\\scalebox{.1}{$' + latex(eqn) + '$}'] #['\\tiny ' + latex(eqn)]
	if not(latexdoc.closed):
		latexdoc.write(cleanlatex(eqn) + '\n\n')
	if not(plaindoc.closed):
		plaindoc.write(str(eqn) + '\n\n')
		
def textadd(string,all=True):
	global output
	print "%.2f" % (walltime(w)) + "s: " + string
	if all:
		output += [latex('') + '\\scalebox{.1}{' + string + '}']
	if not(latexdoc.closed):
		latexdoc.write(string + '\n\n')
	if not(plaindoc.closed):
		plaindoc.write(string + '\n\n')
		

### Finalize output
def end_output(viewpdf=True):	
	if warning:
		textadd('!!! WARNING(S) GENERATED !!!')
	plaindoc.close()
	latexdoc.close()

	if viewpdf:
		view(output)

def save_lagrangian(filename, data):
	if not os.path.exists('lagrangians/'):
		os.makedirs('lagrangians/')
		
	lagfile = open('lagrangians/' + filename + '-plain','w')
	for line in data:
		lagfile.write(str(line))
		lagfile.write('\n\n')
	lagfile.close()
	
	lagfile = open('lagrangians/' + filename + '-tex','w')
	for line in data:
		lagfile.write(cleanlatex(line))
		lagfile.write('\n\n')
	lagfile.close()

