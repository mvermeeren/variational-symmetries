#!/usr/bin/env sage

### SELECT LATTICE EQUATION ###

switch = 'KdV'
# Options: lin,KdV, NLS, SG

### DIMENSION OF MULTI-TIME

numtimes = 3

# Include even-numbered times?
even = True #DO NOT CHANGE

### CHECKS

# Verify commutativity of input equations?
# (Also: check if input equations solve EL for given Lagrangians)
checkcomm = True

# To which order to check multi-time EL equations
elcheckdepth = numtimes

### REQUESTED OUTPUT

# Show output as pdf?
viewpdf = True

# Print EL equations even if they are satisfied?
viewELeqs = False

# Save logs
save = False

### SET EQUATION DEPENDENT PARAMETERS ###

# weights = [order of eqn in t_i for i in [1..numtimes]]

if switch == 'KdV' or switch == 'lin':
	weights = [2*i-1 for i in [1..numtimes]]
	components = 1
if switch == 'NLS':
	weights = [1..numtimes]
	components = 2
	varnames = ["v","\\bar v"]
if switch == 'SG':
	weights = [1,1] + [2*i-1 for i in [2..numtimes-1]]
	components = 1
	varnames = ["v","\\bar v"]

#########################################

### LIST OF INPUT LAGRANGIANS ###

#laglist = [L_11, L_12, ...]
#pde = [{i: v_i == ...} for each component]
#constraints = [equations of non-evolutionary type]

def get_input(switch):
    constraints = []
    ###
    ###
    if switch == 'lin':
        if numtimes == 3:
            pde = [{2: v_2 == v_111, 3: v_3 == v_11111}]
            laglist = [0, v_1*v_2 - v_1*v_111, v_1*v_3 - v_1*v_11111]
        elif numtimes == 4:
            pde = [{2: v_2 == v_111, 3: v_3 == v_11111, 4: v_4 == v_1111111}]
            laglist = [0, v_1*v_2 - v_1*v_111, v_1*v_3 - v_1*v_11111, v_1*v_4 - v_1*v_1111111]
        else:
            pde = [{}]
            laglist = [0]
    ###
    ###
    if switch == 'KdV':
        laglist = [0, 1/2*v_1*v_2 - v_1^3 - 1/2*v_1*v_111, 
            1/2*v_1*v_3 - 5/2*v_1^4 + 5*v_1*v_11^2 - 1/2*v_111^2]
        pde = [{2: v_2 == 3*v_1^2 + v_111, 3: v_3 == 10*v_1^3 + 5*v_11^2 + 10*v_1*v_111 + v_11111}]
    ###
    ###
    if switch == 'NLS':
        k = 1
        laglist = [0, I/2*(v2_*v1_2 - v1_*v2_2) - v2_1*v1_1 - k*(v1_*v2_)^2,
            I/2*(v2_*v1_3 - v1_*v2_3) + I/2*(v2_1*v1_11 - v1_1*v2_11) - 3/2*I*k*(v1_*v2_)*(v2_1*v1_ - v1_1*v2_)]
        pde = [{2: v1_2 == I*v1_11 - 2*I*k*v1_*v2_*v1_, 3: v1_3 == v1_111 - 6*k*v1_*v2_*v1_1},
            {2: v2_2 == -I*v2_11 + 2*I*k*v2_*v1_*v2_, 3: v2_3 == v2_111 - 6*k*v2_*v1_*v2_1}]
    ###
    ###
    if switch == 'SG':
        laglist = [0, 1/2*v_1*v_2 - cos(v_), 1/2*v_1*v_3 - 1/8*v_1^4 + 1/2*v_11^2, 1/2*v_1*v_4 - 1/16*v_1^6 - 5/12*v_1^3*v_111 - 1/2*v_111^2]
        if numtimes == 3:
            pde = [{3: v_3 == v_111 + 1/2*v_1^3}]
        elif numtimes == 4:
            pde = [{3: v_3 == v_111 + 1/2*v_1^3, 4: v_4 == v_11111 + 5/2*v_1^2*v_111 + 5/2*v_1*v_11^2 + 3/8*v_1^5}]
        else:
            pde = [{}]
        constraints = [v_12 == sin(v_)]
    ###
    ###
    return [laglist,pde,constraints]
