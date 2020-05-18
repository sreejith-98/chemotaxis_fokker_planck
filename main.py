import sympy as sy
import math

# Initialization of the variables


#velocity coeffiecients
# x_dot,v = sy.symbols('x_dot v')
# v = sy.symbols('v')

# alpha_coefficients
a1x,a2x,a11,a12,a21,a22 = sy.symbols('a1x a2x a11 a12 a21 a22')

#beta coefficients with x in it
bxx,b1x,b2x = sy.symbols('bxx b1x b2x')

#beta coefficients without x
b11,b12,b22 = sy.symbols('b11 b12 b22')

# d coeffiecients
d1,d2 = sy.symbols('d1 d2')

#omega coefficients
w1,w2 = sy.symbols('w1 w2')

#mean coeffiecients
nx,n1,n2 = sy.symbols('nx n1 n2')

#defining theta_0
t0 = sy.symbols('t0')

#defining the relationships
nx = -w2*b2x/w1
#Define the numerical

# COMPLEETTTEE THIS

v =
d1 =
d2 = 0
a12 = 0



#Solcing the mean equations and defining relations with respect to the solutions
expr1 = n1*(a11-2*w1)+n2*a12-2*w2*b12+a1x*nx
expr2 =n1*(a21-2*w2)+n2*a22-2*w2*b12+a2x*nx

aasd = sy.linsolve([expr1,expr2],(n1,n2))

n1 = sy.simplify(aasd.args[0][0])
n2 = sy.simplify(aasd.args[0][1])

t0 = -w2*n2/w1 #Thus we have relations only of the beta coefficients
bxx = -v*nx/w1


# The actual set of equations
eq1x = -v*n1+a1x*bxx+a11*b1x+a12*b2x-2*w1*b1x
eq2x = -v*n2+a2x*bxx+a21*b1x+a22*b2x-2*w1*b2x
eq11 = 2*b1x*a1x+2*a11*b11+2*(b12*a12)-2*w1*b11+2*d1*t0
eq22 = 2*b2x*a2x+2*a22*b22+2*b12*a21-2*w1*b22+2*d2*t0
eq12 = b2x*a1x+a2x*b1x+a12*b22+a21*b11+a11*b12+a12*b12-2*w1*b12

aasd = sy.nonlinsolve([eq1x,eq2x,eq11,eq22,eq12],(b1x,b2x,b11,b22,b12))

b1x = sy.simplify(aasd2.args[0][0])
b2x = sy.simplify(aasd2.args[0][1])
b11 = sy.simplify(aasd3.args[0][2])
b22 = sy.simplify(aasd3.args[0][3])
b12 = sy.simplify(aasd3.args[0][4])


b1x_f = sy.lambdify((a1x,a2x,a11,a21,a22,w1,w2),b1x)
b2x_f = sy.lambdify((a1x,a2x,a11,a21,a22,w1,w2),b2x)
b11_f = sy.lambdify((a1x,a2x,a11,a21,a22,w1,w2),b11)
b22_f = sy.lambdify((a1x,a2x,a11,a21,a22,w1,w2),b22)
b12_f = sy.lambdify((a1x,a2x,a11,a21,a22,w1,w2),b12)

#Write code to print this equations out

## Runnin the simulations


import numpy as np
#Initial conditions

t_start = 0
t_final = 100
N = 1000
Num = 6 # This is the variable that goes in to the activation
L0 =  #This gives the amplitude of the ligand field
Gamma = #This gives the range of the ligand field

alpha =  #The next four params are from the epsilon calculation
m0 =
k_a =
k_i =

k_r = #v1 params
k_b =

#Fill this initial list
x0 = 0
y_10 =
y_20 =

t = np.linspace(t_start,t_final,N)
x = np.zeros(N)
y_1 = np.zeros(N)
y_2 = np.zeros(N)
a1x_a = np.zeros(N)
a2x_a = np.zeros(N)
a11_a = np.zeros(N)
a22_a = np.zeros(N)
w1_a = np.zeros(N)
w2_a = np.zeros(N)

derv_a_y1 = np.zeros(N)
derv_a_x = np.zeros(N)
activ_a = np.zeros(N)
epsilon_a = np.zeros(N)

#Imposing the initial conditions
x[0] = x0
y_1[0] = y_10
y_2[0] = y_20

def ligand_f(pos):
    return L0*math.exp(pos/Gamma)

def epsilon_f(i):
    return alpha*(m0-y_1[i])-math.log((1+ligand_f(x[i])/k_a)/(1+ligand_f(x[i])/k_i))

def activ_f(i):
    return 1/(1+math.exp(Num*epsilon_a(i)))

def derv_a_x_f(i):
    ligand_derv_x = -((1+ligand_f(x[i])/k_a)/(1+ligand_f(x[i])/k_i))*(-(1+ligand_f(x[i])/k_a)*(1/((1+ligand_f(x[i])/k_i)**2))*(1/k_i)*ligand_f(x[i])*(1/Gamma)+(1/k_a)*(1/(1+ligand_f(x[i])/k_i))*ligand_f(x[i])*(1/Gamma))
    return -1*activ_a[i]*activ_a[i]*Num*math.exp(Num*epsilon_a[i])*ligand_derv_x

def derv_a_y1_f(i):
    return -1*activ_a[i]*activ_a[i]*Num*math.exp(Num*epsilon_a[i])*(-alpha)

def a1x_f(i):
    return -(k_r+k_b)*derv_a_x[i]

def a11_f(i):
    return -(k_r+k_b)*derv_a_y1[i]

def a2x_f(i):
    return k_y*(1-y_2[i])*derv_a_x[i]

def a21_f(i):
    return k_y*(1-y_2[i])*derv_a_y1[i]

def a22_f(i):
    return -k_z-k_y*activ_a[i]

for i in range(N):
    epsilon_a[i] = epsilon_f(i)
    activ_a[i] = activ_f(i)
    derv_a_x[i] = derv_a_x_f(i)
    derv_a_y1[i] = derv_a_y1_f(i)

    a = [a1x_f(i),a2x_f(i),a11_f(i),a21_f(i),a22_f(i),w1_f(i),w2_f(i)]

    a1x_a[i] = a[0]
    a2x_a[i] = a[1]
    a11_a[i] = a[2]
    a21_a[i] = a[3]
    a22_a[i] = a[4]
    w1_a[i] = a[5]
    w2_a[i] = a[6]

    b1x_a[i] = b1x_f(a[0],a[1],a[2],a[3],a[4],a[5],a[6])
    b2x_a[i] = b2x_f(a[0],a[1],a[2],a[3],a[4],a[5],a[6])
    b11_a[i] = b1x_f(a[0],a[1],a[2],a[3],a[4],a[5],a[6])
    b22_a[i] = b2x_f(a[0],a[1],a[2],a[3],a[4],a[5],a[6])
    b12_a[i] = b1x_f(a[0],a[1],a[2],a[3],a[4],a[5],a[6])

    x[i+1] = x_update(i)
    y_1[i+1] = y_1update(i)
    y_2[i+1] = y_2update(i)
