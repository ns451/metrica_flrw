# Relatividad General | 2020-02-04 | t.me/ns451

import os
import codecs
from sympy import *

init_printing(pretty_print=True, wrap_line=False)

half = Rational(1/2)

# variables
k,ρ,p,Λ,G,c = symbols('k ρ p Λ G c')

# coordenadas
t,r,θ,φ = symbols('t r θ φ ')
ξ = [t,r,θ,φ]

# funciones
a = Function('a')(t)

echo = [] # para exportar

# línea de separación
def line(e):
    print('■'*30+' '+e)

def to_Latex(s): # left, long equations: word wrap
    begin = '\\longeq\\hspace*{-3em}$\\displaystyle'
    end = '$\\longeqend'
    echo.append('{}{}\n{}'.format(begin,s,end))

def to_Latex2(s): # por si requiere scale!=100%
    begin = '\\begin{equation}\\scalemath{1}{'
    end = '}\\end{equation}'
    echo.append('{}{}\n{}'.format(begin,s,end))

# muestra símbolo de Christoffel
def show_Γ(e,μ,ν,σ,hide_zeros):
    s = simplify(e)
    #s=e
    if str(s)!='0' and hide_zeros:
        print('\n{} Γ^[{}]_[{}][{}] =\n'.format('■'*5,μ,ν,σ))
        #pprint(s)
        a = '\n\Gamma^{}_{}={}'
        sup = '{'+str(μ)+'}'
        sub = '{'+str(ν)+str(σ)+'}'
        to_Latex(a.format(sup,sub,latex(s)))

# muestra componenente de tensor de Ricci
def show_R(e,α,γ,hide_zeros):
    s = simplify(e)
    #s=e
    if str(s)!='0' and hide_zeros:
        print('\n{} R_[{}][{}] =\n'.format('■'*5,α,γ))
        #pprint(s)
        a = '\nR_{}={}'
        sub = '{'+str(α)+str(γ)+'}'
        to_Latex(a.format(sub,latex(s)))

# muestra tensor, matriz, escalar
def show(e,v):
    print('\n{} {} ='.format('■'*15,v))
    #s = factor(e)
    s = simplify(e)
    pprint(s)
    to_Latex2('\n'+v+'='+latex(s))

# derivada covariante
def dg(α,β,γ):
    return diff(g[α,β], ξ[γ], 1)

# métrica covariante
g = zeros(4,4)
g[0,0] = -1
g[1,1] = a**2/(1-k*r**2)
g[2,2] = a**2*r**2
g[3,3] = a**2*r**2*(sin(θ))**2
show(g,'g_{\\mu\\nu}')

# métrica contravariante
gc = g.inv()
show(gc,'g^{\\mu\\nu}')

# tensor energía-momento (contravariante)
Tc = zeros(4,4)
Tc[0,0] = ρ
Tc[1,1] = g[1,1]*p
Tc[2,2] = g[2,2]*p
Tc[3,3] = g[3,3]*p
show(Tc,'Tc_{\\mu\\nu}')

# métrica contravariante (covariante)
T = Tc.inv()
show(T,'T^{\\mu\\nu}')

# símbolos de Christoffel
Γ = []  # Γ^μ_νσ
for i in range(4):
    Γ.append(zeros(4,4))
    
for μ in range(4):
    line('μ='+str(μ))
    for ν in range(4):
        for σ in range(4):
            aux = 0
            for λ in range(4):
                d1 = dg(λ,ν,σ)
                d2 = dg(λ,σ,ν)
                d3 = dg(ν,σ,λ)
                Γ[μ][ν,σ] += expand(half*gc[μ,λ]*(d1+d2-d3))
            show_Γ(Γ[μ][ν,σ],μ,ν,σ,True)
    #print('########## μ={}, ν={}, σ={}, λ={}'.format(μ,ν,σ,λ))

for μ in range(4):
    show(Γ[μ],'\\Gamma^{'+str(μ)+'}_{\\nu\\sigma}')

# tensor de Riemann
Riemann = [] # R^σ_αβγ
for m in range(4):
    Riemann.append([])
    for n in range(4):
        Riemann[m].append(zeros(4,4))
        
for α in range(4):
    for β in range(4):
        for γ in range(4):
            for σ in range(4):
                t1 = diff(Γ[σ][α,γ], ξ[β], 1)
                t2 = diff(Γ[σ][α,β], ξ[γ], 1)
                Riemann[σ][γ][α,β] += expand(t1-t2)
                for τ in range(4):
                    t3 = Γ[τ][α,γ]*Γ[σ][τ,β]
                    t4 = Γ[τ][α,β]*Γ[σ][τ,γ]
                    Riemann[σ][γ][α,β] += expand(t3-t4)

# tensor de Ricci
Ricci_tensor = zeros(4,4) # R_αγ
for α in range(4):
    for γ in range(4):
        suma = 0
        for σ in range(4):
            suma += expand(Riemann[σ][γ][α,σ])
        Ricci_tensor[α,γ] = suma
        show_R(Ricci_tensor[α,γ],α,γ,True)

show(Ricci_tensor,'R_{\\mu\\nu}')

# escalar de Ricci
Ricci_escalar = 0
for α in range(4):
    for γ in range(4):
        Ricci_escalar += gc[α,γ]*Ricci_tensor[γ,α]
Ricci_escalar = factor(simplify(Ricci_escalar))

show(Ricci_escalar,'R')
#show(Ricci_escalar.subs(k,0),'R(k=0)')

# tensor de Einstein
Einstein_tensor = Ricci_tensor-half*Ricci_escalar*g
show(Einstein_tensor,'G_{\\mu\\nu}')

# ecuaciones de Einstein
LHS = Einstein_tensor+Λ*g
show(LHS,'LHS')

RHS = (8*pi*G/(c**4))*T
show(RHS,'RHS')

EQ=simplify(LHS-RHS)
for i in range (4):
    show(EQ[i,i],'0')

# exportando todo a Latex
fn = os.path.basename(__file__)+'.tex'
with codecs.open(fn,'w', encoding='utf8') as file:
    for s in echo:
        s = s.replace('θ','\\theta')
        s = s.replace('φ','\\phi')
        s = s.replace('μ','\\mu')
        s = s.replace('ν','\\nu')
        s = s.replace('σ','\\sigma')
        s = s.replace('ρ','\\rho')
        s = s.replace('Δ','\\Delta')
        s = s.replace('frac','dfrac')
        s = s.replace('dfrac{d}{d','frac{\partial}{\partial')
        s = s.replace('{\\left(t \\right)}','')
        s = s.replace('\\frac{\\partial}{\\partial t} a','\,\\dot{a}')
        s = s.replace('\dfrac{d^{2}}{d t^{2}} a','\,\\ddot{a}')
        file.write(s+'\n\n')
