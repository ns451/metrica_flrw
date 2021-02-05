# Relatividad General | 2020-02-04 | t.me/ns451

import os
import codecs
from sympy import *

init_printing(pretty_print=True, wrap_line=False)

half = Rational(1/2)

# variables
k = symbols('k')

# coordenadas
t,r,θ,φ, = symbols('t r θ φ ')
ξ = [t,r,θ,φ]

# funciones
a = Function('a')(t)

echo = [] # para exportar

# muestra símbolo de Christoffel
def show_Γ(e,μ,ν,σ,hide_zeros):
    #s = expand(e)
    s=e
    if str(s)!='0' and hide_zeros:
        print('\n{} Γ^[{}]_[{}][{}] =\n'.format('■'*5,μ,ν,σ))
        #pprint(s)
        a = '\n\Gamma^{}_{}={}'
        sup = '{'+str(μ)+'}'
        sub = '{'+str(ν)+str(σ)+'}'
        echo.append(a.format(sup,sub,latex(s)))

# muestra componenente de tensor de Ricci
def show_R(e,α,γ,hide_zeros):
    #s = expand(e)
    s=e
    if str(s)!='0' and hide_zeros:
        print('\n{} R_[{}][{}] =\n'.format('■'*5,α,γ))
        #pprint(s)
        a = '\nR_{}={}'
        sub = '{'+str(α)+str(γ)+'}'
        echo.append(a.format(sub,latex(s)))

# muestra tensor rank 4
def show(e):
    print('\n{} {} ='.format('■'*15,e))
    s = expand(eval(e))
    #pprint(s)
    echo.append('\n'+e+'='+latex(s))

# derivada covariante
def dg(α,β,γ):
    return diff(g[α,β], ξ[γ], 1)

# línea de separación
def line(e):
    print('■'*30+' '+e)

# métrica covariante
g = zeros(4,4)
g[0,0] = 1
g[1,1] = -a**2/(1-k*r**2)
g[2,2] = -a**2*r**2
g[3,3] = -a**2*r**2*(sin(θ))**2
show('g')

# metrica contravariante
gc = g.inv()
show('gc')

# símbolos de Christoffel
Γ = [zeros(4,4)]*4 # Γ^μ_νσ
for μ in range(4):
    line('μ='+str(μ))
    for ν in range(4):
        for σ in range(4):
            for λ in range(4):
                d1 = dg(λ,ν,σ)
                d2 = dg(λ,σ,ν)
                d3 = dg(ν,σ,λ)
                Γ[μ][ν,σ] += expand(half*gc[μ,λ]*(d1+d2-d3))
            show_Γ(Γ[μ][ν,σ],μ,ν,σ,True)

# tensor de Riemann
Riemann = [[zeros(4,4)]*4]*4 # R^σ_αβγ
for α in range(4):
    for β in range(4):
        for γ in range(4):
            for σ in range(4):
                t1 = diff(Γ[σ][α,γ], ξ[β], 1)
                t2 = diff(Γ[σ][α,β], ξ[γ], 1)
                for τ in range(4):
                    t3 = Γ[τ][α,γ]*Γ[σ][τ,β]
                    t4 = Γ[τ][α,β]*Γ[σ][τ,γ]
                    Riemann[σ][γ][α,β] += expand(t1-t2+t3-t4)

# tensor de Ricci
Ricci_tensor = zeros(4,4) # R_αγ
for α in range(4):
    for γ in range(4):
        Ricci_tensor[α,γ] += Riemann[1][γ][α,1]
        show_R(Ricci_tensor[α,γ],α,γ,True)

# escalar de Ricci
Ricci_escalar = 0
for α in range(4):
    Ricci_escalar += Ricci_tensor[α,α]
show('Ricci_escalar')

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
