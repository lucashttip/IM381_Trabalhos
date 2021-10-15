# Elemento de barra com 2 ou 3 nós
#     E(A u')' + p0 = 0     in (0,5)
#     u(0) = 0
#     u(L) = 0

# Pacotes necessários
using Revise
using Plots
using LinearAlgebra
using LaTeXStrings 
using Printf

# Definições para plots bonitos
Plots.default(guidefontsize=16,tickfontsize=14,titlefontsize=18,legendfontsize=11,dpi=200, lw = 1.5,size = (800,600))

# Arquivos necessários
includet("plotfunctions.jl")
includet("functionsmef.jl")
includet("MEFsolver.jl")

# Configurações de saída
salvafigs = false    # Salvar os plots em arquivos png?

# Dados de entrada
L = 5               # Comprimento da barra [m]
A0 = 0.01*0.01      # Área da seção [m^2]
nel = [4, 8, 16, 32]        # Número de elementos a cada iteração (ter que ser vetor)
nnel = [2,2,2,2]        # Número de nós por elemento a cada iteração (tem que ser vetor)
p0 = 1000           # carga [N/m^2]
tipocarga = :const  # Se a carga é constante ou linear (:const, :linear)
E = 2.07e11         # Módulo de elasticidade [Pa]

niter = length(nel)

# Definindo a função analítica e intervalo de avaliação 
N = 200             # Número de pontos de avaliação da função analítica
xa = range(0,L,length = N)
if tipocarga ==:const
    u_a(x) = (L^2*p0*log(3 - 2*x/L))/(2*A0*E*log(3)) + (L*p0*(x - L)/(2*A0*E))
    epsilon_a(x) = (-p0*x - ((L*p0)/log(3)) + ((3*L*p0)/2))/(E*A0*(3-2*x/L))
    sigma_a(x) = epsilon_a(x)*E
elseif tipocarga ==:linear
    u_a(x) = (L^2*p0*log(3-2*x/L))/(2*A0*E*log(3)) + (p0*(3*L*x + x^2 - 4*L^2))/(8*A0*E)
    epsilon_a(x) = ((9*p0*L/8) - ((L*p0)/log(3)) - ((p0*x^2)/(2*L)))/(E*A0*(3-2*x/L))
    sigma_a(x) = epsilon_a(x)*E
end

# Cálculo dos vetores das soluções numéricas
des_a = u_a.(xa)
tens_a = sigma_a.(xa)
def_a = epsilon_a.(xa)

# Definido as funções de forma (usadas no cálculo dos erros)
# Para elemento com 2 nós:
N1(x,le) = 1-x/le
N2(x,le) = x/le

# Para elemento com 3 nós: (funções diferentes)
N1(xi) = xi^2/2 - xi/2
N2(xi) = 1 - xi^2
N3(xi) = xi^2/2 + xi/2

Ni = [N1, N2, N3]

# Inicialização das variáveis de saída do solver
ug = Vector{Any}(undef,niter)       # Vetor dos deslocamentos
coord = Vector{Any}(undef,niter)    # Vetor das coordenadas dos nós
inci = Vector{Any}(undef,niter)     # Vetor das incidências
coorde = Vector{Any}(undef,niter)   # Vetor dos pontos entre dois nós do mesmo elemento (avaliação de tensão e deformação)
epsilon = Vector{Any}(undef,niter)  # Vetor das deformações
sigma = Vector{Any}(undef,niter)    # Vetor das tensões

# Código da solução
for i in 1:niter
    if tipocarga == :const
        coord[i], inci[i], coorde[i], ug[i], epsilon[i], sigma[i] = MEFcargaconstante(nel[i],nnel[i],L,A0,p0,E)
    elseif tipocarga ==:linear
        coord[i], inci[i], coorde[i], ug[i], epsilon[i], sigma[i] = MEFcargalinear(nel[i],nnel[i],L,A0,p0,E)
    end
end

# Definição dos marcadores e das linhas dos plots 
# cada tupla no vetor é a definição de uma curva das soluções numéricas 
# Atualmente são suportadas apenas até 4 curvas numéricas diferentes.
markers = [(:circle, 5,:red), (:star5, 5,:green), (:cross, 5,:orange), (:xcross, 5,:black)]
linhas = [(1, :dash, 0.7,:red),(1, :dash, 0.7,:green),(1, :dash, 0.7,:orange),(1, :dash, 0.7,:black)]

# Plot dos resultados
plt1 = plotbonito(coord,xa,ug,des_a,nel,nnel,"Deslocamento",L"x \ [m]",L"u \ [mm]",1e3,markers,:topleft,linhas)
plt2 = plotbonito(coorde,xa,sigma,tens_a,nel,nnel,"Tensão",L"x \ [m]",L"\sigma \ [MPa]",1e-6,markers,:topright,linhas)
plt3 = plotbonito(coorde,xa,epsilon,def_a,nel,nnel,"Deformação",L"x \ [m]",L"\varepsilon \ [\%]",1e2,markers,:topright,linhas)

display(plt1)
display(plt2)
display(plt3)

# Alocação das variáveis dos erros
eru = Vector{Any}(undef,niter)
umef = Vector{Any}(undef,niter)
xu = Vector{Any}(undef,niter)
ert = Vector{Any}(undef,niter)
tmef = Vector{Any}(undef,niter)
xt = Vector{Any}(undef,niter)
erd = Vector{Any}(undef,niter)
dmef = Vector{Any}(undef,niter)
xd = Vector{Any}(undef,niter)

# Cálculo dos erros em comparação com a solução numérica
for i in 1:niter
    eru[i], umef[i], xu[i] = erroL2(coord[i],inci[i],ug[i],des_a,Ni,xa,L,1)
    ert[i], tmef[i], xt[i] = erroL2(coord[i],inci[i],sigma[i],tens_a,Ni,xa,L,2)
    erd[i], dmef[i], xd[i] = erroL2(coord[i],inci[i],epsilon[i],def_a,Ni,xa,L,2)
end

# linhas = [(1, :solid, 0.7,:red),(1, :solid, 0.7,:green),(1, :dash, 0.7,:orange),(1, :dash, 0.7,:black)]
# plt4 = plotbonito(xt,xa,umef[1:2],des_a,nnel,"Deslocamento",L"x \ [m]",L"u \ [mm]",1e3,[(:none),(:none),(:none),(:none)],:topright,linhas)
# plt4 = plotbonito(xt,xa,tmef[1:2],tens_a,nnel,"Tensão",L"x \ [m]",L"\sigma \ [MPa]",1e-6,[(:none),(:none),(:none),(:none)],:topright,linhas)


# Em caso as figuras sejam salvas
if salvafigs == true
    savefig(plt1,string("./imagens/deslocamento.png"))
    savefig(plt2,string("./imagens/tensao.png"))
    savefig(plt3,string("./imagens/deformacao.png"))
    # savefig(plt4,string("../imagens/tensaopx2.png"))
end
