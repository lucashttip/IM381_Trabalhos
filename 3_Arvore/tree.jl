## IM381 - Tarefa Computacional 3 - MEF - Pórtico
#   Resposta dinâmica de árvore
#   
#   Lucas Agatti Pacheco
#   RA: 235212


using Plots, LinearAlgebra, Revise, PrettyTables
includet("plotfunctions.jl")
includet("auxfunctions.jl")
includet("MEFsolver.jl")

Plots.default(
    guidefontsize=16,
    tickfontsize=14,
    titlefontsize=18,
    legendfontsize=11,
    dpi=200, 
    lw = 1.5,
    size = (800,600))

plotly()

## =================================================
# ==================================================
# ================ Dados de entrada ================
# ==================================================
# ==================================================
salvafigs = true

D = 0.18              # Diametro da base [m]
beta = 1.5          # Parametro comprimento diametro de galho
lambda = 0.5        # Parametro de redução da área
alfad = 20        # Angulo de ramificação [graus]
N = 4               # Número de ramificações
rho = 805           # Massa específica [kg/m³] 
E = 11.3e9            # Módulo de Young [Pa]
v = 0.38            # Coeficiente de Poisson

nelgalho = 1
ngalhos = 2^N - 1
nel = nelgalho*ngalhos
nnel = 2
nnos = (nnel-1)*nel + 1

# Definição dos diametros e dos comprimentos 
Ds = [D*lambda^(i-1) for i in 1:N]
Ls = Ds.^(1/beta)

# Definição das tabelas de material e geométrica
# Tmat
# [Módulo de Young, Massa especifica, coeficiente de poisson]
Tmat = [E rho v]

# Tgeo
# [Área]
Tgeo = [(pi.*Ds.^2)./4 (pi.*(Ds./2).^4)./4]

coord, inci = geramalha(N)

plt = plot()
plt = plotMalha(plt,coord, inci, 0)
plot!(plt,aspect_ratio=:equal)
plot!(plt,xlabel="x [m]", ylabel = "y [m]")
if salvafigs
    savefig(plt,string("./imagens/malha",N,"niveis.png"))
end
display(plt)

plt = plot()
plt = plotMalha(plt,coord, inci, 1)
plot!(plt,aspect_ratio=:equal)
plot!(plt,xlabel="x [m]", ylabel = "y [m]")
if salvafigs
    savefig(plt,string("./imagens/malha",N,"niveisnumerado.png"))
end
display(plt)

## ==============================================
# ===============================================
# ================ Solver do MEF ================
# ===============================================
# ===============================================

freedofs = 4:3*nnos

vals, vecs = MEFSolver(coord, inci, Tmat, Tgeo,freedofs)

## ==================================================
# ===================================================
# ================ Pós-processamento ================
# ===================================================
# ===================================================

# Vibrações livres
nmodos = 10

Freq = sqrt.(real(vals))./(2*pi)

FN = sort(Freq)
I = sortperm(Freq)

fator = 0.1*N

for j in 1:nmodos
    jj = I[j]
    Modo = zeros(3*nnos)
    Modo[freedofs,:] = vecs[:,jj]

    coord2 = coord + fator.*[Modo[1:3:end,1] Modo[2:3:end,1]]
    
    plt1 = plot()
    plt1 = plotMalha(plt1,coord, inci, 0;lw = 1,ls=:dot)
    plt1 = plotMalha(plt1,coord2, inci, 0,"red";lw = 3)
    plt1 = plot!(plt1,aspect_ratio=:equal,title=string("Modo ",j,", Frequência: ", trunc(FN[j],digits=3), " [Hz]"),size = (600,600))
    display(plt1)
    if salvafigs
        savefig(plt1,string("./imagens/",N,"niveismodo",j,".png"))
    end
end

##
header = (["Modo", "Freq. [Hz]"])

dados = [Int.(1:nmodos) Freq[1:nmodos]]

tabela = pretty_table(dados,header = header,  tf = tf_markdown)
# tabela = pretty_table(dados,  tf = tf_markdown)






# Carregamento simples:
# scale = 100

# disp = zeros(size(coord))

# k = 1
# for i in 1:size(disp,1)
#     disp[i,1] = ug[3*i-2]
#     disp[i,2] = ug[3*i-1]
# end

# coord2 = coord + scale.*disp
# plt = plot()
# plt = plotMalha(plt,coord, inci, 0)
# plt = plotMalha(plt,coord2, inci, 0,"red")
# plot(plt,aspect_ratio=:equal)
