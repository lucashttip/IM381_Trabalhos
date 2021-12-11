## Aproximação escoamento potencial
# 
#   k * nabla^2 T = Qo = 20
#   T1 = 0 no lado L1 (abaixo)
#   T3 = 0 no lado L3 (acima)
#   dT2dx = 0 no lado L2 > Condição Natural
#   dT4dx = 10 no lado L4
# 
#   Elemento quadrilateral de 3 nós
# 

using Plots, LinearAlgebra, Revise, FastGaussQuadrature, Statistics
import Makie, CairoMakie

includet("plotfunctions.jl")
includet("MEFsolver.jl")
includet("auxfunctions.jl")
BLAS.set_num_threads(1)

Plots.default(
    guidefontsize=16,
    tickfontsize=14,
    titlefontsize=18,
    legendfontsize=11, 
    dpi=200, 
    lw = 1.5,
    size = (800,600))

## =================================================
# ==================================================
# ================ Dados de entrada ================
# ==================================================
# ==================================================
salvafigs = false

Lx = 10
Ly = 2
nelx = 100
nely = 20

E = 2.07e11
rho = 7800
v = 0.3

F = - 40000

nnosx = nelx + 1
nnosy = nely + 1
nel = Int(nelx*nely)
nnos = Int((nelx + 1)* (nely + 1))

coord = [repeat(range(0,Lx,length=nnosx),outer=(nnosy)) repeat(range(0,Ly,length=nnosy),inner=(nnosx))]


#  inci[ Tmat , Tgeo , no1 , no 2, no3 ]

inci = geramalha(nel,nely,nelx,nnosx)

# Tmat
Tmat = [E rho v]

# Tgeo
Tgeo = [0]

plt = plot()
plt = plotMalha(plt,coord, inci, 1)
plot!(plt,aspect_ratio=:equal)



# ==============================================
# ===============================================
# ================ Solver do MEF ================
# ===============================================
# ===============================================

# forcing = [posição x do nó, posição y do nó, força em x, força em y]
# forcing = [10 2 0 -40000]


# forcing = [no, força em x, força em y]
xi = 8
xf = 10
yf = 2
nosf = (coord[:,1].>=xi .&& coord[:,1].<=xf) .&& coord[:,2].==yf
nosf = findall(nosf)

fno = (xf-xi)*F/size(nosf,1)

forcing = zeros(size(nosf,1),3)
for i in 1:size(nosf,1)
    forcing[i,:] = [nosf[i], 0, fno] 
end



# restricoes = [lado (1 pra x, 2 pra y), poslado, coordfix(0 pra y, 1 pra x)]
restricoes = [1 0 0
            1 0 1]

u = MEFSolver(coord, inci, Tmat, Tgeo, forcing,restricoes)
# ug, kg, fg, freedofs = MEFSolver(coord, inci, Tmat, Tgeo, forcing,restricoes)

## ==================================================
# ===================================================
# ================ Pós-processamento ================
# ===================================================
# ===================================================


disp = [u[1:2:end] u[:2:2:end]]

coord2 = coord + 1e5*disp

plt1 = plot()
plt1 = plotMalha(plt1,coord, inci, 0;lw = 2,ls=:solid)
plot!(plt1,aspect_ratio=:equal,title="Malha original",xlabel = "x [m]", ylabel = "y [m]")
savefig(plt1,"imagens/MalhaOriginal.png")

plt1 = plot()
plt1 = plotMalha(plt1,coord2, inci, 0,"red";lw = 2)
plot!(plt1,aspect_ratio=:equal,title="Malha deformada com ampliação de 100.000x",xlabel = "x [m]", ylabel = "y [m]")
savefig(plt1,"imagens/MalhaDeformada.png") 


plt1 = plot()
plt1 = plotMalha(plt1,coord, inci, 0;lw = 2,ls=:solid)
plt1 = plotMalha(plt1,coord2, inci, 0,"red";lw = 2)
plot!(plt1,aspect_ratio=:equal,title="Comparação das malhas com ampliação de 100.000x",xlabel = "x [m]", ylabel = "y [m]")
savefig(plt1,"imagens/ComparaMalha.png")
## Cálculo das tensões

pontos = pontoselem(coord,inci)
def = calculadeformacao(disp,coord,inci)
tens = calculatensoes(E,v,def)
tvm = vonmises(tens)

fig= plotMesh(pontos[:,1],pontos[:,2],tens[:,1]./1e3,"Tensão normal em x","σxx [kPa]")
Makie.save("imagens/tens/TensX.png",fig)

fig= plotMesh(pontos[:,1],pontos[:,2],tens[:,2]./1e3,"Tensão normal em y","σyy [kPa]")
Makie.save("imagens/tens/TensY.png",fig)

fig= plotMesh(pontos[:,1],pontos[:,2],tens[:,3]./1e3,"Tensão de cisalhamento","σxy [kPa]")
Makie.save("imagens/tens/TensXY.png",fig)

fig= plotMesh(pontos[:,1],pontos[:,2],tvm./1e3,"Von Mises","σvm [kPa]")
Makie.save("imagens/tens/VonMises.png",fig)

fig= plotMesh(pontos[:,1],pontos[:,2],def[:,1]./1e3,"Deformação normal em x","σxx [kPa]")
Makie.save("imagens/def/defX.png",fig)

fig= plotMesh(pontos[:,1],pontos[:,2],def[:,2]./1e3,"Deformação normal em y","σyy [kPa]")
Makie.save("imagens/def/defY.png",fig)

fig= plotMesh(pontos[:,1],pontos[:,2],def[:,3]./1e3,"Deformação de cisalhamento","σxy [kPa]")
Makie.save("imagens/def/defXY.png",fig)


fig= plotMesh(coord[:,1],coord[:,2],disp[:,1]./1e3,"Deslocamentos em x","u [mm]")
Makie.save("imagens/deslocamentosx.png",fig)

fig= plotMesh(coord[:,1],coord[:,2],disp[:,2]./1e3,"Deslocamentos","v [mm]")
Makie.save("imagens/deslocamentosy.png",fig)

fig= plotMesh(coord[:,1],coord[:,2],sqrt.(disp[:,1].^2 + disp[:,2].^2)./1e3,"Deslocamentos","Deslocamento total [mm]")
Makie.save("imagens/deslocamentost.png",fig)


