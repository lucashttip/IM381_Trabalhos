using Statistics

"""
Se index: \n
  1: Plota a malha com os nós e os números dos elementos \n
  2: Plota a malha sem a numeração dos elementos
"""
function plotMalha(plt,coord, inci, index)
    nel = size(inci,1)
    nnos = size(coord,1)
    vnos = collect(1:nnos)
    nnel = size(inci,2) - 2
    dimension = size(coord,2)

    for e in 1:nel
        plt = plot!(plt,[coord[inci[e,3:end],1]; coord[inci[e,3],1]],[coord[inci[e,3:end],2]; coord[inci[e,3],2]], lw = 2, leg = false, color = "black")
        if index == 1
            # plt = annotate!([(mean(coord[inci[e,3:end],1]), mean(coord[inci[e,3:end],2]), (string(e),12,:red))])
            plt = scatter!(plt,[mean(coord[inci[e,3:end],1])], [mean(coord[inci[e,3:end],2])], marker = (12, 0.9, :white),series_annotations =[(string(e),12,:red)])
            plt = scatter!(plt,coord[inci[e,3:end],1], coord[inci[e,3:end],2], marker = (10, 0.9, :black), series_annotations = [text(string(inci[e,3+i]),10,:white, :center) for i in 0:nnel-1])

        end
    end
    
    return plt

end

"""
Faz plots bonitos comparando a solução analítica e as soluções numéricas
"""
function plotbonito(x,xa,y,ya,nel,titulo,nomex,nomey,escala,markers,posleg,linhas)
    plt = plot(xa,ya.*escala,title=titulo,label="analítico",xlabel = nomex, ylabel = nomey,legend=posleg,lw = 2)
    for i in 1:length(y)
        plt = plot!(plt,x[i][:,1],y[i].*escala,marker = markers[i],line = linhas[i],label=string(nel[i]," nós"))
    end
    return plt
end

# """
# Se index: \n
#   1: Plota a malha com os nós e os números dos elementos \n
#   2: Plota a malha sem anumeração dos elementos
# """
# function plotResultado(plt,coord, inci, z, index, cor = true)
#     nel = size(inci,1)
#     nnos = size(coord,1)
#     vnos = collect(1:nnos)
#     nnel = size(inci,2) - 2
#     dimension = size(coord,2)

#     for e in 1:size(inci,1)
#         plt = plot!(plt,coord[inci[e,3:end],1],coord[inci[e,3:end],2], leg = false, color = "black")
#         if index == 1
#             plt = annotate!([(mean(coord[inci[e,3:end],1]), mean(coord[inci[e,3:end],2]), (string(e),10,:red))])
#         end
#     end
    
#     return plt

# end

