using Statistics

"""
plotMalha(plt,coord, inci, index,cor="black";lw = 1.5,ls=:solid)

Se index: \n
  0: Plota a malha sem a numeração dos elementos \n
  1: Plota a malha com os nós e os números dos elementos \n
"""
function plotMalha(plt,coord, inci, index,cor="black";lw = 1.5,ls=:solid)
    nel = size(inci,1)
    nnos = size(coord,1)
    vnos = collect(1:nnos)
    nnel = size(inci,2) - 2
    dimension = size(coord,2)

    for e in 1:nel
        plt = plot!(plt,[coord[inci[e,3:end],1]; coord[inci[e,3],1]],[coord[inci[e,3:end],2]; coord[inci[e,3],2]], lw = lw,ls = ls, leg = false, color = cor)
        if index == 1
            # plt = annotate!([(mean(coord[inci[e,3:end],1]), mean(coord[inci[e,3:end],2]), (string(e),12,:red))])
            plt = scatter!(plt,[mean(coord[inci[e,3:end],1])], [mean(coord[inci[e,3:end],2])], marker = (12, 0.9, :white),series_annotations =[(string(e),12,:red)])
            plt = scatter!(plt,coord[inci[e,3:end],1], coord[inci[e,3:end],2], marker = (10, 0.9, :black), series_annotations = [text(string(inci[e,3+i]),10,:white, :center) for i in 0:nnel-1])

        end
    end
    
    return plt

end

