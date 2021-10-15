function valoremx(x,coord,inci,ug,N)
    nnel = size(inci,2)-2

    # Encontrando o elemento que o ponto pertence
    e = 0
    for i in 1:size(inci,1)
        nos = inci[i,3:end]
        if x >= coord[nos[1],1] && x <= coord[nos[end],1]
            e = i
            break
        end
    end
    if e == 0
        error("x fora do intervalo")
        return
    end

    nos = inci[e,3:end]
    xs = coord[nos,1]
    us = ug[nos]
    # Caso seja um elemento com 2 nós
    if nnel == 2
        le = xs[end]-xs[1]
        xhat = x - xs[1]
        N1v = N[1](xhat,le)
        N2v = N[2](xhat,le)
        u = N1v*us[1] + N2v*us[2]

    # Caso seja um elemento com 3 nós
    elseif nnel == 3
        le = xs[end]-xs[1]
        xhat = x - xs[1]
        csi = 2*xhat/le - 1
        N1v = N[1](csi)
        N2v = N[2](csi)
        N3v = N[3](csi)
        u = N1v*us[1] + N2v*us[2] + N3v*us[3]
    end
    return u
end

function valoremelem(x,coord,inci,sigma)
    nnel = size(inci,2)-2
    e = 0
    p = 0
    for i in 1:size(inci,1)
        nos = inci[i,3:end]
        if x >= coord[nos[1],1] && x <= coord[nos[end],1]
            e = i
            break
        end
    end
    if e == 0
        error("x fora do intervalo")
        return
    end
    nos = inci[e,3:end]
    for i = 1:nnel-1
        if x >= coord[nos[i],1] && x <= coord[nos[i+1],1]
            p = i
            break
        end
    end
    ntrecho = (nnel-1)*(e-1)+p
    return sigma[ntrecho]
end

function erroL2(coord,inci,valg,valanal,N,x,L,tipo)
    n = length(x)
    valmef = zeros(n)

    if tipo == 1
        # Deslocamento (valores nodais)
        for i in 1:n
            valmef[i] = valoremx(x[i],coord,inci,valg,N)
        end
    elseif tipo == 2
        # Deformação e tensão (valores nos elementos)
        for i in 1:n
            valmef[i] = valoremelem(x[i],coord,inci,valg)
        end
    end

    dL = L/n
    erro = sqrt(sum((valanal-valmef).^2*dL))

    return erro, valmef, x
end