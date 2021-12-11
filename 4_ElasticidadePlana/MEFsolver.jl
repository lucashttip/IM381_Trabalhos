function MEFSolver(coord,inci,Tmat,Tgeo,forcing,resticoes)
    
    nel = size(inci,1)
    nnel = size(inci,2)-2
    nnos = size(coord,1)

    ke = zeros(2*nnel,2*nnel)

    kg = zeros(2*nnos,2*nnos)
    fg = zeros(2*nnos,1)
    ug = zeros(2*nnos,1)

    npgauss = 3
    pg, w = gausslegendre(npgauss)

    # Montando a matriz de rigidez e o vetor de forças
    for i in 1:nel
        mat = inci[i,1]
        geo = inci[i,2]
        nos = inci[i,3:end]
        
        E = Tmat[mat,1]
        v = Tmat[mat,3]

        x = coord[nos,1]
        y = coord[nos,2]

        d11 = E/(1-v^2)
        d12 = (v*E)/(1-v^2)
        d13 = 0
        d22 = E/(1-v^2)
        d23 = 0
        d33 = E/(2+2*v)

        D = [d11 d12 d13
            d12 d22 d23
            d13 d23 d33]

        ke = integraelem(pg,w,x,y,D,2*nnel)
        
        # me= A*rho*le*mm

        # Graus de liberdade dos nós
        dof = sort([2*nos; 2*nos.-1])

        # Soma na matriz global e no vetor de carga global
        kg[dof,dof] = kg[dof,dof]+ke

        # fg[dof,1] = fg[dof,1] + fe


    end


    fixeddofs = []
    for i in 1:size(resticoes,1)
        ind = coord[:,resticoes[i,1]].==resticoes[i,2]
        fixednos = findall(ind)
        fixeddofs = [fixeddofs;2*fixednos.-resticoes[i,3]]
        fixeddofs = sort(fixeddofs)
    end

    # fixednodes = [1,22,43,64,85]
    # fixeddofs = sort([2*fixednodes; 2*fixednodes.-1])
    freedofs = setdiff(1:2*nnos,fixeddofs)

    # for i in 1:size(forcing,1)
    #     pos = forcing[i,1:2]
    #     dist = sqrt.((coord[:,1].-pos[1]).^2 + (coord[:,2].-pos[2]).^2)
    #     _,nof = findmin(dist)
    #     fg[2*nof-1] = forcing[3]
    #     fg[2*nof] = forcing[4]
    # end

    for i in 1:size(forcing,1)
        nof = Int(forcing[i,1])
        fg[2*nof-1] = forcing[i,2]
        fg[2*nof] = forcing[i,3]
    end

    # fg[2*nnos] = -40000

    ug[freedofs] = kg[freedofs,freedofs]\fg[freedofs]


    # return ug, kg, fg, freedofs
    return ug
    # return coord, inci, coorde
end

