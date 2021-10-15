function MEFcargaconstante(nel,nnel,L,A0,p0,E)
    nnos = (nnel-1)*nel + 1

    coord = [range(0,5,length = nnos) zeros(nnos)]

    # inci[ Tmat, Tgeo, no1, no2, no3 ]
    inci = zeros(Int64,nel,2+nnel)
    for i in 1:nel
        inci[i,1] = 1
        inci[i,2] = 1
        for j in 1:nnel
            inci[i,2+j] = (nnel-1)*i + 1 - nnel + j
        end
    end

    coorde = zeros((nnel-1)*nel,2)
    for i in 1:nel
        nos = inci[i,3:end]
        x = coord[nos,1]
        y = coord[nos,2]
        for j in 1:nnel-1
            coorde[(nnel-1)*(i-1)+j,1] = (x[j] + x[j+1])/2
            coorde[(nnel-1)*(i-1)+j,2] = (y[j] + y[j+1])/2
        end
    end

    # Tmat - barra de aço
    # Tmat = [E ρ v]
    Tmat = [E 7800 0.3]

    # Tabela de propriedades geométricas
    # O elemento do vetor representa o valor de A_0
    Tgeo = [A0]

    # Plotting
    # plt = plot()
    # plt = plotMalha(plt,coord,inci,1)
    # display(plt)

    ke = zeros(nnel,nnel)
    fe = zeros(nnel,1)
    kg = zeros(nnos,nnos)
    fg = zeros(nnos,1)
    ug = zeros(nnos,1)

    for i in 1:nel
        nos = inci[i,3:end]

        x = coord[nos,1]

        le = x[end]-x[1]

        # Definição da matriz de rigidez
        if nnel == 2
            k0 = [1 -1
            -1 1]
            ke = ((E*A0*(3*L-x[1]-x[2]))/(L*le))*k0

            fe = [p0*le/2, p0*le/2]
        elseif nnel == 3
            k11 = (105*L - 37*x[1] - 36*x[2] + 3*x[3])
            k12 = 4*(-30*L + 11*x[1] + 8*x[2] + x[3])
            k13 = (15*L - 7*x[1] + 4*x[2] - 7*x[3])
            k22 = 16*(15*L - 3*x[1] - 4*x[2] - 3*x[3])
            k23 = 4*(-30*L + x[1] + 8*x[2] + 11*x[3])
            k33 = (105*L + 3*x[1] - 36*x[2] - 37*x[3])

            k0 = [k11 k12 k13
            k12 k22 k23
            k13 k23 k33]
            ke = ((E*A0)/(15*L*le))*k0
            # Defininção do vetor de carga
            fe = p0*[le/6, 2*le/3, le/6]
        end

        # Graus de liberdade dos nós
        dof = nos

        # Soma na matriz global e no vetor de carga global
        kg[dof,dof] = kg[dof,dof]+ke
        fg[dof,1] = fg[dof,1] + fe


    end

    freedofs = 2:nnos-1
    ug[1] = 0
    ug[end] = 0
    ug[freedofs] = kg[freedofs,freedofs]\fg[freedofs]

    sigma = zeros((nnel-1)*nel,1)
    epsilon = zeros((nnel-1)*nel,1)

    for i in 1:nel
        nos = inci[i,3:end]
        x = coord[nos,1]

        for j in 1:nnel-1
            le = x[j+1]-x[j]
            u1 = ug[nos[j]]
            u2 = ug[nos[j+1]]

            ntrecho = (nnel-1)*(i-1)+j

            epsilon[ntrecho] = (u2-u1)/le
            sigma[ntrecho] = E*epsilon[ntrecho]
        end
    end

    return coord, inci, coorde, ug, epsilon, sigma
    # return coord, inci, coorde
end

function MEFcargalinear(nel,nnel,L,A0,p0,E)
    nnos = (nnel-1)*nel + 1

    coord = [range(0,5,length = nnos) zeros(nnos)]

    # inci[ Tmat, Tgeo, no1, no2, no3 ]
    inci = zeros(Int64,nel,2+nnel)
    for i in 1:nel
        inci[i,1] = 1
        inci[i,2] = 1
        for j in 1:nnel
            inci[i,2+j] = (nnel-1)*i + 1 - nnel + j
        end
    end

    coorde = zeros((nnel-1)*nel,2)
    for i in 1:nel
        nos = inci[i,3:end]
        x = coord[nos,1]
        y = coord[nos,2]
        for j in 1:nnel-1
            coorde[(nnel-1)*(i-1)+j,1] = (x[j] + x[j+1])/2
            coorde[(nnel-1)*(i-1)+j,2] = (y[j] + y[j+1])/2
        end
    end

    # Tmat - barra de aço
    # Tmat = [E ρ v]
    Tmat = [E 7800 0.3]

    # Tabela de propriedades geométricas
    # O elemento do vetor representa o valor de A_0
    Tgeo = [A0]

    # Plotting
    # plt = plot()
    # plt = plotMalha(plt,coord,inci,1)
    # display(plt)

    ke = zeros(nnel,nnel)
    fe = zeros(nnel,1)
    kg = zeros(nnos,nnos)
    fg = zeros(nnos,1)
    ug = zeros(nnos,1)

    for i in 1:nel
        nos = inci[i,3:end]

        x = coord[nos,1]

        le = x[end]-x[1]

        # Definição da matriz de rigidez
        if nnel == 2
            k0 = [1 -1
            -1 1]
            ke = ((E*A0*(3*L-x[1]-x[2]))/(L*le))*k0

            fe = (p0*le/(6*L))*[2*x[1] + x[2], x[1] + 2*x[2]]
        elseif nnel == 3
            k11 = (105*L - 37*x[1] - 36*x[2] + 3*x[3])
            k12 = 4*(-30*L + 11*x[1] + 8*x[2] + x[3])
            k13 = (15*L - 7*x[1] + 4*x[2] - 7*x[3])
            k22 = 16*(15*L - 3*x[1] - 4*x[2] - 3*x[3])
            k23 = 4*(-30*L + x[1] + 8*x[2] + 11*x[3])
            k33 = (105*L + 3*x[1] - 36*x[2] - 37*x[3])

            k0 = [k11 k12 k13
            k12 k22 k23
            k13 k23 k33]
            ke = ((E*A0)/(15*L*le))*k0
            # Defininção do vetor de carga
            fe = (p0*le/(30*L))*[4*x[1]+2*x[2]-x[3], 2*x[1]+16*x[2]+2*x[3], -x[1]+2*x[2]+ 4*x[3]]
        end

        # Graus de liberdade dos nós
        dof = nos

        # Soma na matriz global e no vetor de carga global
        kg[dof,dof] = kg[dof,dof]+ke
        fg[dof,1] = fg[dof,1] + fe


    end

    freedofs = 2:nnos-1
    ug[1] = 0
    ug[end] = 0
    ug[freedofs] = kg[freedofs,freedofs]\fg[freedofs]

    sigma = zeros((nnel-1)*nel,1)
    epsilon = zeros((nnel-1)*nel,1)

    for i in 1:nel
        nos = inci[i,3:end]
        x = coord[nos,1]

        for j in 1:nnel-1
            le = x[j+1]-x[j]
            u1 = ug[nos[j]]
            u2 = ug[nos[j+1]]

            ntrecho = (nnel-1)*(i-1)+j

            epsilon[ntrecho] = (u2-u1)/le
            sigma[ntrecho] = E*epsilon[ntrecho]
        end
    end

    return coord, inci, coorde, ug, epsilon, sigma
    # return coord, inci, coorde
end
