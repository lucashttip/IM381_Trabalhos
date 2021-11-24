function MEFSolver(coord,inci,Tmat,Tgeo,freedofs)
    
    nel = size(inci,1)
    nnel = size(inci,2)-2
    nnos = size(coord,1)

    ke = zeros(3*nnel,3*nnel)
    me = zeros(3*nnel,3*nnel)

    kg = zeros(3*nnos,3*nnos)
    mg = zeros(3*nnos,3*nnos)
    fg = zeros(3*nnos,1)
    ug = zeros(3*nnos,1)

    for i in 1:nel
        mat = inci[i,1]
        geo = inci[i,2]
        nos = inci[i,3:end]

        x1 = coord[nos[1],1]
        x2 = coord[nos[2],1]
        y1 = coord[nos[1],2]
        y2 = coord[nos[2],2]

        E = Tmat[mat,1]
        rho = Tmat[mat,2]
        v = Tmat[mat,3]
        A = Tgeo[geo,1]
        I = Tgeo[geo,2]

        le = sqrt((x2-x1)^2 + (y2 - y1)^2)

        c = (x2 - x1)/le
        s = (y2 - y1)/le

        # Construindo a matriz de rigidez do elemento
        k11 = A*c^2 + s^2*((12*I)/le^2)
        k21 = (A - ((12*I)/le^2))*c*s
        k31 = s*((6*I)/le)
        k41 = - (A*c^2 + ((12*I)/le^2)*s^2)
        k51 = - (A - ((12*I)/le^2))*c*s
        k61 = s*((6*I)/le)

        k22 = A*s^2 + c^2*((12*I)/le^2)
        k32 = -c*((6*I)/le)
        k42 = (-A + ((12*I)/le^2))*c*s
        k52 = - (A*s^2 + ((12*I)/le^2)*c^2)
        k62 = -c*((6*I)/le)

        k33 = 4*I
        k43 = - s*((6*I)/le)
        k53 = c*((6*I)/le)
        k63 = 2*I

        k44 = A*c^2 + s^2*((12*I)/le^2)
        k54 = (A - ((12*I)/le^2))*c*s
        k64 = - s*((6*I)/le)

        k55 = A*s^2 + c^2*((12*I)/le^2)
        k65 = c*((6*I)/le)

        k66 = 4*I


        ke = (E/le)*[k11 k21 k31 k41 k51 k61
                    k21 k22 k32 k42 k52 k62
                    k31 k32 k33 k43 k53 k63
                    k41 k42 k43 k44 k54 k64
                    k51 k52 k53 k54 k55 k65
                    k61 k62 k63 k64 k65 k66]


        m11 = 1/3
        m21 = 0
        m31 = 0
        m41 = 1/6
        m51 = 0
        m61 = 0

        m22 = 13/35
        m32 = 11*le/210
        m42 = 0
        m52 = 9/70
        m62 = -13*le/420
        m33 = le^2/105
        m43 = 0
        m53 = 13*le/420
        m63 = -le^2/140
        m44 = 1/3
        m54 = 0
        m64 = 0
        m55 = 13/35
        m65 = - 11*le/210
        m66 = le^2/105

        R = [c s 0 0 0 0
            -s c 0 0 0 0
            0 0 1 0 0 0
            0 0 0 c s 0
            0 0 0 -s c 0
            0 0 0 0 0 1]

        mm = [m11 m21 m31 m41 m51 m61
            m21 m22 m32 m42 m52 m62
            m31 m32 m33 m43 m53 m63
            m41 m42 m43 m44 m54 m64
            m51 m52 m53 m54 m55 m65
            m61 m62 m63 m64 m65 m66]

        me= A*rho*le*R'*mm*R

        # me= A*rho*le*mm

        # Graus de liberdade dos nós
        dof = [3*nos[1]-2, 3*nos[1]-1, 3*nos[1], 3*nos[2]-2, 3*nos[2]-1, 3*nos[2]]

        # Soma na matriz global e no vetor de carga global
        kg[dof,dof] = kg[dof,dof]+ke
        mg[dof,dof] = mg[dof,dof] + me
        # fg[dof,1] = fg[dof,1] + fe


    end

    # freedofs = 4:3*nnos
    # fg[4] = 10000

    # ug[freedofs] = kg[freedofs,freedofs]\fg[freedofs]


    # nmodos = 4
    vals, vecs = eigen(kg[freedofs,freedofs], mg[freedofs,freedofs])

    # Freq = sqrt.(vals)/(2*pi)

    # FN = sort(Freq)
    # I = sortperm(Freq)


    # for j in 1:nmodos
    #     jj = I[j]
    #     Modo = zeros(2*nnos)
    #     Modo[freedofs] = vecs[:,jj]
    #     plt = plot(coord[:,1],Modo[1:2:end],title=string("Frequencia Natural: ", trunc(FN[j],digits=3), " [Hz]"))
    #     display(plt)
    # end

    # sigma = zeros((nnel-1)*nel,1)
    # epsilon = zeros((nnel-1)*nel,1)

    # for i in 1:nel
    #     nos = inci[i,3:end]
    #     x = coord[nos,1]

    #     for j in 1:nnel-1
    #         le = x[j+1]-x[j]
    #         u1 = ug[nos[j]]
    #         u2 = ug[nos[j+1]]

    #         ntrecho = (nnel-1)*(i-1)+j

    #         epsilon[ntrecho] = (u2-u1)/le
    #         sigma[ntrecho] = E*epsilon[ntrecho]
    #     end
    # end

    return vals, vecs
    # return coord, inci, coorde
end

