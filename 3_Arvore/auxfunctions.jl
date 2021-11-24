function geramalha(N)

    nelgalho = 1
    ngalhos = 2^N - 1
    nel = nelgalho*ngalhos
    nnel = 2
    nnos = (nnel-1)*nel + 1

    coord = zeros(nnos,2)
    inci = zeros(Int,nel,2+nnel)

    coord[1,:] = [0,0]
    coord[2,:] = [0,Ls[1]]
    inci[1,:] = [1, 1, 1, 2]
    ef = 1
    nof = 2

    # Loop para criar a malha
    for i in 1:N-1
        for j in 1:2^(i-1)
            e0 = 2^(i-1) - 1 + j
            no0 = inci[e0,3]
            no1 = inci[e0,4]

            x1 = coord[no1,1]
            y1 = coord[no1,2]

            alfa0 = angulo(coord[no1,:]-coord[no0,:],[1,0])
            alfa1 = alfa0-alfad
            alfa2 = alfa0+alfad

            e1 = ef+1
            e2 = ef+2
            no12 = nof+1
            no22 = nof+2

            y12 = y1 + Ls[i+1]*sind(alfa1)
            y22 = y1 + Ls[i+1]*sind(alfa2)
            x12 = x1 + Ls[i+1]*cosd(alfa1)
            x22 = x1 + Ls[i+1]*cosd(alfa2)

            inci[e1,:] = [1,i+1, no1,no12]
            inci[e2,:] = [1,i+1, no1,no22]

            coord[no12,:] = [x12,y12]
            coord[no22,:] = [x22,y22]

            ef = e2
            nof = no22

        end

    end
    return coord, inci
end

function angulo(a, b)
    return acosd(clamp(dot(a,b)/(norm(a)*norm(b)), -1, 1))
end