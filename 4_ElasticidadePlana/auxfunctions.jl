function geramalha(nel,nely,nelx,nnosx)
    inci = zeros(Int64,nel,6)
    inci[:,1:2] .= 1

    e = 1;
    c = 1;
    for j in 1:nely
        for i = 1:nelx
            inci[e,3:6] = [c, c+1, c+nnosx+1, c+nnosx];
            e = e+1
            c = c+1;
        end
    c = c+1;
    end 
    return inci
end


N(csi,eta) = 1/4*[1 - csi - eta + csi*eta
                1 + csi - eta - csi*eta
                1 + csi + eta + csi*eta
                1 - csi + eta - csi*eta]

dNdcsi(csi,eta) = 1/4*[-1 + eta
                    1 - eta
                    1 + eta
                    -1 - eta]

dNdeta(csi,eta) = 1/4*[-1 + csi
                    -1 - csi
                    1 + csi
                    1 - csi]

dNdcsieta(csi,eta) = [dNdcsi(csi,eta)'; dNdeta(csi,eta)']

J(csi,eta,x,y) = dNdcsieta(csi,eta)*[x y]

function Bcsieta(csi,eta,jac)
    

    dN = inv(jac)*dNdcsieta(csi,eta)

    dNdx = dN[1,:]
    dNdy = dN[2,:]


    B = [dNdx[1] 0 dNdx[2] 0 dNdx[3] 0 dNdx[4] 0
        0 dNdy[1] 0 dNdy[2] 0 dNdy[3] 0 dNdy[4]
        dNdy[1] dNdx[1] dNdy[2] dNdx[2] dNdy[3] dNdx[3] dNdy[4] dNdx[4]]

    return B
end

function calculaponto(csi,eta,jac,D)
    B = Bcsieta(csi,eta,jac)
    m = B'*D*B*det(jac)

    return m
end

function integraelem(pontos,pesos,x,y,D,n)
    ke = zeros(n,n)
    npgauss = size(pontos,1)

    for i in 1:npgauss
        for j in 1:npgauss
            csi = pontos[i]
            eta = pontos[j]
            jac = J(csi,eta,x,y) 
            ke = ke + pesos[i]*pesos[j]*calculaponto(pontos[i],pontos[j],jac,D)
        end
    end
    return ke
end

function pontoselem(coord,inci)
    nel = size(inci,1)
    
    pontos = zeros(nel,2)


    for i in 1:nel
        nos = inci[i,3:end]
        x = mean(coord[nos,1])
        y = mean(coord[nos,2])
        pontos[i,1] = x
        pontos[i,2] = y
    end
    return pontos

end

function calculadeformacao(disp,coord,inci)
    nel = size(inci,1)

    def = zeros(nel,3)

    for i in 1:nel
        nos = inci[i,3:end]
        dxx = (disp[nos[2],1] - disp[nos[1],1])/(coord[nos[2],1] - coord[nos[1],1])
        dyy = (disp[nos[4],2] - disp[nos[1],2])/(coord[nos[4],2] - coord[nos[1],2])
        dxy = ((disp[nos[2],2] - disp[nos[1],2])/(coord[nos[2],1] - coord[nos[1],1])) + ((disp[nos[4],1] - disp[nos[1],1])/(coord[nos[4],2] - coord[nos[1],2]))
        
        def[i,1] = dxx
        def[i,2] = dyy
        def[i,3] = dxy
    end

    return def
end

function calculatensoes(E,v,def)
    nel = size(def,1)

    tens = zeros(nel,3)

    for i in 1:nel
        defelem = def[i,:]
        tenselem = (E/(1-v^2))*[1 v 0
                                v 1 0
                                0 0 (1-v)/2]*defelem
        tens[i,:] = tenselem
    end

    return tens
end

function vonmises(tens)
    tvm = zeros(size(tens,1))
    for i in 1:size(tens,1)
        tvm[i] = sqrt((tens[i,1]-tens[i,2])^2 + 3*tens[i,3]^2)
    end
    return tvm
end