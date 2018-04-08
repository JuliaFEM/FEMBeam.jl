# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBeam.jl/blob/master/LICENSE

""" Beam implementation for JuliaFEM. """
module FEMBeam

using FEMBase

"""
    fembeam(l,E,I,A,ro,ne,BCs,qt,qn,F)


l is lenght of the beam
E is Young's modulus
I is Moment of inertia
A is Cross section area
ro is Density
ne is Number of elements
BCs is Boundary conditions, for example with [1,2,3] displacements 1,2 and 3 =0
qt is Vertical Uniform load
qn is Horizontal uniform load
F is point force
"""
function fembeam(l,E,I,A,ro,ne,BCs,qt,qn,F)
    le=l/ne     # Lenght of element
    nn=ne+1     # Number of nodes
    nd=3*ne+3   # Number of DOFs
    Gp=[-1/3*sqrt(5+2*sqrt(10/7)) -1/3*sqrt(5-2*sqrt(10/7)) 0 1/3*sqrt(5-2*sqrt(10/7)) 1/3*sqrt(5+2*sqrt(10/7))]
    w=[(322-13*sqrt(70))/900 (322+13*sqrt(70))/900 128/225 (322+13*sqrt(70))/900 (322-13*sqrt(70))/900]
    # Integration of the truss elments stiffness matrix
    detJ=2/le
    function tkint(w,xi)
        dN1=-1/2
        dN2=1/2
        dN=[dN1 dN2]
        w*dN'*dN
    end
    tk=0
    for i = 1:size(Gp,2)
        tk +=tkint(w[i],Gp[i])
    end
    tk=detJ*A*E*tk
    # Integration of the 4 DOF beam elments stiffness matrix
    detJ=(2/le)^3
    function bkint(w,xi)
        d2N1=(-1.0)*(1 - xi) + 0.5*(2 + xi)
        d2N2=(-1/2)*le*(1 - xi) + (1/4)*le*(1 + xi)
        d2N3=(-1.0)*(1 + xi) + 0.5*(2 - xi)
        d2N4=(1/4)*le*(-1 + xi) + (1/2)*le*(1 + xi)
        d2N=[d2N1 d2N2 d2N3 d2N4]
        w*d2N'*d2N
    end
    bk=0
    for i = 1:size(Gp,2)
        bk +=bkint(w[i],Gp[i])
    end
    bk=E*I*detJ*bk
    # Assembly of 6 DOF truss-beam stiffness matrix
    k=zeros(6,6)
    k[1,1]=tk[1,1]; k[1,4]=tk[1,2]; k[4,1]=tk[2,1]; k[4,4]=tk[2,2]
    k[2,2]=bk[1,1];k[2,3]=bk[1,2];k[2,5]=bk[1,3];k[2,6]=bk[1,4]
    k[3,2]=bk[2,1];k[3,3]=bk[2,2];k[3,5]=bk[2,3];k[3,6]=bk[2,4]
    k[5,2]=bk[3,1];k[5,3]=bk[3,2];k[5,5]=bk[3,3];k[5,6]=bk[3,4]
    k[6,2]=bk[4,1];k[6,3]=bk[4,2];k[6,5]=bk[4,3];k[6,6]=bk[4,4]
    # Global stiffness matrix K
    K=zeros(nd,nd)
    for a in 0:ne-1
    K_temp = zeros((nd,nd))
    for i in 1:6
        for j in 1:6
            K_temp[a*3+i,a*3+j]=k[j,i]
        end
    end
    K +=K_temp
    end
    # Integration of the truss mass matrix
    detJ=le/2
    function tmint(w,xi)
        N1s=1+(-1/2)*(1+xi)
        N2s=(1/2)*(1 + xi)
        Ns=[N1s N2s]
        w*Ns'*Ns
    end
    tm=0
    for i = 1:size(Gp,2)
        tm +=tmint(w[i],Gp[i])
    end
    tm=ro*A*tm*detJ
    # Integration of the beam mass matrix
    detJ=(le/2)
    function bmint(w,xi)
        N1=1/4*(1-xi)^2*(2+xi)
        N2=le/8*(1-xi)^2*(xi+1)
        N3=1/4*(1+xi)^2*(2-xi)
        N4=le/8*(1+xi)^2*(xi-1)
        N=[N1 N2 N3 N4]
        w*N'*N
    end
    bm=0
    for i = 1:size(Gp,2)
        bm +=bmint(w[i],Gp[i])
    end
    bm=ro*A*bm*detJ
    # Assembly of the 6 DOF truss-beam mass matrix
    m=zeros(6,6)
    m[1,1]=tm[1,1];m[1,4]=tm[1,2];m[4,1]=tm[2,1];m[4,4]=tm[2,2]
    m[2,2]=bm[1,1];m[2,3]=bm[1,2];m[2,5]=bm[1,3];m[2,6]=bm[1,4]
    m[3,2]=bm[2,1];m[3,3]=bm[2,2];m[3,5]=bm[2,3];m[3,6]=bm[2,4]
    m[5,2]=bm[3,1];m[5,3]=bm[3,2];m[5,5]=bm[3,3];m[5,6]=bm[3,4]
    m[6,2]=bm[4,1];m[6,3]=bm[4,2];m[6,5]=bm[4,3];m[6,6]=bm[4,4]
    # Global mass matrix M
    M=zeros(nd,nd)
    for a in 0:ne-1
    M_temp = zeros((nd,nd))
    for i in 1:6
        for j in 1:6
            M_temp[a*3+i,a*3+j]=m[j,i]
        end
    end
    M +=M_temp
    end
    # Integration of the equivalent forces vector
    # For truss element
    detJ=le/2
    function tfqeint(w,xi)
        N1s=1+(-1/2)*(1+xi)
        N2s=(1/2)*(1 + xi)
        Ns=[N1s N2s]
        w*Ns'
    end
    tfqe=zeros(2,1)
    for i = 1:size(Gp,2)
        tfqe +=tfqeint(w[i],Gp[i])
    end
    tfqe=qt*tfqe*detJ
    # For 4 DOF beam element
    detJ=le/2
    function bfqeint(w,xi)
        N1=1/4*(1-xi)^2*(2+xi)
        N2=le/8*(1-xi)^2*(xi+1)
        N3=1/4*(1+xi)^2*(2-xi)
        N4=le/8*(1+xi)^2*(xi-1)
        N=[N1 N2 N3 N4]
        w*N'
    end
    bfqe=zeros(4,1)
    for i = 1:size(Gp,2)
        bfqe +=bfqeint(w[i],Gp[i])
    end
    bfqe=qn*bfqe*detJ
    # Assembly of the 6 DOF beam element equivalent forces vector
    fql=zeros(6,1)
    fql[1,1],fql[4,1]=tfqe[1,1],tfqe[2,1]
    fql[2,1],fql[3,1],fql[5,1],fql[6,1]=bfqe[1,1],bfqe[2,1],bfqe[3,1],bfqe[4,1]
    # Assembly of the global equivalent forces vector with uniformly distributed
    # load over the whole beam
    fq=zeros(nd)
    fq[1:3]    =fql[1:3]
    fq[nd-2:nd]=fql[4:6]
    for i in 1:ne-1
        fq[3*i+1]=fql[1]+fql[4]
        fq[3*i+2]=fql[2]+fql[5]
        fq[3*i+3]=fql[3]+fql[6]
    end
    # Point forces vector
    f=zeros(nd)
    f[nd-1,1]=F
    # Adding equivalent forces vector to point forces vector
    f +=fq
    # Lagrange multiplier method
    A0=zeros(size(BCs,1),nd)
    b=zeros(size(BCs,1))
    for i = BCs
        A0[i,i]=1
    end
    f=[f; b]
    K=[K    A0';
        A0  zeros(size(BCs,1),size(BCs,1))]
    # Solving static problem
    U=K\f
    # Collecting Lagrange multipliers from U vector
    la=U[nd+1:end]
    # Removing Lagrange multipliers from U vector
    U_temp=U
    K_temp=K
    U=0;K=0
    U=U_temp[1:nd]
    K=K_temp[1:nd,1:nd]
    return U,K,M,la
end

end
