# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBeam.jl/blob/master/LICENSE

""" Beam implementation for JuliaFEM. """
module FEMBeam

using FEMBase

"""
    fembeam(l,E,I,A,ro,ne,BCs,qt,qn,F)

Function integrates
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
    Gp=[-sqrt(3/7+2/7*sqrt(6/5)) -sqrt(3/7-2/7*sqrt(6/5)) sqrt(3/7-2/7*sqrt(6/5)) sqrt(3/7+2/7*sqrt(6/5))]
    w=[(18-sqrt(30))/36 (18+sqrt(30))/36 (18+sqrt(30))/36 (18-sqrt(30))/36]
    # Integration of the truss elments stiffness matrix
    detJ=2/le
    function tkint(w,xi)
    w*
    [0.25 -0.25; -0.25 0.25]
    end
    tk=0
    for i = 1:size(Gp,2)
        tk +=tkint(w[i],Gp[i])
    end
    tk=detJ*A*E*tk
    # Integration of the 4 DOF beam elments stiffness matrix
    detJ=(2/le)^3
    function bkint(w,xi)
    w*
    [(-3/2 + (3/2)*(1 + xi))^2 (-3/2 + (3/2)*(1 + xi))*(le*(-1 + (1/2)*(1 + xi)) + (1/4)*le*(1 + xi)) (-3/2 + (3/2)*(1 + xi))*(-(1 + xi) + (1/2)*(3 - (1 + xi))) (-3/2 + (3/2)*(1 + xi))*((1/2)*le*(-1 + (1/2)*(1 + xi)) + (1/2)*le*(1 + xi));
    (-3/2 + (3/2)*(1 + xi))*(le*(-1 + (1/2)*(1 + xi)) + (1/4)*le*(1 + xi)) (le*(-1 + (1/2)*(1 + xi)) + (1/4)*le*(1 + xi))^2 (-(1 + xi) + (1/2)*(3 - (1 + xi)))*(le*(-1 + (1/2)*(1 + xi)) + (1/4)*le*(1 + xi)) ((1/2)*le*(-1 + (1/2)*(1 + xi)) + (1/2)*le*(1 + xi))*(le*(-1 + (1/2)*(1 + xi)) + (1/4)*le*(1 + xi));
    (-3/2 + (3/2)*(1 + xi))*(-(1 + xi) + (1/2)*(3 - (1 + xi))) (-(1 + xi) + (1/2)*(3 - (1 + xi)))*(le*(-1 + (1/2)*(1 + xi)) + (1/4)*le*(1 + xi)) (-(1 + xi) + (1/2)*(3 - (1 + xi)))^2 (-(1 + xi) + (1/2)*(3 - (1 + xi)))*((1/2)*le*(-1 + (1/2)*(1 + xi)) + (1/2)*le*(1 + xi));
    (-3/2 + (3/2)*(1 + xi))*((1/2)*le*(-1 + (1/2)*(1 + xi)) + (1/2)*le*(1 + xi)) ((1/2)*le*(-1 + (1/2)*(1 + xi)) + (1/2)*le*(1 + xi))*(le*(-1 + (1/2)*(1 + xi)) + (1/4)*le*(1 + xi)) (-(1 + xi) + (1/2)*(3 - (1 + xi)))*((1/2)*le*(-1 + (1/2)*(1 + xi)) + (1/2)*le*(1 + xi)) ((1/2)*le*(-1 + (1/2)*(1 + xi)) + (1/2)*le*(1 + xi))^2]
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
    # Massmatrix

    ## working on this

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
    # equivalent forces vector ## I didn't get the integration working with this
    fql=[qt*le/2;      # x
         qn*le/2;      # y
         qn*le^2/12;   # M
         qt*le/2;      # x
         qn*le/2;      # y
        -qn*le^2/12]   # M
    fq=zeros(nd)
    fq[1:3]    =fql[1:3]
    fq[nd-2:nd]=fql[4:6]
    for i in 1:ne-1
        fq[3*i+1]=fql[1]+fql[4]
        fq[3*i+2]=fql[2]+fql[5]
        fq[3*i+3]=fql[3]+fql[6]
    end
    f=zeros(nd)
    f[nd-1,1]=F
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
