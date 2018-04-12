# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBeam.jl/blob/master/LICENSE

""" Beam implementation for JuliaFEM. """
module FEMBeam

using FEMBase

"""
Function integrates stiffness matrix, mass matrix and forces vector for
6 DOF Euler-Bernoulli beam element.

    fembeam(X1,X2,E,I,A,ro,ne,qt,qn,f)

X1 = beams left node coordinates
X2 = beams right node coordinates
E = Young's modulus
I = Moment of inertia
A = Cross section area
ro = Density
qt = Tangential uniformly distributed load
qn = Normal uniformly distributed load
f = Point forces vector in global coordinates
"""
function fembeam(X1,X2,E,I,A,ro,qt,qn,f)
    le=norm(X2-X1)                      # Lenght of element
    a=atan((X2[2]-X1[2])/(X2[1]-X1[1])) # Rotation angle of the element
    nn=2                            # Number of nodes
    nd=3*1+3                           # Number of DOFs
    Gp=[-1/3*sqrt(5+2*sqrt(10/7)) -1/3*sqrt(5-2*sqrt(10/7)) 0 1/3*sqrt(5-2*sqrt(10/7)) 1/3*sqrt(5+2*sqrt(10/7))] # Five Gauss integration points
    w=[(322-13*sqrt(70))/900 (322+13*sqrt(70))/900 128/225 (322+13*sqrt(70))/900 (322-13*sqrt(70))/900] # Five Gauss integration weights
    # Rotation matrix
    B=[cos(a) sin(a) 0  0     0       0;
      -sin(a) cos(a) 0  0     0       0;
       0      0      1  0     0       0;
       0      0      0  cos(a) sin(a) 0;
       0      0      0 -sin(a) cos(a) 0;
       0      0      0  0      0      1]
    # Integration of the truss elements stiffness matrix
    detJ=2/le
    function tkint(w,xi)
        dN1=-1/2
        dN2=1/2
        dN=[dN1 dN2]
        w*dN'*dN
    end
    tk=zeros(2,2)
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
    bk=zeros(4,4)
    for i = 1:size(Gp,2)
        bk +=bkint(w[i],Gp[i])
    end
    bk=E*I*detJ*bk
    # Assembly of 6 DOF truss-beam stiffness matrix k
    k=zeros(6,6)
    k[1,1]=tk[1,1]; k[1,4]=tk[1,2]; k[4,1]=tk[2,1]; k[4,4]=tk[2,2]
    k[2,2]=bk[1,1];k[2,3]=bk[1,2];k[2,5]=bk[1,3];k[2,6]=bk[1,4]
    k[3,2]=bk[2,1];k[3,3]=bk[2,2];k[3,5]=bk[2,3];k[3,6]=bk[2,4]
    k[5,2]=bk[3,1];k[5,3]=bk[3,2];k[5,5]=bk[3,3];k[5,6]=bk[3,4]
    k[6,2]=bk[4,1];k[6,3]=bk[4,2];k[6,5]=bk[4,3];k[6,6]=bk[4,4]
    # Rotation
    k=B'*k*B
    # Integration of the truss mass matrix tm
    detJ=le/2
    function tmint(w,xi)
        N1s=1+(-1/2)*(1+xi)
        N2s=(1/2)*(1 + xi)
        Ns=[N1s N2s]
        w*Ns'*Ns
    end
    tm=zeros(2,2)
    for i = 1:size(Gp,2)
        tm +=tmint(w[i],Gp[i])
    end
    tm=ro*A*tm*detJ
    # Integration of the beam mass matrix bm
    detJ=(le/2)
    function bmint(w,xi)
        N1=1/4*(1-xi)^2*(2+xi)
        N2=le/8*(1-xi)^2*(xi+1)
        N3=1/4*(1+xi)^2*(2-xi)
        N4=le/8*(1+xi)^2*(xi-1)
        N=[N1 N2 N3 N4]
        w*N'*N
    end
    bm=zeros(4,4)
    for i = 1:size(Gp,2)
        bm +=bmint(w[i],Gp[i])
    end
    bm=ro*A*bm*detJ
    # Assembly of the 6 DOF truss-beam mass matrix m
    m=zeros(6,6)
    m[1,1]=tm[1,1];m[1,4]=tm[1,2];m[4,1]=tm[2,1];m[4,4]=tm[2,2]
    m[2,2]=bm[1,1];m[2,3]=bm[1,2];m[2,5]=bm[1,3];m[2,6]=bm[1,4]
    m[3,2]=bm[2,1];m[3,3]=bm[2,2];m[3,5]=bm[2,3];m[3,6]=bm[2,4]
    m[5,2]=bm[3,1];m[5,3]=bm[3,2];m[5,5]=bm[3,3];m[5,6]=bm[3,4]
    m[6,2]=bm[4,1];m[6,3]=bm[4,2];m[6,5]=bm[4,3];m[6,6]=bm[4,4]
    # Rotation
    m=B'*m*B
    # Integration of the equivalent forces vector
    # For truss element tfq
    detJ=le/2
    function tfqint(w,xi)
        N1s=1+(-1/2)*(1+xi)
        N2s=(1/2)*(1 + xi)
        Ns=[N1s N2s]
        w*Ns'
    end
    tfq=zeros(2,1)
    for i = 1:size(Gp,2)
        tfq +=tfqint(w[i],Gp[i])
    end
    tfq=qn*tfq*detJ
    # For 4 DOF beam element bfqe
    detJ=le/2
    function bfqint(w,xi)
        N1=1/4*(1-xi)^2*(2+xi)
        N2=le/8*(1-xi)^2*(xi+1)
        N3=1/4*(1+xi)^2*(2-xi)
        N4=le/8*(1+xi)^2*(xi-1)
        N=[N1 N2 N3 N4]
        w*N'
    end
    bfq=zeros(4,1)
    for i = 1:size(Gp,2)
        bfq +=bfqint(w[i],Gp[i])
    end
    bfq=qt*bfq*detJ # tarkista onko qn ja qt oikeissa paikoissa.
    # Assembly of the 6 DOF beam element equivalent forces vector fqe
    fqe=zeros(6,1)
    fqe[1,1],fqe[4,1]=tfq[1,1],tfq[2,1]
    fqe[2,1],fqe[3,1],fqe[5,1],fqe[6,1]=bfq[1,1],bfq[2,1],bfq[3,1],bfq[4,1]
    # Rotation
    fqe=B'*fqe
    # Adding equivalent forces vector to point forces vector
    f +=fqe
    return k,m,f
end

end
