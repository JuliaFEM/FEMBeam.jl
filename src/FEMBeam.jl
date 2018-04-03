# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMBeam.jl/blob/master/LICENSE

""" Beam implementation for JuliaFEM. """
module FEMBeam

using FEMBase

"""
    fembeam(l,E,I,A,ro,ne,BCs,qt,qn,F)

Function integrates
l=          # lenght of the beam
E=          # Young's modulus
I=          # Moment of inertia
A=          # Cross section area
ro=         # Density
ne=         # Number of elements
BCs=        # Boundary conditions, for example with [1,2,3] displacements 1,2 and 3 =0
qt=         # Vertical Uniform load
qn=         # Horizontal uniform load
F=
"""
function fembeam(l,E,I,A,ro,ne,BCs,qt,qn,F)
le=l/ne     # Lenght of element
nn=ne+1     # Number of nodes
nd=3*ne+3   # Number of DOFs
# Local stiffness and mass matrices
Gp=[-sqrt(3/5) 0 sqrt(3/5)] ; w=[5/9 8/9 5/9] # Gaussian points and weights
detJ=(2/le)^3 # Jacobian ## Truss dofs should have different detJ ?
function Gauss3(w,ks) # This function and the for loop after this does the integration
    w*[0.25*A/I 0 0 -0.25*A/I  0 0;
    0 (-3/2 + (3/2)*(1 + ks))^2 (-3/2 + (3/2)*(1 + ks))*(le*(-1 + (1/2)*(1 + ks)) + (1/4)*le*(1 + ks)) 0 (-3/2 + (3/2)*(1 + ks))*(-(1 + ks) + (1/2)*(3 - (1 + ks))) (-3/2 + (3/2)*(1 + ks))*((1/2)*le*(-1 + (1/2)*(1 + ks)) + (1/2)*le*(1 + ks));
    0 (-3/2 + (3/2)*(1 + ks))*(le*(-1 + (1/2)*(1 + ks)) + (1/4)*l*(1 + ks)) (le*(-1 + (1/2)*(1 + ks)) + (1/4)*le*(1 + ks))^2 0 (-(1 + ks) + (1/2)*(3 - (1 + ks)))*(le*(-1 + (1/2)*(1 + ks)) + (1/4)*le*(1 + ks)) (le*(-1 + (1/2)*(1 + ks)) + (1/4)*le*(1 + ks))*((1/2)*le*(-1 + (1/2)*(1 + ks)) + (1/2)*le*(1 + ks));
     -0.25*A/I 0 0 0.25*A/I 0 0;
    0 (-3/2 + (3/2)*(1 + ks))*(-(1 + ks) + (1/2)*(3 - (1 + ks))) (-(1 + ks) + (1/2)*(3 - (1 + ks)))*(le*(-1 + (1/2)*(1 + ks)) + (1/4)*le*(1 + ks)) 0 (-(1 + ks) + (1/2)*(3 - (1 + ks)))^2 (-(1 + ks) + (1/2)*(3 - (1 + ks)))*((1/2)*le*(-1 + (1/2)*(1 + ks)) + (1/2)*le*(1 + ks));
    0  (-3/2 + (3/2)*(1 + ks))*((1/2)*le*(-1 + (1/2)*(1 + ks)) + (1/2)*le*(1 + ks)) (le*(-1 + (1/2)*(1 + ks)) + (1/4)*le*(1 + ks))*((1/2)*le*(-1 + (1/2)*(1 + ks)) + (1/2)*le*(1 + ks)) 0 (-(1 + ks) + (1/2)*(3 - (1 + ks)))*((1/2)*le*(-1 + (1/2)*(1 + ks)) + (1/2)*le*(1 + ks)) ((1/2)*le*(-1 + (1/2)*(1 + ks)) + (1/2)*le*(1 + ks))^2]
end
k=zeros(6,6)
for i =1:size(Gp,2)
    k += Gauss3(w[i],Gp[i])
end
k=k*detJ*E*I
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
# Massmatrix ## I didnt get the integration working with this
m = ro*A*le/6*[2       0          0            1      0             0;
               0       156/70     22*le/70     0      54/70        -13*le/70;
               0       22*le/70   4*le^2/70    0      13*le/70     -3*le^2/70;
               1       0          0            2      0             0;
               0       54/70      13*le/70     0      156/70       -22*le/70;
               0      -13*le/70  -3*le^2/70    0     -22*le/70      4*le^2/70]
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
#
U=K\f
return U,K
end

end
