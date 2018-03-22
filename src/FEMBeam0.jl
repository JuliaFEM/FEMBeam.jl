# Element stiffness matrice
# Number of elements ne,number of nodes nn, number of dofs nd
l=2.0
E=210.0e9
I=50/100^4
A=14/100^2
ro=1
ne=1
BCs=[1 , 2 , 3] # DOFs 1,2 and 3
qt=0
qn=-5.0e3
F=-0.0e3
le=l/ne
nn=ne+1
nd=3*ne+3 # kommentti
# Local stiffness and mass matrices
Gp=[-sqrt(3/5) 0 sqrt(3/5)] ; w=[5/9 8/9 5/9] # Gaussian points and weights
detJ=(2/le)^3
function Gauss3(w,ks) w*
    [0.25*A/I 0 0 -0.25*A/I  0 0;
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
# equivalent forcevector
fql=[qt*le/2;      # x
     qn*le/2;      # y
     qn*le^2/12;   # M
     qt*le/2;      # x
     qn*le/2;      # y
    -qn*le^2/12]   # M
fq=zeros(nd)
fq[1:3]=fql[1:3]
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
U=K\f
# Printing out displacement v from the right end of the beam
println("comsolin vertailuarvo: -0.095238")
println("                Julia: ", U[end-3-1]) # tulostetaan palkin pään pystysiirtymä
