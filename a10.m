
m = 300;
J  = 4;
D = 400e-3;
d = 50e-3;
L = 1.2;
E = 2e11;
mu = 0.3;
Sl = 1.5e-4;
Ll = 8;
El = 1.6e10;
G = E/2/(1+mu)

Kl = Sl*El/Ll
Kh = pi*d^4/32 * G / L

M = [J,0;
     0,m]

K = [Kh + Kl*D^2/4, -Kl*D/2;
     -Kl*D/2,       Kl]

[A,Lam] = eig(inv(M)*K) 
f = Lam.^(1/2)./(2*pi)
