

m = 200;
Jb = 12;
i = 4;
Ll = 8;
Sl = 1e-4;
El = 2e10;

Lh = 1.6;
dh = 80e-3;
Eh = 2e11;
mu=0.3
Gh = Eh/2/(1+mu)

Rb = 300e-3;

Kh = pi*dh^4/32 * Gh / Lh
Kl = El*Sl/Ll


M = [Jb,0; 0, m]

K = [Kh*i^2 + Kl*Rb^2, -Kl*Rb;
     -Kl*Rb,  Kl]


[A,Lam] = eig(inv(M)*K) 
f = Lam.^(1/2)./(2*pi)

inv(M)*K

%2.76521e5 + 0.125000e5 
%2.76521e5*0.125000e5  + 0.62500e5*-0.0375000e5
