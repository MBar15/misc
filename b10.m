

Jp = 4;
Jh = 320;
Lm = 2;
dm = 80e-3;
Lh = 4;
dh= 100e-3;
E = 2e11;
mu = 0.3
i = 8

G = E/2/(1+mu)
Ktm = pi*dm^4/32  * G   / Lm
Kth = pi*dh^4/32  * G   / Lh


M = [Jp,0;
      0, Jh]

K = [Ktm+Kth/i^2, -Kth/i;
     -Kth/i,       Kth]

[A,D] = eig(inv(M)*K)

F = D.^(1/2)./(2*pi)

inv(M)*K

%3.9403e4 + 0.0590e4
%3.9403e4*0.0590e4 + -1*-0.5900e4*-0.0074e4
%syms x

%a = x^2 + (3.9403e4 + 0.0590e4)*x -0.5900e4*0.0074e4 +  3.9403e4 * 0.0590e4== 0
%xx = vpa(solve(a,x))