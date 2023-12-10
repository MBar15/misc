clc
clear all
close all

syms omega

Mzv = 150;
ns = 750;
nzv = 550;
Mps = 648;
i = 7.2;

omegas = 2*pi*ns/60
omegazv = 2*pi*nzv/60

s = (omegas-omega)/omegas 
szv = (omegas-omegazv)/omegas 

Mm = (2*Mzv)/(s/szv + szv/s)

f = Mps/i == Mm

omega1 = vpa(solve(f, omega))

dMps = 0
dMm = diff(Mm,omega)

dMm = vpa(subs(dMm, omega1))


