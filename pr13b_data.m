clc
clear all
close all


Mz = 120
Mn = 60
nn = 200
kv = 2
z1 = 17
z2 = 96
omegan = 2*pi*nn/60

am = Mz 
km = (Mz-Mn)/omegan

D = sqrt(km^2 - 4*am*kv*(z1/z2)^3)

rooty = [-kv*(z1/z2)^3 -km am]

omegaMS = roots(rooty)

omegaVS = omegaMS*z1/z2

Q = -km - 2*kv*(z1/z2)^3*omegaMS 






















