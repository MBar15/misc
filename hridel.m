clc;
clear all;
rho = 7850; % kg/m^3

r20 = 10e-3; % m
r24 = 12e-3; % m

l60 = 60e-3; % m
l30 = 30e-3; % m
l10 = 10e-3; % m

% hmotnost
mV20x60 = r20^2 * pi * l60 * rho 
mV24x30 = r24^2 * pi * l30 * rho
mV24x10 = r24^2 * pi * l10 * rho 

mK80x30x30 = 80e-3 * 30e-3 * 30e-3 * rho 

% x
% valec = 1/2 * m * r
IV20x60 = 1/2 * mV20x60 * r20^2
IV24x30 = 1/2 * mV24x30 * r24^2 + mV24x30 * (40e-3)^2
IV24x10 = 1/2 * mV24x10 * r24^2

% kostka 1/12 * m * (a^2 + b^2)
IK80x30x30 = 1/12 * mK80x30x30 * ((80e-3)^2 + (30e-3)^2) + mK80x30x30*(20e-3)^2

Ix = IV20x60*2 + IV24x30*2 + IV24x10 * 2 + IK80x30x30*4 + 1/2 * mV24x30 * r24^2
Ix * 1e12
mass = mV20x60*2 + mV24x30*3 + mV24x10*2 + mK80x30x30*4