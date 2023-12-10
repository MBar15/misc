clc
clear all
close all


aps = 200
bps = 75
i = 5
am = 75
bm = 2.5
cm = -0.0625


rooty = [cm bm-bps/i^2 am-aps/i]

omegaMS = roots(rooty)

dMn = bm+2*cm*20
dMp = bps/i 

Q = dMn - dMp  %stabilita