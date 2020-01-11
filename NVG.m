function G = NVG(T, L, b, Ea, Vx)
q = 1.60217657.*10^-19;%Coulombs
h= 6.62606957.*10^-31;%m^2kg/s
G0 = 5.5*10^3;%2*q^2./h;
kb = 8.6173324*10^-5; %eV/K
G = G0*exp(-b*L - Vx - 2*Ea/(kb*T));