function [na,Sme,Hme,Tme,L] = RandlGMinClean(h1,T,l1,l2,lmin,ls,b1,Ea,Vx,Tr)
%% Initial Variables
% h1 is number of hexagon rows there are
H = h1*3^(1/2); % H is the hieght of the grid  
W = 3*H; % W is the exact width needed to fulfill 3:1 ratio
w = round(W); % w is the width - using one bond length as the counting measure 
b = w + 1; % b is the number of atoms in one of the short rows
c = b + 1; % c is the number of atoms in one of the long rows
na = b + (b + c)*h1; % na is the number of atoms
xm1 = 2*(w + 1) + 1; % xm1 is the maximum x position at the longer rows
xm2 = 2*(w + 1); % xm2 is the maximum x position at the shorter rows
xm2f = xm2 + 4;
xm1f = xm2f + 1;
ym = h1*2 + 1; % ym is the maximum y position
ymax = 2*ym - 1;
ymaxf = ymax + 2;
yf = ymaxf; % yf is the final y position
hp = 6.62606957.*10^-31; %plancks constant
EA = Ea*1.6022e-19;
k = 8.6173324*10^-5;
L = zeros(xm1f,ymaxf,2);
%% Generating random Length and Probabilities
G = exp(-b1*lmin);
for i = 2:xm2f
    L(i,1,1) = inf;
    L(i,1,2) = 0;
end
for j = 1:4:ymaxf
    L(2,j,1) = inf;
    L(2,j,2) = 0;
    L(xm2f,j,1) = inf;
    L(xm2f,j,2) = 0;
end
for j = 3:4:ymaxf
    L(3,j,1) = inf;
    L(3,j,2) = 0;
    L(xm2f - 1,j,1) = inf;
    L(xm2f - 1,j,2) = 0;
end
for i = 1: xm2f
    L(i,2,1) = inf;
    L(i,2,2) = 0;
end
for j = 4:2:ymaxf
    L(1:2,j,1) = inf;
    L(1:2,j,2) = 0;
end
for j = 4:4:ymaxf
    L(xm2f - 1:xm2f,j,1) = inf;
    L(xm2f - 1:xm2f,j,2) = 0;
end
for i = 3:(xm2f - 2)
    for j = 4:2:ymax
        L(i,j,1) = normrnd(l1,ls);
        if(L(i,j,1) <= 0)
            L(i,j,1) = 0;
        end
        if (L(i,j,1) <= l2 && L(i,j,1) >= lmin)
            L(i,j,2) = G;
        else
            L(i,j,2) = exp(-b1*L(i,j,1));
        end
    end
end

for i = 5:2:(xm2f - 2)
    for j = 3:4:ymaxf
        L(i,j) = normrnd(l1,ls);
        if(L(i,j) <= 0)
            L(i,j) = 0;
        end
        if (L(i,j,1) <= l2 && L(i,j,1) >= lmin)
            L(i,j,2) = G;
        else
            L(i,j,2) = exp(-b1*L(i,j,1));
        end
    end
end

for i = 4:2:(xm2f - 2)
    for j = 5:4:ymaxf
        L(i,j) = normrnd(l1,ls);
        if(L(i,j) <= 0)
            L(i,j) = 0;
        end
        if (L(i,j,1) <= l2 && L(i,j,1) >= lmin)
            L(i,j,2) = G;
        else
            L(i,j,2) = exp(-b1*L(i,j,1));
        end
    end
end     
%% For Loop
parfor t = 1:Tr % Tr is the number electron trials
    xi = 4 + 2*round((1 + w)*rand); % xi is the randomized initial x position
    yi = 3; % yi is the initial y position
    s = 0; % s is total hop attempts
    i = xi; % i and j are the position variables for x and y respectively
    j = yi;
    h = 0;
    t_hop = 0;
    %% Grid Section  
    while (j ~= yf) 
        g = rand; % g is the random number
        G1R = PVG(T,Ea,Vx,L(i,j + 1,2));
        G1L = PVG(T,Ea,Vx,L(i - 1,j + 1,2));
        G2R = NOVG(T,Ea,Vx,L(i + 1,j,2));
        G2L = NOVG(T,Ea,Vx,L(i - 1,j,2));
        G3R = NVG(T,Ea,Vx,L(i,j - 1,2));
        G3L = NVG(T,Ea,Vx,L(i - 1,j - 1,2));
        [i1,j1] = Ad6PR(i,j,g,G1R,G1L,G2R,G2L,G3R,G3L);
   % increasing increment     
        if (i1 ~= i || j1 ~= j)
            di = i1 - i;
            dj = j1 - j;
            h = h + 1;
            t_hop = hp*exp(Ea/(k*T) - b*L(i + di,j+dj))/(2*EA) + t_hop;
        end
        i = i1;
        j = j1;
        s = s + 1;
    end
    T_hop = t_hop(i);
    S(t) = s;
    HN(t) = h;
end
Tme = mean(T_hop);
Sme = mean(S);
Hme = mean(HN);