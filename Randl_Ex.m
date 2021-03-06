function [na, L_Ex, L_O, Sme_ex, Hme_ex, SmeStd_ex, HmeStd_ex, C_ex, G_ex, Sme, Hme, SmeStd, HmeStd, C, G] = Randl_Ex(h1, T, lm, l2, lmin, ls, l_ptc, b_alk, b_ptc, Ea_alk, Ea_ptc,Vx,Tr)
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
ym = h1*2 + 1;
ymax = 2*ym - 1;% ymax is the maximum y position
yf = ymax; % yf is the final y position
L = zeros(xm1, ymax, 2);% matrix of ligand length
L_O = zeros(xm1, ymax);
L_Ex = zeros(xm1, ymax);
M = zeros(xm1, ymax, Tr);
M_ex = zeros(xm1, ymax, Tr);
C = zeros(xm1, ymax);
C_ex = zeros(xm1, ymax);
G = zeros(xm1, ymax, 2);
G_ex = zeros(xm1, ymax, 2);

%% Generating random Length and Probabilities

for i = 1:xm2
    for j = 2:2:ymax
        L(i, j, 1) = normrnd(lm, ls);
        if(L(i, j, 1) <= 1)
            L(i, j, 2) = 1;
            L(i, j, 1) = 1;
        end
        if (L(i, j, 1) <= l2 && L(i, j, 1) >= lmin)
            L(i, j, 2) = l_ptc;
            G_ex(i, j, 1)= PVG(T, L(i, j, 2), b_ptc, Ea_ptc, Vx);
            G_ex(i, j, 2)= NVG(T, L(i, j, 2), b_ptc, Ea_ptc, Vx);
        else
            L(i, j, 2) = L(i, j, 1);
            G_ex(i, j, 1)= PVG(T, L(i, j, 2), b_alk, Ea_alk, Vx);
            G_ex(i, j, 2)= NVG(T, L(i, j, 2), b_alk, Ea_alk, Vx);
        end
        
        G(i, j, 1)= PV(T, L(i, j, 1), b_alk, Ea_alk, Vx);
        G(i, j, 2)= NV(T, L(i, j, 1), b_alk, Ea_alk, Vx);
    end
end

for i = 3:2:xm2
    for j = 1:4:ymax
        L(i, j, 1) = normrnd(lm, ls);
        if(L(i, j, 1) <= 1)
            L(i, j, 2) = 1;
            L(i, j, 1) = 1;
        end
        if (L(i, j, 1) <= l2 && L(i, j, 1) >= lmin)
            L(i, j, 2) = l_ptc;
            G_ex(i, j, 1)= NOVG(T, L(i, j, 2), b_ptc, Ea_ptc, Vx);
            G_ex(i, j, 2)= NOVG(T, L(i, j, 2), b_ptc, Ea_ptc, Vx);
        else
            L(i, j, 2) = L(i, j, 1);
            G_ex(i, j, 1)= NOVG(T, L(i, j, 2), b_alk, Ea_alk, Vx);
            G_ex(i, j, 2)= NOVG(T, L(i, j, 2), b_alk, Ea_alk, Vx);
        end
        
        G(i, j, 1)= NOV(T, L(i, j, 1), b_alk, Ea_alk, Vx);
        G(i, j, 2)= NOV(T, L(i, j, 1), b_alk, Ea_alk, Vx);
    end
end

for i = 2:2:xm1
    for j = 3:4:ymax
        L(i, j, 1) = normrnd(lm, ls);
        
        if(L(i, j, 1) <= 1)
            L(i, j, 1) = 1;
            L(i, j ,2) = 1;
        end
        
        if (L(i, j, 1) <= l2 && L(i, j, 1) >= lmin)
            L(i, j, 2) = l_ptc;
            G_ex(i, j, 1)= NOVG(T, L(i, j, 2), b_ptc, Ea_ptc, Vx);
            G_ex(i, j, 2)= NOVG(T, L(i, j, 2), b_ptc, Ea_ptc, Vx);
        else
            L(i, j, 2) = L(i, j, 1);
            G_ex(i, j, 1)= NOVG(T, L(i, j, 2), b_alk, Ea_alk, Vx);
            G_ex(i, j, 2)= NOVG(T, L(i, j, 2), b_alk, Ea_alk, Vx);
        end
        
        G(i, j, 1)= NOV(T, L(i, j, 1), b_alk, Ea_alk, Vx);
        G(i, j, 2)= NOV(T, L(i, j, 1), b_alk, Ea_alk, Vx);
    end
end


%% Main Loop for exchanged ligand
for t = 1:Tr % Tr is the number electron trials
    xi = 2*round(1 + w*rand); % xi is the randomized initial x position
    yi = 1; % yi is the initial y position
    s = 0; % s is total hop attempts
    i = xi; % i and j are the position variables for x and y respectively
    j = yi;
    h = 0;
    Temp_ex = zeros(xm1, ymax);
    Temp_ex(xi, 1) = Temp_ex(xi, 1) + 1;
    
    %% Grid Section
    while (j ~= yf)
        g = rand; % g is the random number
   % Corner
        if (i == 2 && j == 1)
            G1R = G_ex(i, j + 1, 1);
            G1L = G_ex(i - 1, j + 1, 1);
            G2R = G_ex(i + 1, j, 1);
            [i1, j1] = Ad3C1PR(i, j, g, G1R, G1L, G2R);

            if g <= G1R
             Temp_ex(i, j + 1) = Temp_ex(i, j + 1) + 1;
            elseif  (G1R < g && g <= (G1R + G1L))
             Temp_ex(i - 1, j + 1) = Temp_ex(i - 1, j + 1) + 1;
            elseif  ((G1R + G1L) < g && g <= (G1R + G1L + G2R))
             Temp_ex(i + 1, j) = Temp_ex(i + 1, j) + 1;
            end
            
             Temp_ex(i1, j1) = Temp_ex(i1, j1) + 1;

        elseif (j == 1 && i == xm2)
            G1R = G_ex(i, j + 1, 1);
            G1L = G_ex(i - 1, j + 1, 1);
            G2L = G_ex(i - 1, j, 1);
            [i1, j1] = Ad3C2PR(i, j, g, G1R, G1L, G2L);

             if g <= G1R
              Temp_ex(i, j + 1) = Temp_ex(i, j + 1) + 1;
             elseif  (G1R < g && g <= (G1R + G1L))
              Temp_ex(i - 1, j + 1) = Temp_ex(i - 1, j + 1) + 1;
             elseif  ((G1R + G1L) < g && g <= (G1R + G1L + G2L))
              Temp_ex(i - 1, j) = Temp_ex(i - 1, j) + 1;
             end
             
              Temp_ex(i1, j1) = Temp_ex(i1, j1) + 1;

    % Bottom side
        elseif j == 1
            G1R = G_ex(i, j + 1, 1);
            G1L = G_ex(i - 1, j + 1, 1);
            G2R = G_ex(i + 1, j, 1);
            G2L = G_ex(i - 1, j, 1);
            [i1, j1] = Ad4BPR(i, j, g, G1R, G1L, G2R, G2L);

            if g <= G1R
             Temp_ex(i, j + 1) = Temp_ex(i, j + 1) + 1;
             elseif  (G1R < g && g <= (G1R + G1L))
             Temp_ex(i - 1, j + 1) = Temp_ex(i - 1, j + 1) + 1;
             elseif  ((G1R + G1L) < g && g <= (G1R + G1L + G2R))
             Temp_ex(i + 1, j) = Temp_ex(i + 1, j) + 1;
             elseif  ((G1R + G1L + G2R) < g && g <= (G1R + G1L + G2R + G2L))
             Temp_ex(i - 1, j) = Temp_ex(i - 1, j) + 1;
            end
            
             Temp_ex(i1, j1) = Temp_ex(i1, j1) + 1;

    % on the sides of the longer rows
        elseif i == 1
            G1R = G_ex(i, j + 1, 1);
            G2R = G_ex(i + 1, j, 1);
            G3R = G_ex(i, j - 1, 2);
            [i1, j1] = Ad3SCNPR(i, j, g, G1R, G2R, G3R);
            if g <= G1R
             Temp_ex(i, j + 1) = Temp_ex(i, j + 1) + 1;
             elseif  (G1R < g && g <= (G1R + G2R))
             Temp_ex(i + 1, j) = Temp_ex(i + 1, j) + 1;
             elseif  ((G1R + G2R) < g && g <= (G1R + G2R + G3R))
             Temp_ex(i, j - 1) = Temp_ex(i, j - 1) + 1;
            end
             Temp_ex(i1, j1) = Temp_ex(i1, j1) + 1;

        elseif i == xm1
             G1L = G_ex(i - 1, j + 1, 1);
             G2L = G_ex(i - 1, j, 1);
             G3L = G_ex(i - 1, j - 1, 2);
             [i1,j1] = Ad3SCFPR(i, j, g, G1L, G2L, G3L);

             if g <= G1L
             Temp_ex(i - 1,j + 1) = Temp_ex(i - 1,j + 1) + 1;
             elseif  (G1L < g && g <= (G1L + G2L))
             Temp_ex(i - 1, j) = Temp_ex(i - 1, j) + 1;
             elseif  ((G1L + G2L) < g && g <= (G1L + G2L + G3L))
             Temp_ex(i - 1, j - 1) = Temp_ex(i - 1, j - 1) + 1;
            end
             Temp_ex(i1, j1) = Temp_ex(i1, j1) + 1;

    % on the sides of the shorter rows
        elseif i == 2
            G1R = G_ex(i, j + 1, 1);
            G1L = G_ex(i - 1, j + 1, 1);
            G2R = G_ex(i + 1, j, 1);
            G3R = G_ex(i, j - 1, 2);
            G3L = G_ex(i - 1, j - 1, 2);
            [i1, j1] = Ad5SCNPR(i, j, g, G1R, G1L, G2R, G3R, G3L);

            if g <= G1R
             Temp_ex(i, j + 1) = Temp_ex(i, j + 1) + 1;
             elseif  (G1R < g && g <= (G1R + G1L))
             Temp_ex(i - 1, j + 1) = Temp_ex(i - 1, j + 1) + 1;
             elseif  ((G1R + G1L) < g && g <= (G1R + G1L + G3R))
             Temp_ex(i, j - 1) = Temp_ex(i, j - 1) + 1;
             elseif  ((G1R + G1L + G3R) < g && g <= (G1R + G1L + G3R + G3L))
             Temp_ex(i - 1, j - 1) = Temp_ex(i - 1, j - 1) + 1;
             elseif  ((G1R + G1L + G3R + G3L) < g && g <= (G1R + G1L + G3R + G3L + G2R))
             Temp_ex(i + 1, j) = Temp_ex(i + 1, j) + 1;
            end
              Temp_ex(i1, j1) = Temp_ex(i1, j1) + 1;

        elseif i == xm2
            G1R = G_ex(i, j + 1, 1);
            G1L = G_ex(i - 1, j + 1, 1);
            G2L = G_ex(i - 1, j, 1);
            G3R = G_ex(i, j - 1, 2);
            G3L = G_ex(i - 1, j - 1, 2);
            [i1, j1] = Ad5SCFPR(i, j, g, G1R, G1L, G2L, G3R, G3L);

            if g <= G1R
             Temp_ex(i, j + 1) = Temp_ex(i, j + 1) + 1;
             elseif  (G1R < g && g <= (G1R + G1L))
             Temp_ex(i - 1,j + 1) = Temp_ex(i - 1,j + 1) + 1;
             elseif  ((G1R + G1L) < g && g <= (G1R + G1L + G3R))
             Temp_ex(i, j - 1) = Temp_ex(i, j - 1) + 1;
             elseif  ((G1R + G1L + G3R) < g && g <= (G1R + G1L + G3R + G3L))
             Temp_ex(i - 1, j - 1) = Temp_ex(i - 1, j - 1) + 1;
             elseif  ((G1R + G1L + G3R + G3L) < g && g <= (G1R + G1L + G3R + G3L + G2L))
             Temp_ex(i - 1, j) = Temp_ex(i - 1, j) + 1;
            end
              Temp_ex(i1, j1) = Temp_ex(i1, j1) + 1;

    % Elsewhere on the grid
        else
            G1R = G_ex(i, j + 1, 1);
            G1L = G_ex(i - 1, j + 1, 1);
            G2R = G_ex(i + 1, j, 1);
            G2L = G_ex(i - 1, j, 1);
            G3R = G_ex(i, j - 1, 2);
            G3L = G_ex(i - 1, j - 1, 2);
            [i1,j1] = Ad6PR(i, j, g, G1R, G1L, G2R, G2L, G3R, G3L);

            if g <= G1R
             Temp_ex(i, j + 1) = Temp_ex(i, j + 1) + 1;
             elseif  (G1R < g && g <= (G1R + G1L))
             Temp_ex(i - 1, j + 1) = Temp_ex(i - 1, j + 1) + 1;
             elseif  ((G1R + G1L) < g && g <= (G1R + G1L + G3R))
             Temp_ex(i, j - 1) = Temp_ex(i, j - 1) + 1;
             elseif  ((G1R + G1L + G3R) < g && g <= (G1R + G1L + G3R + G3L))
             Temp_ex(i - 1,j - 1) = Temp_ex(i - 1, j - 1) + 1;
             elseif  ((G1R + G1L + G3R + G3L) < g && g <= (G1R + G1L + G3R + G3L + G2R))
             Temp_ex(i + 1, j) = Temp_ex(i + 1, j) + 1;
             elseif  ((G1R + G1L + G3R + G3L + G2R) < g && g <= (G1R + G1L + G3R + G3L + G2R + G2L))
             Temp_ex(i - 1, j) = Temp_ex(i - 1, j) + 1;
            end
              Temp_ex(i1, j1) = Temp_ex(i1, j1) + 1;

        end
    % increasing increment
        if (i1 ~= i && j1 ~= j)
            h = h + 1;
        end
        i = i1;
        j = j1;
        s = s + 1;
    % failsafe break - to stop program from running too long
        if s == 1000000000
            s = 1000000000;
            break
        end
    end
    S_ex(t) = s;
    HN_ex(t) = h;
    M_ex(:, :, t) = Temp_ex(:, :);
end

%% output
Sme_ex = mean(S_ex);
Hme_ex = mean(HN_ex);
C_ex = sum(M_ex, 3);
SmeStd_ex =std(S_ex);
HmeStd_ex =std(HN_ex);



%% Main Loop for no exchanged NanoArray
for t = 1:Tr % Tr is the number electron trials
    h = 0;
    s = 0; % s is total hop attempts
    xi = 2*round(1 + w*rand); % xi is the randomized initial x position
    yi = 1;
    i = xi; % i and j are the position variables for x and y respectively
    j = yi;
    Temp = zeros(xm1, ymax);
    Temp(xi, 1) = Temp(xi, 1) + 1;

    %% Grid Section
    while (j ~= yf)
        g = rand; % g is the random number
   % Corner
        if (i == 2 && j == 1)
            G1R = G(i, j + 1, 1);
            G1L = G(i - 1, j + 1, 1);
            G2R = G(i + 1, j, 1);
            [i1, j1] = Ad3C1PR(i, j, g, G1R, G1L, G2R);

            if g <= G1R
             Temp(i, j + 1) = Temp(i, j + 1) + 1;
            elseif  (G1R < g && g <= (G1R + G1L))
             Temp(i - 1, j + 1) = Temp(i - 1, j + 1) + 1;
            elseif  ((G1R + G1L) < g && g <= (G1R + G1L + G2R))
             Temp(i + 1, j) = Temp(i + 1, j) + 1;
            end
            
             Temp(i1, j1) = Temp(i1, j1) + 1;

        elseif (j == 1 && i == xm2)
            G1R = G(i, j + 1, 1);
            G1L = G(i - 1, j + 1, 1);
            G2L = G(i - 1, j, 1);
            [i1, j1] = Ad3C2PR(i, j, g, G1R, G1L, G2L);

             if g <= G1R
              Temp(i, j + 1) = Temp(i, j + 1) + 1;
             elseif  (G1R < g && g <= (G1R + G1L))
              Temp(i - 1, j + 1) = Temp(i - 1, j + 1) + 1;
             elseif  ((G1R + G1L) < g && g <= (G1R + G1L + G2L))
              Temp(i - 1, j) = Temp(i - 1, j) + 1;
             end
             
              Temp(i1, j1) = Temp(i1, j1) + 1;

    % Bottom side
        elseif j == 1
            G1R = G(i, j + 1, 1);
            G1L = G(i - 1, j + 1, 1);
            G2R = G(i + 1, j, 1);
            G2L = G(i - 1, j, 1);
            [i1, j1] = Ad4BPR(i, j, g, G1R, G1L, G2R, G2L);

            if g <= G1R
             Temp(i, j + 1) = Temp(i, j + 1) + 1;
             elseif  (G1R < g && g <= (G1R + G1L))
             Temp(i - 1, j + 1) = Temp(i - 1, j + 1) + 1;
             elseif  ((G1R + G1L) < g && g <= (G1R + G1L + G2R))
             Temp(i + 1, j) = Temp(i + 1, j) + 1;
             elseif  ((G1R + G1L + G2R) < g && g <= (G1R + G1L + G2R + G2L))
             Temp(i - 1, j) = Temp(i - 1, j) + 1;
            end
            
             Temp(i1, j1) = Temp(i1, j1) + 1;

    % on the sides of the longer rows
        elseif i == 1
            G1R = G(i, j + 1, 1);
            G2R = G(i + 1, j, 1);
            G3R = G(i, j - 1, 2);
            [i1, j1] = Ad3SCNPR(i, j, g, G1R, G2R, G3R);
            if g <= G1R
             Temp(i, j + 1) = Temp(i, j + 1) + 1;
             elseif  (G1R < g && g <= (G1R + G2R))
             Temp(i + 1, j) = Temp(i + 1, j) + 1;
             elseif  ((G1R + G2R) < g && g <= (G1R + G2R + G3R))
             Temp(i, j - 1) = Temp(i, j - 1) + 1;
            end
             Temp(i1, j1) = Temp(i1, j1) + 1;

        elseif i == xm1
             G1L = G(i - 1, j + 1, 1);
             G2L = G(i - 1, j, 1);
             G3L = G(i - 1, j - 1, 2);
             [i1, j1] = Ad3SCFPR(i, j, g, G1L, G2L, G3L);

             if g <= G1L
             Temp(i - 1, j + 1) = Temp(i - 1, j + 1) + 1;
             elseif  (G1L < g && g <= (G1L + G2L))
             Temp(i - 1, j) = Temp(i - 1, j) + 1;
             elseif  ((G1L + G2L) < g && g <= (G1L + G2L + G3L))
             Temp(i - 1, j - 1) = Temp(i - 1, j - 1) + 1;
            end
             Temp(i1, j1) = Temp(i1, j1) + 1;

    % on the sides of the shorter rows
        elseif i == 2
            G1R = G(i, j + 1, 1);
            G1L = G(i - 1, j + 1, 1);
            G2R = G(i + 1, j, 1);
            G3R = G(i, j - 1, 2);
            G3L = G(i - 1, j - 1, 2);
            [i1, j1] = Ad5SCNPR(i, j, g, G1R, G1L, G2R, G3R, G3L);

            if g <= G1R
             Temp(i, j + 1) = Temp(i, j + 1) + 1;
             elseif  (G1R < g && g <= (G1R + G1L))
             Temp(i - 1, j + 1) = Temp(i - 1, j + 1) + 1;
             elseif  ((G1R + G1L) < g && g <= (G1R + G1L + G3R))
             Temp(i, j - 1) = Temp(i, j - 1) + 1;
             elseif  ((G1R + G1L + G3R) < g && g <= (G1R + G1L + G3R + G3L))
             Temp(i - 1, j - 1) = Temp(i - 1, j - 1) + 1;
             elseif  ((G1R + G1L + G3R + G3L) < g && g <= (G1R + G1L + G3R + G3L + G2R))
             Temp(i + 1, j) = Temp(i + 1, j) + 1;
            end
              Temp(i1, j1) = Temp(i1, j1) + 1;

        elseif i == xm2
            G1R = G(i, j + 1, 1);
            G1L = G(i - 1, j + 1, 1);
            G2L = G(i - 1, j, 1);
            G3R = G(i, j - 1, 2);
            G3L = G(i - 1, j - 1, 2);
            [i1, j1] = Ad5SCFPR(i, j, g, G1R, G1L, G2L, G3R, G3L);

            if g <= G1R
             Temp(i, j + 1) = Temp(i, j + 1) + 1;
             elseif  (G1R < g && g <= (G1R + G1L))
             Temp(i - 1, j + 1) = Temp(i - 1, j + 1) + 1;
             elseif  ((G1R + G1L) < g && g <= (G1R + G1L + G3R))
             Temp(i, j - 1) = Temp(i, j - 1) + 1;
             elseif  ((G1R + G1L + G3R) < g && g <= (G1R + G1L + G3R + G3L))
             Temp(i - 1, j - 1) = Temp(i - 1, j - 1) + 1;
             elseif  ((G1R + G1L + G3R + G3L) < g && g <= (G1R + G1L + G3R + G3L + G2L))
             Temp(i - 1, j) = Temp(i - 1, j) + 1;
            end
              Temp(i1, j1) = Temp(i1, j1) + 1;

    % Elsewhere on the grid
        else
            G1R = G(i, j + 1, 1);
            G1L = G(i - 1, j + 1, 1);
            G2R = G(i + 1, j, 1);
            G2L = G(i - 1, j, 1);
            G3R = G(i, j - 1, 2);
            G3L = G(i - 1, j - 1, 2);
            [i1, j1] = Ad6PR(i, j, g, G1R, G1L, G2R, G2L, G3R, G3L);

            if g <= G1R
             Temp(i, j + 1) = Temp(i, j + 1) + 1;
             elseif  (G1R < g && g <= (G1R + G1L))
             Temp(i - 1, j + 1) = Temp(i - 1, j + 1) + 1;
             elseif  ((G1R + G1L) < g && g <= (G1R + G1L + G3R))
             Temp(i, j - 1) = Temp(i, j - 1) + 1;
             elseif  ((G1R + G1L + G3R) < g && g <= (G1R + G1L + G3R + G3L))
             Temp(i - 1, j - 1) = Temp(i - 1, j - 1) + 1;
             elseif  ((G1R + G1L + G3R + G3L) < g && g <= (G1R + G1L + G3R + G3L + G2R))
             Temp(i + 1, j) = Temp(i + 1, j) + 1;
             elseif  ((G1R + G1L + G3R + G3L + G2R) < g && g <= (G1R + G1L + G3R + G3L + G2R + G2L))
             Temp(i - 1, j) = Temp(i - 1, j) + 1;
            end
              Temp(i1, j1) = Temp(i1, j1) + 1;

        end
        
        % increasing increment
        if (i1 ~= i && j1 ~= j)
            h = h + 1;
        end
        i = i1;
        j = j1;
        s = s + 1;
    % failsafe break - to stop program from running too long
        if s == 1000000000
            s = 1000000000;
            break
        end
    end
    
    S(t) = s;
    HN(t) = h;
    M(:, :, t) = Temp(:, :);
end

%% OUTPUT
Sme = mean(S);
Hme = mean(HN);
C = sum(M, 3);
SmeStd=std(S);
HmeStd=std(HN);

for x = 1: xm1
    for y = 1:ymax
        L_Ex(x,y) = L(x,y,2);
        L_O(x,y) = L(x,y,1);
    end
end
