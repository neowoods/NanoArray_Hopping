function [na,Sme,Hme,L, SmeStd, HmeStd, Matrix_ave, G] = Randl(h1,T,lm,ls,Ea,Vx,Tr)
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
ym = h1*2 + 1; % ym is the maximum y position
ymax = 2*ym - 1;
yf = ymax; % yf is the final y position
L = zeros(xm1,ymax);
M = zeros(xm1, ymax);
G = zeros(xm1, ymax, 2);
Matrix_sum = zeros(xm1, ymax);

%% Generating random Length and Conductivity Probabilities
for i = 1:xm2
    for j = 2:2:ymax
        L(i,j) = normrnd(lm,ls);
        if(L(i,j) <= 0)
            L(i,j) = 0;
        end
        G(i,j,1)= PV(T,L(i,j),Ea,Vx);
        G(i,j,2)= NV(T,L(i,j),Ea,Vx);
    end
end

for i = 3:2:xm2
    for j = 1:4:ymax
        L(i,j) = normrnd(lm,ls);
        if(L(i,j) <= 0)
            L(i,j) = 0;
        end
        G(i,j,1)= NOV(T,L(i,j),Ea,Vx);
        G(i,j,2)= NOV(T,L(i,j),Ea,Vx);
    end
end

for i = 2:2:xm1
    for j = 3:4:ymax
        L(i,j) = normrnd(lm,ls);
        if(L(i,j) <= 0)
            L(i,j) = 0;
        end
        G(i,j,1)= NOV(T,L(i,j),Ea,Vx);
        G(i,j,2)= NOV(T,L(i,j),Ea,Vx);
    end
end 


%% Main Loop
for t = 1:Tr % Tr is the number electron trials
    h = 0;
    s = 0; % s is total hop attempts
    xi = 2*round(1 + w*rand); % xi is the randomized initial x position
    yi = 1;
    i = xi; % i and j are the position variables for x and y respectively
    j = yi;
%% Grid Section  
    while (j ~= yf) 
        g = rand; % g is the random number
   % Corner
        if (i == 2 && j == 1)
            G1R = PV(T,L(i,j + 1),Ea,Vx);
            G1L = PV(T,L(i - 1,j + 1),Ea,Vx);
            G2R = NOV(T,L(i + 1,j),Ea,Vx);
            [i1,j1] = Ad3C1PR(i,j,g,G1R,G1L,G2R);
            if g <= G1R
             M(i,j+1) = M(i,j+1) + 1;
            elseif  (G1R < g && g <= (G1R + G1L))
             M(i-1,j+1) = M(i-1,j+1) + 1;
            elseif  ((G1R + G1L) < g && g <= (G1R + G1L + G2R))
             M(i + 1,j) = M(i + 1,j) + 1;
            end
             M(i1,j1) = M(i1,j1)+1;
             
        elseif (j == 1 && i == xm2)
            G1R = PV(T,L(i,j + 1),Ea,Vx);
            G1L = PV(T,L(i - 1,j + 1),Ea,Vx);
            G2L = NOV(T,L(i - 1,j),Ea,Vx);
            [i1,j1] = Ad3C2PR(i,j,g,G1R,G1L,G2L);
             if g <= G1R
              M(i,j+1) = M(i,j+1) + 1;
             elseif  (G1R < g && g <= (G1R + G1L))
              M(i - 1,j + 1) = M(i - 1,j + 1) + 1;
             elseif  ((G1R + G1L) < g && g <= (G1R + G1L + G2L))
              M(i - 1,j) = M(i - 1,j) + 1;
             end
              M(i1,j1) = M(i1,j1)+1;
              
    % Bottom side
        elseif j == 1
            G1R = PV(T,L(i,j + 1),Ea,Vx);
            G1L = PV(T,L(i - 1,j + 1),Ea,Vx);
            G2R = NOV(T,L(i + 1,j),Ea,Vx);
            G2L = NOV(T,L(i - 1,j),Ea,Vx);
            [i1,j1] = Ad4BPR(i,j,g,G1R,G1L,G2R,G2L);
            if g <= G1R
             M(i,j+1) = M(i,j+1) + 1;
             elseif  (G1R < g && g <= (G1R + G1L))
             M(i - 1,j + 1) = M(i - 1,j + 1) + 1;
             elseif  ((G1R + G1L) < g && g <= (G1R + G1L + G2R))
             M(i + 1,j) = M(i + 1,j) + 1;
             elseif  ((G1R + G1L + G2R) < g && g <= (G1R + G1L + G2R + G2L))
             M(i - 1,j) = M(i - 1,j) + 1;    
            end
              M(i1,j1) = M(i1,j1)+1;
              
    % on the sides of the longer rows
        elseif i == 1
            G1R = PV(T,L(i,j + 1),Ea,Vx);
            G2R = NOV(T,L(i + 1,j),Ea,Vx);
            G3R = NV(T,L(i,j - 1),Ea,Vx);
            [i1,j1] = Ad3SCNPR(i,j,g,G1R,G2R,G3R);
            if g <= G1R
             M(i,j+1) = M(i,j+1) + 1;
             elseif  (G1R < g && g <= (G1R + G2R))
             M(i + 1,j) = M(i + 1,j) + 1;
             elseif  ((G1R + G2R) < g && g <= (G1R + G2R + G3R))
             M(i,j - 1) = M(i,j - 1) + 1;
            end
             M(i1,j1) = M(i1,j1)+1;
             
        elseif i == xm1
             G1L = PV(T,L(i - 1,j + 1),Ea,Vx);
             G2L = NOV(T,L(i - 1,j),Ea,Vx);
             G3L = NV(T,L(i - 1,j - 1),Ea,Vx);
            [i1,j1] = Ad3SCFPR(i,j,g,G1L,G2L,G3L);
            
             if g <= G1L
             M(i - 1,j + 1) = M(i - 1,j + 1) + 1;
             elseif  (G1L < g && g <= (G1L + G2L))
             M(i - 1,j) = M(i - 1,j) + 1;
             elseif  ((G1L + G2L) < g && g <= (G1L + G2L + G3L))
             M(i - 1,j - 1) = M(i - 1,j - 1) + 1;
            end
             M(i1,j1) = M(i1,j1)+1;
             
    % on the sides of the shorter rows
        elseif i == 2
            G1R = PV(T,L(i,j + 1),Ea,Vx);
            G1L = PV(T,L(i - 1,j + 1),Ea,Vx);
            G2R = NOV(T,L(i + 1,j),Ea,Vx);
            G3R = NV(T,L(i,j - 1),Ea,Vx);
            G3L = NV(T,L(i - 1,j - 1),Ea,Vx);
            [i1,j1] = Ad5SCNPR(i,j,g,G1R,G1L,G2R,G3R,G3L);
            if g <= G1R
             M(i,j+1) = M(i,j+1) + 1;
             elseif  (G1R < g && g <= (G1R + G1L))
             M(i - 1,j + 1) = M(i - 1,j + 1) + 1;
             elseif  ((G1R + G1L) < g && g <= (G1R + G1L + G2R))
             M(i + 1,j) = M(i + 1,j) + 1;
             elseif  ((G1R + G1L + G2R) < g && g <= (G1R + G1L + G2R + G3R))
             M(i,j - 1) = M(i,j - 1) + 1;
             elseif  ((G1R + G1L + G2R + G3R) < g && g <= (G1R + G1L + G2R + G3R + G3L))
             M(i - 1,j - 1) = M(i - 1,j - 1) + 1;
            end
              M(i1,j1) = M(i1,j1)+1;
              
        elseif i == xm2 
            G1R = PV(T,L(i,j + 1),Ea,Vx);
            G1L = PV(T,L(i - 1,j + 1),Ea,Vx);
            G2L = NOV(T,L(i - 1,j),Ea,Vx);
            G3R = NV(T,L(i,j - 1),Ea,Vx);
            G3L = NV(T,L(i - 1,j - 1),Ea,Vx);
            [i1,j1] = Ad5SCFPR(i,j,g,G1R,G1L,G2L,G3R,G3L);
            if g <= G1R
             M(i,j+1) = M(i,j+1) + 1;
             elseif  (G1R < g && g <= (G1R + G1L))
             M(i - 1,j + 1) = M(i - 1,j + 1) + 1;
             elseif  ((G1R + G1L) < g && g <= (G1R + G1L + G2L))
             M(i - 1,j) = M(i - 1,j) + 1;
             elseif  ((G1R + G1L + G2L) < g && g <= (G1R + G1L + G2L + G3R))
             M(i,j - 1) = M(i,j - 1) + 1;
             elseif  ((G1R + G1L + G2L + G3R) < g && g <= (G1R + G1L + G2L + G3R + G3L))
             M(i - 1,j - 1) = M(i - 1,j - 1) + 1;
            end
              M(i1,j1) = M(i1,j1)+1;
              
    % Elsewhere on the grid
        else
            G1R = PV(T,L(i,j + 1),Ea,Vx);
            G1L = PV(T,L(i - 1,j + 1),Ea,Vx);
            G2R = NOV(T,L(i + 1,j),Ea,Vx);
            G2L = NOV(T,L(i - 1,j),Ea,Vx);
            G3R = NV(T,L(i,j - 1),Ea,Vx);
            G3L = NV(T,L(i - 1,j - 1),Ea,Vx);
            [i1,j1] = Ad6PR(i,j,g,G1R,G1L,G2R,G2L,G3R,G3L);
            if g <= G1R
             M(i,j+1) = M(i,j+1) + 1;
             elseif  (G1R < g && g <= (G1R + G1L))
             M(i - 1,j + 1) = M(i - 1,j + 1) + 1;
             elseif  ((G1R + G1L) < g && g <= (G1R + G1L + G2R))
             M(i + 1,j) = M(i + 1,j) + 1;
             elseif  ((G1R + G1L + G2R) < g && g <= (G1R + G1L + G2R + G2L))
             M(i - 1,j) = M(i - 1,j) + 1;
             elseif  ((G1R + G1L + G2R + G2L) < g && g <= (G1R + G1L + G2R + G2L + G3R))
             M(i,j - 1) = M(i,j - 1) + 1;
             elseif  ((G1R + G1L + G2R + G2L + G3R) < g && g <= (G1R + G1L + G2R + G2L + G3R + G3L))
             M(i - 1,j - 1) = M(i - 1,j - 1) + 1;    
            end
              M(i1,j1) = M(i1,j1)+1;
              
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
    Matrix_sum = Matrix_sum + M;
end

%% OUTPUT
Sme = mean(S);
Hme = mean(HN);
Matrix_ave = Matrix_sum;
SmeStd=std(S);
HmeStd=std(HN);

