function [i,j] = Ad3C2PR(i,j,g,G1R,G1L,G2L)
if (G1R + G1L + G2L) <= 1
    
    if g <= G1R
    i = i + 1;
    j = j + 2;
   
    elseif  (G1R < g && g <= (G1R + G1L))
    j = j + 2;
    i = i - 1;
    
    elseif  ((G1R + G1L) < g && g <= (G1R + G1L + G2L))
    j = j;
    i = i - 2;
    
    else
    j = j;
    i = i;

    end
    
else
    
    if g <= G1R/(G1R + G1L + G2L)
    i = i + 1;
    j = j + 2;
   
    elseif  (G1R/(G1R + G1L + G2L) < g && g <= (G1R + G1L)/(G1R + G1L + G2L))
    j = j + 2;
    i = i - 1;
    
    elseif  ((G1R + G1L)/(G1R + G1L + G2L) < g && g <= 1)
    j = j;
    i = i - 2;
    
    else
    j = j;
    i = i;
    
    end
    
end