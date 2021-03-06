function [i,j] = Ad4BPR(i,j,g,G1R,G1L,G2R,G2L)
if (G1R + G1L + G2R + G2L) <= 1
    if g <= G1R
    i = i + 1;
    j = j + 2;
   
    elseif  (G1R < g &&  g <= (G1R + G1L))
    j = j + 2;
    i = i - 1;
    
    elseif  ((G1R + G1L) < g && g <= (G1R + G1L + G2R))
    j = j;
    i = i + 2;
    
    elseif  ((G1R + G1L + G2R) < g && g <= (G1R + G1L + G2R + G2L))
    j = j;
    i = i - 2;
    
    else
    j = j;
    i = i;
    
    end
    
else
    if g <= G1R/(G1R + G1L + G2R + G2L)
    i = i + 1;
    j = j + 2;
   
    elseif  (G1R/(G1R + G1L + G2R + G2L) < g &&  g <= (G1R + G1L)/(G1R + G1L + G2R + G2L))
    j = j + 2;
    i = i - 1;
    
    elseif  ((G1R + G1L)/(G1R + G1L + G2R + G2L) < g && g <= (G1R + G1L + G2R)/(G1R + G1L + G2R + G2L))
    j = j;
    i = i + 2;
    
    elseif  ((G1R + G1L + G2R)/(G1R + G1L + G2R + G2L) < g && g <= 1)
    j = j;
    i = i - 2;
    
    else
    j = j;
    i = i;
    
    end
    
end