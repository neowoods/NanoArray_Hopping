function [i,j] = Ad3SCFPR(i,j,g,G1L,G2L,G3L)
if (G1L + G2L + G3L) <= 1
    if g <= G1L
    i = i - 1;
    j = j + 2;
   
    elseif  (G1L < g && g <= (G1L + G2L))
    j = j;
    i = i - 2;
    
    elseif  ((G1L + G2L) < g && g <= (G1L + G2L + G3L))
    j = j - 2;
    i = i - 1;
    
    else
    j = j;
    i = i;
    
    end
    
else
    if g <= G1L/(G1L + G2L + G3L)
    i = i - 1;
    j = j + 2;
   
    elseif  (G1L/(G1L + G2L + G3L) < g && g <= (G1L + G2L)/(G1L + G2L + G3L))
    j = j;
    i = i - 2;
    
    elseif  ((G1L + G2L)/(G1L + G2L + G3L) < g && g <= 1)
    j = j - 2;
    i = i - 1;
    
    else
    j = j;
    i = i;
    
    end
    
end