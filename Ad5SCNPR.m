function [i,j] = Ad5SCNPR(i,j,g,G1R,G1L,G2R,G3R,G3L)
if (G1R + G1L + G2R + G3R + G3L) <= 1
    if g <= G1R 
    i = i + 1;
    j = j + 2;
    
    elseif  (G1R < g && g <= (G1R + G1L))
    j = j + 2;
    i = i - 1;
    
    elseif  ((G1R + G1L) < g && g <= (G1R + G1L + G3R))
    j = j - 2;
    i = i + 1;
    
    elseif  ((G1R + G1L + G3R) < g && g <= (G1R + G1L + G3R + G3L))
    j = j - 2;
    i = i - 1;
    
    elseif  ((G1R + G1L + G3R + G3L) < g && g <= (G1R + G1L + G3R + G3L + G2R))
    j = j;
    i = i + 2;
    
    else
    j = j;
    i = i;
    
    end

else 
    if g <= G1R/(G1R + G1L + G2R + G3R + G3L)
    i = i + 1;
    j = j + 2;
    
    elseif  (G1R/(G1R + G1L + G2R + G3R + G3L) < g && g <= (G1R + G1L)/(G1R + G1L + G2R + G3R + G3L))
    j = j + 2;
    i = i - 1;
    
    elseif  ((G1R + G1L)/(G1R + G1L + G2R + G3R + G3L) < g && g <= (G1R + G1L + G3R)/(G1R + G1L + G2R + G3R + G3L))
    j = j - 2;
    i = i + 1;
    
    elseif  ((G1R + G1L + G3R)/(G1R + G1L + G2R + G3R + G3L) < g && g <= (G1R + G1L + G3R + G3L)/(G1R + G1L + G2R + G3R + G3L))
    j = j - 2;
    i = i - 1;
    
    elseif  ((G1R + G1L + G3R + G3L)/(G1R + G1L + G2R + G3R + G3L) < g && g <= 1)
    j = j;
    i = i + 2;
    
    else
    j = j;
    i = i;
    
    end
    
end