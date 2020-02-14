function [i,j] = Ad3SCNPR(i,j,g,G1R,G2R,G3R)
if (G1R + G2R + G3R) <= 1
    if g <= G1R
    i = i + 1;
    j = j + 2;
   
    elseif  (G1R < g && g <= (G1R + G2R))
    j = j;
    i = i + 2;
    
    elseif  ((G1R + G2R) < g && g <= (G1R + G2R + G3R))
    j = j - 2;
    i = i + 1;
    
    else
    j = j;
    i = i;
    
    end
    
else
    if g <= G1R/(G1R + G2R + G3R)
    i = i + 1;
    j = j + 2;
   
    elseif  (G1R/(G1R + G2R + G3R) < g && g <= (G1R + G2R)/(G1R + G2R + G3R))
    j = j;
    i = i + 2;
    
    elseif  ((G1R + G2R)/(G1R + G2R + G3R) < g && g <= 1)
    j = j - 2;
    i = i + 1;
    
    else
    j = j;
    i = i;
    
    end
    
end