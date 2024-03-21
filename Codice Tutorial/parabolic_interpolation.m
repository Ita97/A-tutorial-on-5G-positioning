%% Parabolic Interpolation

function [a, b, c, max] = parabolic_interpolation(corr, index)

    R_max = corr(index);
    R_r = corr(index+1);
    R_l = corr(index-1);
    
    a = (R_l+R_r)/2-R_max;
    b = (R_r-R_l)/2;
    c = R_max;
    
    max=-b/(2*a);

end

