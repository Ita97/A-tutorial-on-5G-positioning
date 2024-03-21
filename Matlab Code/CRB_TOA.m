%% Cramer Rao Bound for ranging measurements
function lb = CRB_TOA(SNR,  BW, H)
    % SNR - lineare
    % BW - Hz
    c = physconst("LightSpeed");
    lb = (8*pi^2*BW^2*SNR)^(-1);
    lb = c*sqrt(lb);
    lb = trace(lb^2./(H'*H));
end