%% CRB for AoA measurements
function lb = CRB_AOA(SNR, N, d, BW, angles, H)
    c = physconst("LightSpeed");
%     lb = (3*lambda)/(4*pi^2*SNR*d^2*(sin(angle))^2*(N-1)*N*(2*N-1));
    lb = sqrt(3)*c./(sqrt(2*SNR*N*(N^(2)-1))*BW*pi*d*cos(angles));
    R = diag(lb.^2);
    lb = trace(inv(H'*(R)^(-1)*H));
end