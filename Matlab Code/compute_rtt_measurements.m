function [rho, s] = compute_rtt_measurements(ul_toa, dl_toa, gNBs)
    rtt = [];
    estimated_dist = [];
    s = [];
    rho = [];
    
    c = physconst('LightSpeed');          % speed of light
    numgNBs = size(gNBs,1);
    for gNBIdx=1:numgNBs
        delayULDL = randsrc(1,1,0:1e-9:1e-6);
    
        t0 = 0;
        t1 = t0 + ul_toa(gNBIdx);
        t2 = t1 + delayULDL; 
        t3 = t2 + dl_toa(gNBIdx); 
    
        rtt(gNBIdx) = (t3- t0)-(t2-t1);
        % rtt(gNBIdx) = dl_estimated_toa(gNBIdx)+ul_estimated_toa(gNBIdx)-2*chInfo.ChannelFilterDelay/ofdmInfo.SampleRate;
        estimated_dist(gNBIdx) = rtt(gNBIdx)*c/2;
        s = [s; gNBs(gNBIdx).getXYZ()];
        rho = [rho; estimated_dist(gNBIdx)];
        
    end
end