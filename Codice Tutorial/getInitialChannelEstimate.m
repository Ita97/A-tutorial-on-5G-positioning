function estChannelGrid = getInitialChannelEstimate(carrier,nTxAnts,propchannel)
% Obtain channel estimate before first transmission. This can be used to
% obtain a precoding matrix for the first slot.

    ofdmInfo = nrOFDMInfo(carrier);
    
    chInfo = info(propchannel);
    maxChDelay = ceil(max(chInfo.PathDelays*propchannel.SampleRate)) + chInfo.ChannelFilterDelay;
    
    % Temporary waveform (only needed for the sizes)
    tmpWaveform = zeros((ofdmInfo.SampleRate/1000/carrier.SlotsPerSubframe)+maxChDelay,nTxAnts);
    
    % Filter through channel    
    [~,pathGains,sampleTimes] = propchannel(tmpWaveform);
    
    % Perfect timing synch    
    pathFilters = getPathFilters(propchannel);
    offset = nrPerfectTimingEstimate(pathGains,pathFilters);
    
    % Perfect channel estimate
    estChannelGrid = nrPerfectChannelEstimate(carrier,pathGains,pathFilters,offset,sampleTimes);
    
end