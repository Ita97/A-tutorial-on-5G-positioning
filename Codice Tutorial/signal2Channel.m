function rxWaveform = signal2Channel(downlink, TxPower, txWaveform, channel, sampleDelay, numgNBs, los, N0, Nfft, NSizeGrid) 
    chInfo = info(channel{1});
    numRxAnt = chInfo.NumReceiveAntennas;
    numTxAnt = chInfo.NumTransmitAntennas;
    maxChDelay = max(cell2mat(sampleDelay))+ chInfo.ChannelFilterDelay;
    rx = cell(1,numgNBs);
    rxWaveform = zeros(length(txWaveform{1}) + maxChDelay,numRxAnt);
    TxAmp = db2mag(TxPower-30)*sqrt(Nfft^2/(NSizeGrid*12*numTxAnt));
    

    for gNBIdx = 1:numgNBs
        rx{gNBIdx}=zeros(length(txWaveform{gNBIdx}) + maxChDelay,numRxAnt);
        channel{gNBIdx}.ChannelFiltering = false;
        [pathGains, ~] = channel{gNBIdx}();
        pg = permute(pathGains,[2 1 3 4]); % first dimension is the number of paths

        % Amplitude of the tx signal
        txWaveform{gNBIdx} = txWaveform{gNBIdx}*TxAmp;
    
        % in LOS cases sum the first to paths, they correspond to the LOS ray
        if los(gNBIdx)
            pg = [sum(pg(1:2,:,:,:)); pg(3:end,:,:,:)];
            sampleDelay{gNBIdx} = sampleDelay{gNBIdx}(2:end);
        end
        pg = abs(pg);
    
        % instantaneus channel gain
        
        for j = 1: numRxAnt
            for i = 1: numTxAnt
                for path = 1:length(sampleDelay{gNBIdx})
                    rx{gNBIdx}(:,j) = rx{gNBIdx}(:,j)+[zeros(sampleDelay{gNBIdx}(path),1); txWaveform{gNBIdx}(:,i); ...
                    zeros(maxChDelay-sampleDelay{gNBIdx}(path),1)]*pg(path,1,i,j);
                end
            end
        end
    
        % add noise
        noise = N0*complex(randn(size(rx{gNBIdx})),randn(size(rx{gNBIdx})));
        %SF_noise = 10.^(sigma_SF*complex(randn(size(rx{gNBIdx})),randn(size(rx{gNBIdx})))/10);
    
        rx{gNBIdx} = rx{gNBIdx}+noise;    
        
        % Sum waveforms from all gNBs
        rxWaveform = rxWaveform + rx{gNBIdx};
    end


    if ~downlink
        rxWaveform = rx;
    end
end