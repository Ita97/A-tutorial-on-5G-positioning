function [rho, s] = toa_estimation(downlink, rxWaveform, refGrid, gNBs, gNBs_idx, carrier, ofdmInfo, reflections_order, interpolation, plot)
    numgNBs = size(gNBs,1);
    corr = cell(1,numgNBs);
    delayEst = zeros(1,numgNBs);
    estimated_toa = zeros(1,numgNBs);
    maxCorr = zeros(1,numgNBs);
    max_index = zeros(1,numgNBs);
    s=[];
    rho=[];
    delays=[];
    for gNBIdx = 1:numgNBs
        % REFGRID
        % K-by-N-by-P, 
        % K: number of subcarriers
        % N: number of OFDM symbols
        % P: number of reference signal ports
        if downlink
            % Extract correlation data samples spanning about 1/14 ms for normal
            % cyclic prefix and about 1/12 ms for extended cyclic prefix (this
            % truncation is to ignore noisy side lobe peaks in correlation outcome)
            [~,mag] = nrTimingEstimate(carrier{gNBIdx},rxWaveform,refGrid{gNBIdx});
            corr{gNBIdx} = mag(1:(ofdmInfo.Nfft*carrier{gNBIdx}.SubcarrierSpacing/15));
        else
            [~,mag] = nrTimingEstimate(carrier,rxWaveform{gNBIdx},refGrid{gNBIdx});
            corr{gNBIdx} = mag(1:(ofdmInfo.Nfft*carrier.SubcarrierSpacing/15));
        end
        if sum(corr{gNBIdx})==0
            warning('No correlation!')
            continue;
        end
        if reflections_order == 0
            % return the max value and the index
            [maxCorr(gNBIdx), max_index(gNBIdx)] = max(corr{gNBIdx}); 
        else
            % use first peak methods
            [maxCorr(gNBIdx), max_index(gNBIdx)] = first_arrival(corr{gNBIdx}, 2); 
        end
        if interpolation && max_index(gNBIdx)>1
            [~, ~, ~, delay_offset] = parabolic_interpolation(corr{gNBIdx}, max_index(gNBIdx));
            delays = [delays; delay_offset/ofdmInfo.SampleRate];
        else
            delay_offset = 0;
        end

        delayEst(gNBIdx) = find(corr{gNBIdx} == maxCorr(gNBIdx),1)-1+delay_offset;
        estimated_toa(gNBIdx) = delayEst(gNBIdx)/ofdmInfo.SampleRate;
        s = [s; gNBs(gNBIdx).xyz_cartesian];
        rho = [rho; estimated_toa(gNBIdx)];
        
        
    end

    if plot
        % plot correlations
        plotPRSCorr(corr,ofdmInfo.SampleRate, gNBs_idx);
    end
    
end