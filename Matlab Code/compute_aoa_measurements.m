function [rho, s] = compute_aoa_measurements(rxWaveform, gNBs)

    numgNBs = size(gNBs, 1);

    estimator = phased.MUSICEstimator2D(...
        'SensorArray',gNBs(1).antennaArray,...
        'OperatingFrequency',gNBs(1).fc,...
        'NumSignalsSource','Property',...
        'NumSignals',1,...
        'DOAOutputPort',true, ...
        'ForwardBackwardAveraging', true, ...
        'AzimuthScanAngles',-60:.1:60,...
        'ElevationScanAngles',-45:.1:45);
    estimated_doa = zeros(numgNBs, 2);
    s = zeros(numgNBs, 3);

    for gNBIdx = 1:numgNBs
        % estimate AoA
        [~,doa] = estimator(rxWaveform{gNBIdx});
        plotSpectrum(estimator)
        estimator.reset()

        estimated_doa(gNBIdx,1:2)=doa(:,1)'+gNBs(gNBIdx).orientation;
        estimated_doa(gNBIdx,1) = estimated_doa(gNBIdx,1)+180;
        if abs(estimated_doa(gNBIdx,1)) > 180
            estimated_doa(gNBIdx,1) = estimated_doa(gNBIdx,1)-sign(estimated_doa(gNBIdx,1))*360;
        end
        estimated_doa(gNBIdx,2) = -estimated_doa(gNBIdx, 2);
        if abs(estimated_doa(gNBIdx,2)) > 90
            estimated_doa(gNBIdx,2) = estimated_doa(gNBIdx,2)-sign(estimated_doa(gNBIdx,2))*180;
        end
        
        s(gNBIdx,:) =[gNBs(gNBIdx).getXYZ];
    end
    rho = estimated_doa.*pi/180;
end