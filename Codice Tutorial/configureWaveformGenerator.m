function cfgDL = configureWaveformGenerator(prm,txBurst)
% Configure an nrDLCarrierConfig object to be used by nrWaveformGenerator
% to generate the SS burst waveform.

    cfgDL = nrDLCarrierConfig;
    cfgDL.SCSCarriers{1}.SubcarrierSpacing = prm.SCS;
    if (prm.SCS==240)
        cfgDL.SCSCarriers = [cfgDL.SCSCarriers cfgDL.SCSCarriers];
        cfgDL.SCSCarriers{2}.SubcarrierSpacing = prm.SubcarrierSpacingCommon;
        cfgDL.BandwidthParts{1}.SubcarrierSpacing = prm.SubcarrierSpacingCommon;
    else
        cfgDL.BandwidthParts{1}.SubcarrierSpacing = prm.SCS;
    end
    cfgDL.PDSCH{1}.Enable = false;
    cfgDL.PDCCH{1}.Enable = false;
    cfgDL.ChannelBandwidth = prm.ChannelBandwidth;
    cfgDL.FrequencyRange = prm.FreqRange;
    cfgDL.NCellID = prm.NCellID;
    cfgDL.NumSubframes = 5;
    cfgDL.WindowingPercent = 0;
    cfgDL.SSBurst = txBurst;

end