function prm = validateParams(prm)
% Validate user specified parameters and return updated parameters
%
% Only cross-dependent checks are made for parameter consistency.

    if strcmpi(prm.FreqRange,'FR1')
        if prm.CenterFreq > 7.125e9 || prm.CenterFreq < 410e6
            error(['Specified center frequency is outside the FR1 ', ...
                   'frequency range (410 MHz - 7.125 GHz).']);
        end
        if strcmpi(prm.SSBlockPattern,'Case D') ||  ...
           strcmpi(prm.SSBlockPattern,'Case E')
            error(['Invalid SSBlockPattern for selected FR1 frequency ' ...
                'range. SSBlockPattern must be one of ''Case A'' or ' ...
                '''Case B'' or ''Case C'' for FR1.']);
        end
        if ~((length(prm.SSBTransmitted)==4) || ...
             (length(prm.SSBTransmitted)==8))
            error(['SSBTransmitted must be a vector of length 4 or 8', ...
                   'for FR1 frequency range.']);
        end
        if (prm.CenterFreq <= 3e9) && (length(prm.SSBTransmitted)~=4)
            error(['SSBTransmitted must be a vector of length 4 for ' ...
                   'center frequency less than or equal to 3GHz.']);
        end
        if (prm.CenterFreq > 3e9) && (length(prm.SSBTransmitted)~=8)
            error(['SSBTransmitted must be a vector of length 8 for ', ...
                   'center frequency greater than 3GHz and less than ', ...
                   'or equal to 7.125GHz.']);
        end
    else % 'FR2'
        if prm.CenterFreq > 52.6e9 || prm.CenterFreq < 24.25e9
            error(['Specified center frequency is outside the FR2 ', ...
                   'frequency range (24.25 GHz - 52.6 GHz).']);
        end
        if ~(strcmpi(prm.SSBlockPattern,'Case D') || ...
                strcmpi(prm.SSBlockPattern,'Case E'))
            error(['Invalid SSBlockPattern for selected FR2 frequency ' ...
                'range. SSBlockPattern must be either ''Case D'' or ' ...
                '''Case E'' for FR2.']);
        end
        if length(prm.SSBTransmitted)~=64
            error(['SSBTransmitted must be a vector of length 64 for ', ...
                   'FR2 frequency range.']);
        end
    end

    prm.NumTx = prod(prm.TxArraySize);
    prm.NumRx = prod(prm.RxArraySize);
    if prm.NumTx==1 || prm.NumRx==1
        error(['Number of transmit or receive antenna elements must be', ...
               ' greater than 1.']);
    end
    prm.IsTxURA = (prm.TxArraySize(1)>1) && (prm.TxArraySize(2)>1);
    prm.IsRxURA = (prm.RxArraySize(1)>1) && (prm.RxArraySize(2)>1);

    if ~( strcmpi(prm.RSRPMode,'SSSonly') || ...
          strcmpi(prm.RSRPMode,'SSSwDMRS') )
        error(['Invalid RSRP measuring mode. Specify either ', ...
               '''SSSonly'' or ''SSSwDMRS'' as the mode.']);
    end

    % Select SCS based on SSBlockPattern
    switch lower(prm.SSBlockPattern)
        case 'case a'
            scs = 15;
            cbw = 10;
            scsCommon = 15;
        case {'case b', 'case c'}
            scs = 30;
            cbw = 25;
            scsCommon = 30;
        case 'case d'
            scs = 120;
            cbw = 100;
            scsCommon = 120;
        case 'case e'
            scs = 240;
            cbw = 200;
            scsCommon = 120;
    end
    prm.SCS = scs;
    prm.ChannelBandwidth = cbw;
    prm.SubcarrierSpacingCommon = scsCommon;

end