function validateCarriers(carrier)
%   validateCarriers(CARRIER) validates the carrier properties of all
%   gNBs.

    % Validate physical layer cell identities
    cellIDs = [carrier(:).NCellID];
    if numel(cellIDs) ~= numel(unique(cellIDs))
        error('nr5g:invalidNCellID','Physical layer cell identities of all of the carriers must be unique.');
    end

    % Validate subcarrier spacings
    scsVals = [carrier(:).SubcarrierSpacing];
    if ~isscalar(unique(scsVals))
        error('nr5g:invalidSCS','Subcarrier spacing values of all of the carriers must be same.');
    end

    % Validate cyclic prefix lengths
    cpVals = {carrier(:).CyclicPrefix};
    if ~all(strcmpi(cpVals{1},cpVals))
        error('nr5g:invalidCP','Cyclic prefix lengths of all of the carriers must be same.');
    end

    % Validate NSizeGrid values
    nSizeGridVals = [carrier(:).NSizeGrid];
    if ~isscalar(unique(nSizeGridVals))
        error('nr5g:invalidNSizeGrid','NSizeGrid of all of the carriers must be same.');
    end

    % Validate NStartGrid values
    nStartGridVals = [carrier(:).NStartGrid];
    if ~isscalar(unique(nStartGridVals))
        error('nr5g:invalidNStartGrid','NStartGrid of all of the carriers must be same.');
    end

end