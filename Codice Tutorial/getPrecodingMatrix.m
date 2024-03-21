function wtx = getPrecodingMatrix(carrier,pdsch,hestGrid,prgbundlesize)
% Calculate precoding matrices for all PRGs in the carrier that overlap
% with the PDSCH allocation
    
    % Maximum CRB addressed by carrier grid
    maxCRB = carrier.NStartGrid + carrier.NSizeGrid - 1;
    
    % PRG size
    if nargin==4 && ~isempty(prgbundlesize)
        Pd_BWP = prgbundlesize;
    else
        Pd_BWP = maxCRB + 1;
    end
    
    % PRG numbers (1-based) for each RB in the carrier grid
    NPRG = ceil((maxCRB + 1) / Pd_BWP);
    prgset = repmat((1:NPRG),Pd_BWP,1);
    prgset = prgset(carrier.NStartGrid + (1:carrier.NSizeGrid).');
    
    [~,~,R,P] = size(hestGrid);
    wtx = zeros([pdsch.NumLayers P NPRG]);
    for i = 1:NPRG
    
        % Subcarrier indices within current PRG and within the PDSCH
        % allocation
        thisPRG = find(prgset==i) - 1;
        thisPRG = intersect(thisPRG,pdsch.PRBSet(:) + carrier.NStartGrid,'rows');
        prgSc = (1:12)' + 12*thisPRG';
        prgSc = prgSc(:);
       
        if (~isempty(prgSc))
        
            % Average channel estimate in PRG
            estAllocGrid = hestGrid(prgSc,:,:,:);
            Hest = permute(mean(reshape(estAllocGrid,[],R,P)),[2 3 1]);

            % SVD decomposition
            [~,~,V] = svd(Hest);
            wtx(:,:,i) = V(:,1:pdsch.NumLayers).';

        end
    
    end
    
    wtx = wtx / sqrt(pdsch.NumLayers); % Normalize by NumLayers

end