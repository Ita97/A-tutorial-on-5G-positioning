function [prsGrid, DLdataGrid] = PRSResourceGeneration(prs, pdsch, cell_carrier, numTxAnt, nFrames, precoding)
    if nargin < 6
        precoding = 0;
    end
    numgNBs = numel(cell_carrier);
    totSlots = nFrames*cell_carrier{1}.SlotsPerFrame;
    DLdataGrid = cell(numgNBs,1);
    prsGrid = cell(numgNBs,1);
    wtx = cell(1,numgNBs);

    for slotIdx = 0:totSlots-1 
        [prsSym,prsInd] = deal(cell(1,numgNBs));
        for gNBIdx = 1:numgNBs
            if iscell(precoding)
                % MIMO precoding
                wtx{gNBIdx} = precoding{gNBIdx}; 
            end

            cell_carrier{gNBIdx}.NSlot = slotIdx;

            % Create an empty resource grid spanning one slot in time domain
            slotGrid = nrResourceGrid(cell_carrier{gNBIdx},numTxAnt);

            % Generate PRS symbols and indices
            prsSym{gNBIdx} = nrPRS(cell_carrier{gNBIdx},prs(gNBIdx));
            prsInd{gNBIdx} = nrPRSIndices(cell_carrier{gNBIdx},prs(gNBIdx));

            % Map PRS resources to slot grid
            slotGrid(prsInd{gNBIdx}) = prsSym{gNBIdx};
            prsGrid{gNBIdx} = [prsGrid{gNBIdx} slotGrid];

            % Transmit data in slots in which the PRS is not transmitted by any of
            % the gNBs (to control the hearability problem)
            dataSlotGrid = nrResourceGrid(cell_carrier{gNBIdx},numTxAnt);
            if all(cellfun(@isempty,prsInd))

                % Generate PDSCH indices
                [pdschInd,pdschInfo] = nrPDSCHIndices(cell_carrier{gNBIdx},pdsch(gNBIdx));

                % Generate random data bits for transmission
                data = randi([0 1],pdschInfo.G,1);

                % Generate PDSCH symbols
                pdschSym = nrPDSCH(cell_carrier{gNBIdx},pdsch(gNBIdx),data);
                if iscell(precoding)
                    [pdschAntSym,pdschAntInd] = hPRGPrecode(size(slotGrid),...
                        cell_carrier{gNBIdx}.NStartGrid,pdschSym,pdschInd,wtx{gNBIdx});
                end

                % Generate demodulation reference signal (DM-RS) indices and symbols
                dmrsInd = nrPDSCHDMRSIndices(cell_carrier{gNBIdx},pdsch(gNBIdx));
                dmrsSym = nrPDSCHDMRS(cell_carrier{gNBIdx},pdsch(gNBIdx));
                if iscell(precoding)
                    [dmrsAntSym,dmrsAntInd] = hPRGPrecode(size(slotGrid),...
                        cell_carrier{gNBIdx}.NStartGrid,dmrsSym,dmrsInd,wtx{gNBIdx});  
                end

                % Map PDSCH and its associated DM-RS to slot grid
                if iscell(precoding)
                    dataSlotGrid(pdschAntInd) = pdschAntSym;
                    dataSlotGrid(dmrsAntInd) = dmrsAntSym;
                else
                    dataSlotGrid(pdschInd) = pdschSym;
                    dataSlotGrid(dmrsInd) = dmrsSym;
                end
            end
            DLdataGrid{gNBIdx} = [DLdataGrid{gNBIdx} dataSlotGrid];
                
        end
    end
