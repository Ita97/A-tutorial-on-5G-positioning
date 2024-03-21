function [srsGrid, ULdataGrid] = SRSResourceGeneration(srs, pusch,ue_carrier, numUEAnt, nFrames, precoding_or_numgNBs)
    if iscell(precoding_or_numgNBs)
        precoding = precoding_or_numgNBs;
        numgNBs = numel(precoding);
    else
        numgNBs = precoding_or_numgNBs;
        precoding = 0;
    end
    totSlots = nFrames*ue_carrier.SlotsPerFrame;
    ULdataGrid = cell(numgNBs,1);
    srsGrid = cell(numgNBs,1);
    wtx = cell(1,numgNBs);

    for slotIdx = 0:totSlots-1 
        ue_carrier.NSlot = slotIdx;
        [srsSym,srsInd] = deal(cell(1,numgNBs));
        [srsAntSym,srsAntInd] = deal(cell(1,numgNBs));
    
        for gNBIdx = 1:numgNBs
            
            if iscell(precoding)
                % MIMO precoding
                wtx{gNBIdx} = precoding{gNBIdx}; 
            end
    
            % Create an empty resource grid spanning one slot in time domain
            slotGrid = nrResourceGrid(ue_carrier,numUEAnt);
    
            % Generate SRS symbols and indices
            srsSym{gNBIdx} = nrSRS(ue_carrier,srs);
            srsInd{gNBIdx} = nrSRSIndices(ue_carrier,srs);
            if iscell(precoding)
                [srsAntSym{gNBIdx},srsAntInd{gNBIdx}] = hPRGPrecode(size(slotGrid),...
                    ue_carrier.NStartGrid,srsSym{gNBIdx},srsInd{gNBIdx},wtx{gNBIdx});

                % Map SRS resources to slot grid
                slotGrid(srsAntInd{gNBIdx}) = srsAntSym{gNBIdx};
            else
                slotGrid(srsInd{gNBIdx}) = srsSym{gNBIdx};
            end
    
            srsGrid{gNBIdx} = [srsGrid{gNBIdx} slotGrid];
    
            dataSlotGrid = nrResourceGrid(ue_carrier,numUEAnt);
            if all(cellfun(@isempty,srsAntInd))

                % Generate PUSCH indices
                [puschInd,puschInfo] = nrPUSCHIndices(ue_carrier,pusch);
    
                % Generate random data bits for transmission
                data = randi([0 1],puschInfo.G,1);
              
                % Generate PUSCH symbols
                puschSym = nrPUSCH(ue_carrier,pusch,data);
                if iscell(precoding)
                    [puschAntSym,puschAntInd] = hPRGPrecode(size(slotGrid),...
                        ue_carrier.NStartGrid,puschSym,puschInd,wtx{gNBIdx});
                end
    
                % Generate demodulation reference signal (DM-RS) indices and symbols
                dmrsInd = nrPUSCHDMRSIndices(ue_carrier,pusch);
                dmrsSym = nrPUSCHDMRS(ue_carrier,pusch);
                if iscell(precoding)
                    [dmrsAntSym,dmrsAntInd] = hPRGPrecode(size(slotGrid),...
                        ue_carrier.NStartGrid,dmrsSym,dmrsInd,wtx{gNBIdx});  
                end
                
                % Map PUSCH and its associated DM-RS to slot grid
                if iscell(precoding)
                    dataSlotGrid(puschAntInd) = puschAntSym;
                    dataSlotGrid(dmrsAntInd) = dmrsAntSym;
                else
                    dataSlotGrid(puschInd) = puschSym;
                    dataSlotGrid(dmrsInd) = dmrsSym;
                end
            end
            ULdataGrid{gNBIdx} = [ULdataGrid{gNBIdx} dataSlotGrid];
                
        end
    end
