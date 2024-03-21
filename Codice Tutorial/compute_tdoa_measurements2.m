%% TDOA measurement computation
function [rho, s, curveX, curveY]=compute_tdoa_measurements2(available_gNBs, toa, previous_est)
    c = physconst('LightSpeed');
    numgNBs = size(available_gNBs,1);
    max_nrho=0;
    rhos = cell(numgNBs,1);
    ss = cell(numgNBs,1);
    curves = cell(numgNBs,2);
    rstdVals = zeros(length(toa));
    for jj = 1:length(toa)
        for ii = 1:length(toa)
            rstdVals(ii,jj) = toa(ii)*c - toa(jj)*c;
        end
    end
    
    if ismethod(available_gNBs, 'getXYZ')
        available_gNBs = [available_gNBs(:).getXYZ();];
    end
    for jj = 1:numgNBs
        s = [];
        rho = [];
        curveX = {};
        curveY = {};
        cellIdx = 1;
        %s(cellIdx,:) = available_gNBs(jj).getXYZ();
        s(cellIdx,:) = available_gNBs(jj,:);

        for ii = [1:(jj-1) jj+1:numgNBs]
            rstd = rstdVals(ii,jj); % Delay distance
            % Establish gNBs for which delay distance is applicable by
            % examining detected cell identities
           
            %[x,y] = getRSTDCurve(available_gNBs(ii).getXYZ(),available_gNBs(jj).getXYZ(),rstd); 
            [x,y] = getRSTDCurve(available_gNBs(ii,:),available_gNBs(jj,:),rstd); 

            if isreal(x) && isreal(y)
                curveX{1,cellIdx} = x; %#ok<*SAGROW> 
                curveY{1,cellIdx} = y;
    
                % Get the gNB numbers corresponding to the current hyperbola curve
                gNBNums{cellIdx} = [jj ii];
                rho(cellIdx) = rstd;
                cellIdx = cellIdx + 1;
                %s(cellIdx,:) = available_gNBs(ii).getXYZ();
                s(cellIdx,:) = available_gNBs(ii,:);
            end
                  
            
        end
        if max_nrho<length(rho)
            max_nrho=length(rho);
        end
        ss{jj}=s;
        rhos{jj}=rho;
        curves{jj,1}=curveX;
        curves{jj,2}=curveY;
        
    end
    best_conf=[];
    for i = 1:numgNBs
        if length(rhos{i})==max_nrho
            best_conf = [best_conf i];
        end
    end

    if nargin>2
        index = 0;
        best_HDOP = inf;
        for i = 1:min(3,length(best_conf))
            j = best_conf(i);
            H=compute_H(previous_est,ss{j},3,'TDOA');
            C=(H'*H)^(-1);
            d=diag(C);
            GDOP = sqrt(sum(d.^2));
            if best_HDOP>GDOP
                best_HDOP = GDOP;
                index = j;
            end
        end
    else
        index=best_conf(1);
    end
    if index>0
        rho = rhos{index};
        s = ss{index};
        curveX = curves{index,1};
        curveY = curves{index,2};
    else
        s = [];
        rho = [];
        curveX = {};
        curveY = {};
    end
end