% Get NOAA tide Data and syn it to the ecomapper data
% Jacob Anderson
% 8/15/2019



function emData = GetTideData( stationID, emData, plotResults)

if nargin < 3
    plotResults = false;
end


% Get Tide Data -----------------------------------------------------------
times = cat(1, emData.timeStamp);

start = datetime(times(1));         % Beging of time duration
stop  = datetime(times(end));       % End of time duration

tide = NOAA_tides( stationID, start, stop, plotResults);


% Coralate with endata ----------------------------------------------------
for jj = 1:size(emData,1)
    
    for ii = 1: numel(emData(jj).timeStamp)
        
        timeDiff = abs( tide.DateTime - emData(jj).timeStamp(ii) );
        
        indx = timeDiff == min(timeDiff(:));
        
        emData(jj).tide(ii) = tide.Height(indx);
        
    end
    
    if plotResults
        figure('Name','Get Tide Data - Results', 'NumberTitle','off')
        plot(emData(jj).timeStamp, emData(jj).tide)
        xlabel('Date-Time')
        ylabel('Tide Height [m]')
    end
end
end