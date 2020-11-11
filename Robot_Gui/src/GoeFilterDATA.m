function DATA = GoeFilterDATA(DATA, polygon)

for ii = 1: length(DATA)
    
    Data = DATA{ii};
    
    dataPoints = cat(1, Data.vehicle);                           % Extrace data waypointds from the data structure for visulizations
    inPoints = inpolygon( dataPoints(:,2), dataPoints(:,1) , polygon(:,1), polygon(:,2));
    
    outPoints = ~inPoints;
    
    for jj = 1: length(Data)
        
        if outPoints(jj)
        Data(jj).vehicle = [];
        Data(jj).date    = [];
        Data(jj).bathy   = [];
        Data(jj).wqData  = [];
        end
        
    end
    DATA{ii} = Data;
    
end

end
