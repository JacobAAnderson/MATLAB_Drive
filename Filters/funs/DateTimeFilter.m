%% Filter to filter out Erronius dates and times


function inPut = DateTimeFilter(inPut)

[m,~] = size(inPut);
steps = inPut(2:m) - inPut(1:m-1);
step = median(steps);
dev = std(steps);

if dev < abs(step)
    fprintf('Date and Time appeare to be fine\n  Median incrament: %s \tStandard Deveation: %s\n\n',step,dev)
    return
end

fprintf('Correcting Date-time data')
pOut = 0;

for ii = 2:m
    
    if abs(inPut(ii)-inPut(ii-1)) > abs(step*3)
        inPut(ii) = inPut(ii-1) + step;
        pOut = pOut+1;
      %  fprintf('Cerrection on line %d:  %s \n',ii,inPut(ii-1) + step)
    end
    
end

fprintf(' --> Number of corrected points: %d\n\n',pOut)

end