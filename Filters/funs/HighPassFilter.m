%% Hi pass filter

function DataOut = HighPassFilter( DataIn, std_theshold)
DataOut = DataIn;
r = 5;

[m,n] = size(DataIn);

pOut = zeros(1,n);
for j = (r+1): m-(r+1)               % Cycle through the data set.
    % Reassing values that are "std_theshold" greater than the STD to be the average value

    s = std(DataIn(j-r:j+r, :));
    u = mean(DataIn(j-r:j+r, :));

    out = abs(DataIn(j,:)-u) > std_theshold*s;

    pOut = pOut + out;
    
    DataOut(j,out)  = u(out);
    DataOut(j,~out) = DataIn(j,~out);
end

fprintf('Number of points filtered out: ')
disp(pOut)
end
