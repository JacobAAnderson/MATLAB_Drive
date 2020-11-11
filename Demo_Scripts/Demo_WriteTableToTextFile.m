%% Demo: Write a tabel to a text file
clc

fileName = 'Table.txt';


% Generate Table
LastName      = {'Smith';'Johnson';'Williams';'Jones';'Brown'};
Age           = [38;43;38;40;49];
Height        = [71;69;64;67;64];
Weight        = [176;163;131;133;119];
BloodPressure = [124 93; 109 77; 125 83; 117 75; 122 80];

T = table( Age, Height, Weight, BloodPressure,...
    'RowNames',LastName);


% Write table to text file
writetable(T,fileName,    ...
    'WriteRowNames',true, ...
    'Delimiter', '\t')  


% Display the text output
type Table.txt