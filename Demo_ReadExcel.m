% Demo Read Excel Sheet
% Jacob Anderson
% Sept 8, 2020

clear all
close all
clc


% Read in Excel Sheet
[~,~,raw] = xlsread('/Users/jake/Desktop/demo.xlsx');

% ---- Clean up data -----
for ii = 2: size(raw,1)
    
    % -- Format phone number --
    phone = 'xxx-xxx-xxxx';
    ind = 1;
    
    number = raw{ii,2};
    
    % check if the number was read in as double
    if isa(number, 'double') 
        number = sprintfc('%d',number); 
        number = number{1};
    end
    
    for num = number
        
        switch num
            case {'(',')','-',' '}, continue, 
            otherwise, phone(ind) = num;
        end
        
        if     ind == 3, ind = 5;
        elseif ind == 7, ind = 9;
        else,            ind = ind+1;
        end
        
    end
    
    disp(phone)
    
    raw{ii,2} = phone;
    
end





%% Concatinate Text
combinedTxt = [raw(2:end,:);raw(2:end,:)];

fprintf('\n\nCombined Excel Data:\n\n')
disp(combinedTxt)

% Instanciate Table
sz = size(combinedTxt);
varTypes = repmat({'string'}, 1, sz(2));
varNames = strtrim(raw(1,:));

T = table('Size',sz,'VariableTypes',varTypes,'VariableNames',varNames);


% Copy text data to table
for ii = 1:sz(1), T(ii,:) = combinedTxt(ii,:); end


% Get rid of duplicate entries
T = unique(T);

fprintf('\n\nCombined data after removing duplicates:\n\n')
disp(T)

% Write Table to a new excel sheet
writetable(T, '/Users/jake/Desktop/demo_out.xlsx')