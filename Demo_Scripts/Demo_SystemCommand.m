% Demo: Exicute system Command
clc

command = './a.out';                                    % Comand to Exicute
filePath = '/home/jake/C++/Demo_isam/';                 % Location of the Exicutible


currentDir = pwd;                                       % Get Current directory to come back to

cd( filePath ) 
[status,cmdout] = system(command,'-echo');              % Do exicutable

cd( currentDir )                                        % Go back to the origonal directory


TF = isstrprop(cmdout,'cntrl');                         % Find control charatures
delim = cmdout(TF);

newStr = split(cmdout,delim(1));                        % Sepperate the output string

in = cellfun(@isempty,newStr);                          % Get rid of empty cells
newStr(in) = [];

fprintf('\n\n\n Out Put:\n\n')                          % Display results
disp(newStr);