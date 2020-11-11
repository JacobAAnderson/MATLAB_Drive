%% Use WEBREAD to compile a list of Wikipedia pages
% Final example script shown in the "Reading Web Pages, Part 1: Using
% webread" video on the "Stuart's MATLAB Videos" blog.

% Copyright 2015 The MathWorks, Inc.

%% Read a random page
txt = webread('http://en.wikipedia.org/wiki/Special:Random');
% Find its own name
pg = regexp(txt,'(?<="canonical" href="https://en.wikipedia.org/wiki/)[^:"\s]+','match')

%% Extract the first link
pg{2,1} = findfirstlink(txt)

%% Surf until the philosophy page is reached
while ~isequal(pg{end},'Philosophy')
    txt = webread(['http://en.wikipedia.org/wiki/',pg{end}]);
    pg{end+1,1} = findfirstlink(txt)
end

%% Collect into a table
phd = table(pg,((length(pg)-1):-1:0)','VariableNames',{'Page','Steps'})
