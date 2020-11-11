%% Use WEBREAD to follow Wikipedia pages
% Initial example script shown in the "Reading Web Pages, Part 1: Using
% webread" video on the "Stuart's MATLAB Videos" blog.

% Copyright 2015 The MathWorks, Inc.

%% Read a random page
txt = webread('http://en.wikipedia.org/wiki/Special:Random');
% Find its own name
pg = regexp(txt,'(?<="canonical" href="https://en.wikipedia.org/wiki/)[^:"\s]+','match');
pg = pg{1}

%% Extract the first link
pg = findfirstlink(txt)

%% Surf until the philosophy page is reached
while ~isequal(pg,'Philosophy')
    txt = webread(['http://en.wikipedia.org/wiki/',pg]);
    pg = findfirstlink(txt)
end
