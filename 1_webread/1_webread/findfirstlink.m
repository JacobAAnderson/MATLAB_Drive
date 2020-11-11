% Finds first link in the main text of a Wikipedia page
%
% PG = FINDFIRSTLINK(TXT) parses the html of a Wikipedia page in the char
% array TXT and returns a string (char array) PG containing the name of the
% page of the first link in the main text of the page TXT (not in italics
% or inside parentheses).

% Copyright 2015 The MathWorks, Inc.

function pg = findfirstlink(txt)
% Remove all tables, accounting for possible nesting
% Score all characters in txt as being in a table or not
p = zeros(size(txt));
% Add 1 at the beginning of each new table, and subtract 1 at the end
p(regexp(txt,'<table')) = 1;
p(regexp(txt,'</table')) = -1;
% Cumulative sum => level of nesting (how many tables a given character is
% within). Remove anything above 0.
txt(cumsum(p)>0) = [];

% Now find the beginning of the text
% Remove anything before the first h1 header
[~,idx] = regexp(txt,'<h1.*?<p>','once');
txt(1:(idx-3)) = [];
% Look for the last new paragraph (<p>) tag before the first bold (<b>)
% tag. First bold tag should correspond to the definition of the subject.
% (Occasionally the first link is before the bolded definition, hence go
% back to the beginning of that paragraph.)
idx = regexp(txt,'<b>','once');
idx = regexp(txt(1:idx),'<p>');
txt(1:idx(end)) = [];
% txt now starts at the main text

% Find first link not in parentheses. However, can't just delete everything
% in parentheses first because many links themselves contain parens (eg
% "Property_(philosophy)").
% Calculate level of parenthetical nesting
p = cumsum((txt == '(') - (txt == ')'));
% Find wiki page links
[x,y] = regexp(txt,'<a href\s*=\s*"/wiki/([^:"\s]+)"','start','tokens');
% Take the first link that isn't within parentheses (ie starts where p = 0)
pg = y{find(~p(x),1,'first')}{1};
