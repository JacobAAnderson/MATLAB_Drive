% Query Copernicus Open Access Hub
% for info --> https://scihub.copernicus.eu/
%          --> https://scihub.copernicus.eu/dhus/#/home

clc
close all
clear all


searchUrl = 'https://scihub.copernicus.eu/dhus';

% query = '/search?q=*';  % General Query, returns list of products

% Area inside a polygon: POLYGON((P1Lon P1Lat, P2Lon P2Lat, ..., PnLon PnLat))
% query = '/search?q=footprint:"Intersects(POLYGON((-130.0 49.0,-130.0 45.0,-123.0 45.0,-123.0 49.0,-130.0 49.0)))"'; 
% query = '/search?q=footprint:"Intersects(41.40338, 2.17403)"';  % Point ( lat, lon)

% start = '&start=0';
% rows  = '&rows=100';

query = '/search?q=platformname:Sentinel-2';

options = weboptions('Username','jacob_anderson',...
                     'Password','meccu2-nuvFeq-gedcuw');

data = webread([searchUrl,query], options);

% data = webread(url,options);

disp(data)


%% Pars XML, maybe???



% %"simple" java code to create a document from said string
% xmlDocument = javax.xml.parsers.DocumentBuilderFactory.newInstance().newDocumentBuilder.parse(java.io.StringBufferInputStream(data));
% 
% %"intuitive" methods to explore the xmlDocument
% nodeList = xmlDocument.getElementsByTagName('volume');
% numberOfNodes = nodeList.getLength();
% 
% firstNode = nodeList.item(0);
% firstNodeContent = firstNode.getTextContent;
% 
% disp(firstNodeContent);  %Returns '256'
