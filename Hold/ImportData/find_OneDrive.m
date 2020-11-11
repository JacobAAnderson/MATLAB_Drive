clc
url = 'https://fortlewiscollege-my.sharepoint.com/personal/gnomelab_fortlewis_edu/_layouts/15/onedrive.aspx?id=%2Fpersonal%2Fgnomelab_fortlewis_edu%2FDocuments%2FRaw%20Data%2FEcomapper%2FLake_Nighthorse%2FLog%20files';
%url = 'http://sensorbuoy.appspot.com/books/';
data = webread(url)
%data = urlread(url);

disp(data)
% hostname = char(getHostName(java.net.InetAddress.getLocalHost))