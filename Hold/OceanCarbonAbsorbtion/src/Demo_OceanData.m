% Ocean Color Data From NASA =====================================================================================================
close all
clear all
clc


%% Build URL for Ocean Data -----------------------------------------------------------
url = 'https://oceandata.sci.gsfc.nasa.gov:443/opendap/MODISA/L3SMI';

date_start = datetime('01/01/2017', 'InputFormat', 'MM/dd/yyyy');
date_end   = datetime('02/01/2017', 'InputFormat', 'MM/dd/yyyy'); %datetime(date);

years = year(date_start): year(date_end);

indata = date_start: date_end;
ind = day(indata)==1;
firstDay = indata(ind);
firstDay = day(firstDay,'dayofyear');

ind(1) = [];
lastDay = indata(ind);
lastDay = day(lastDay,'dayofyear');


dates(size(years,2)*size(lastDay,2),:) = sprintf('%d/%03d/A%d%03d%d%03d', years(end), firstDay(end), years(end), firstDay(end), years(end), lastDay(end) );;

ind = 1;
for ii = 1: size(years,2)
for jj = 1: size(lastDay,2)
    dates(ind,:) = sprintf('%d/%03d/A%d%03d%d%03d', years(ii), firstDay(jj), years(ii), firstDay(jj), years(ii), lastDay(jj) );
    ind = ind + 1;

end
end


types = {'L3m_MO_GSM_chl_gsm_4km.nc';
         'L3m_MO_GSM_bbp_443_gsm_4km.nc';
         'L3m_MO_KD490_Kd_490_4km.nc';
         'L3m_MO_PAR_par_4km.nc'};


openDap_url{size(dates,1)*size(types,1)} = '';
     
ind = 1;
for ii = 1: size(dates,1)
for jj = 1: size(types,1)
        openDap_url{ind} = strcat(url, '/', dates(ii,:), '.', types{jj});
        ind = ind + 1;
end
end
        
 
%% Make the Call ====================================================================================================================
ocd = OceanColorData;

for ii = 1 : size(openDap_url,2)
ocd = ocd.AddDataSet(openDap_url{ii}, true);
end


% clear openDap_url


%% Plot Stuff
ocd.Plot('par',2)
