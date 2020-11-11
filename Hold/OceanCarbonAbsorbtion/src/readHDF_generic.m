function result=readHDF_generic(filename, prodName)
% function [lat,lon,result]=readHDF_generic(filename)

%filename=input('which file?? ','s');
SD_id = hdfsd('start',filename,'read');
%[ndatasets,nglobal_attr,status] = hdfsd('fileinfo',SD_id);
[numdata,numdescr]=hdfsd('fileinfo',SD_id);
sds_id=hdfsd('select',SD_id,0);
[name,numdim,dimvector,type,numdescr]=hdfsd('getinfo',sds_id);
% for the 4km daily scenes, name1='Mapped - Composited'
% for the 4km 8-day composites, name1='Mapped-Composite'   - an old typo 
if nargin==3
   name1=prodName;       %%% -- product name needs to be known a priori
else
   name1=name;
end
%name1=prodName;       %%% -- product name needs to be known a priori
%%name1='DataSet';         chl from robert should be stored by this name
%%name1='mld';             this should be name for MLD data

id_u = hdfsd('nametoindex',SD_id,name1);
sds_id1 = hdfsd('select',SD_id, id_u);
[name,rank,dimsizes,data_type,nattrs,status] = hdfsd('getinfo',sds_id1);

% read the data data array
[data,status] = hdfsd('readdata',sds_id1,[0,0],[1,1],dimsizes);

% close access to sds array and to HDF file
status = hdfsd('endaccess',sds_id1);
status = hdfsd('end',SD_id);

% Now scale and plot the data

data=double(data);  % convert int16 array into a standard Matlab double format
data=data';           % reverse row column format (hdf convention)
%data(data==255)=NaN;

%scale to geophysical chl units

result=data;
