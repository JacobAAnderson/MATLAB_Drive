% This function extracts all the data in a netCD file and places it in the base work station

function finfo = Extract_NetCF_File(file)
% Get Info on the file
finfo    = ncinfo(file);
varnames = {finfo.Variables(:).Name};

% fileData = finfo.Attributes(16).Value;
% fileData = strsplit(fileData,',');
% 
% disp('______________________________________________')
% disp('NetCD file Contents:')
% 
% for ii = 1: length(fileData)
%     disp(fileData{ii})
% end

disp(' ')
disp('______________________________________________')
disp('Geting Data From:')
disp(file)
disp(' ')
disp('Imported varables: ')


% Get Data out of the file
ncid = netcdf.open(file);
% varids = netcdf.inqVarIDs(ncid);

for ii = 1: length(varnames)
    
    varid = netcdf.inqVarID(ncid,varnames{ii});    
    data  = netcdf.getVar(ncid,varid);
    assignin('base',varnames{ii}, data);
    
    disp(varnames{ii})
end

netcdf.close(ncid)


disp('______________________________________________')
disp(' ')

end
