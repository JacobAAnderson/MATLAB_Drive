% Phyto Server
% Jacob Anderson
% 01/29/2020


path = '/data1/oceancolor/MODISA/Mapped/R2018/Monthly_Climatology/par/A20021822017212.L3m_MC_PAR_par_9km.nc';
disp(path);

[filePath, name, ext] = fileparts(path);

fprintf('\nBegin Importing Data From: %s\n', filePath)

files = fullfile(filePath,['*',ext]);
datafiles = dir(fullfile(files));

fprintf('\n %d  %s Files Found\n',size(datafiles,1),ext);


filePath = cellstr(repmat(filePath, size(datafiles,1),1));                   % Collect file names
files = fullfile(filePath, {datafiles.name}');

disp(datafiles.name)


x = ncread(path, 'par');

whos('x')

disp('HaZaAAA')

imshow(x)

