% Load C library

% Library is in : /usr/local/include/isam/isam.h
close all
clear all
clc

if not(libisloaded('libtiff'))
    addpath(genpath('/home/jake/Libraries/'))
%    addpath(genpath('/usr/include/c++/5'))
%    addpath(genpath('/usr/local/include'))
    loadlibrary('libtiff','tiffio.h')                % --> look for tiffio
    loadlibrary('libbsb','bsb.h')
end

libfunctions('libtiff')


% if not(libisloaded('libbsb'))
%     addpath(genpath('/home/jake/Libraries/'))
%     loadlibrary('libbsb','bsb.h')
% end


% sc = libstruct('c_struct')
% 
% set(sc,'p1',100,'p2',150,'p3',200)
% 
% get(sc)