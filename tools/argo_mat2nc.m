% This script loads the Argo profile data in .mat format
% (WOD13_PFL_2000.mat) and converts it to netCDF format
%
% Qing Li, 20180322

close all; clear variables;

% data directory
datadir = '/Volumes/Qing_Work/data/obs/Argo';
% input data
filename = 'WOD13_PFL_2000.mat';
indata = [datadir '/' filename];
% output data
[~, fname, ~] = fileparts(filename);
outdata = [datadir '/' fname '.nc'];

% load data
dat = load(indata);

% process time
refdate_str = '1990-01-01 00:00:00';
time_unsorted = dat.DateNum - datenum(refdate_str);

% sort time dimension
[time, idx_sort] = sort(time_unsorted);

% output variable name and meta data
nnt = numel(dat.DateNum);
nnz = numel(dat.Depth);

% one dimensional variables
var1dListOut = {'Basin', 'FloatNum', 'Latitude', 'Longitude',...
                'MaxDepth', 'SalnCastFlag'};          
var1dLongName = {'Basin index number of the Argo cast location',...
                 'Argo float number',...
                 'Latitude',...
                 'Longitude',...
                 'Longitude',...
                 'Maximum depth this Argo cast',...
                 'Quality Control flag of entire salinity profile '};
var1dUnits = {'unitless', 'unitless', 'degree_north', 'degree_east',...
              'm', 'unitless'};
var1dType = {'NC_INT', 'NC_INT', 'NC_DOUBLE', 'NC_DOUBLE', 'NC_DOUBLE', 'NC_INT'};

% two dimensional variables
var2dListOut = {'Salinity', 'Temperature', 'SalnFlag', 'TempFlag'};
var2dLongName = {'Salinity',...
                 'Temperature',...
                 'Quality Control flag for salinity',...
                 'Quality Control flag for temperature'};
var2dUnits = {'psu', 'degC', 'unitless', 'unitless'};
var2dType = {'NC_DOUBLE', 'NC_DOUBLE', 'NC_INT', 'NC_INT'};

nvar1d = numel(var1dListOut);
nvar2d = numel(var2dListOut);         

date_str = datestr(now,'yyyymmdd');

% fill values
fv_int = netcdf.getConstant('NC_FILL_INT');
fv_double = netcdf.getConstant('NC_FILL_DOUBLE');

% create netCDF
ncid = netcdf.create(outdata,'NETCDF4');
% Define dimensions
castDimId = netcdf.defDim(ncid, 'cast', nnt);
zDimId = netcdf.defDim(ncid, 'depth', nnz);

% Define variables
varid_z  = netcdf.defVar(ncid, 'depth', 'NC_INT', zDimId);
netcdf.putAtt(ncid, varid_z, 'long_name', 'depth');
netcdf.putAtt(ncid, varid_z, 'units', 'm');

varid_time  = netcdf.defVar(ncid, 'time', 'NC_DOUBLE', castDimId);
netcdf.putAtt(ncid, varid_time, 'long_name', 'time');
netcdf.putAtt(ncid, varid_time, 'units', ['days since ' refdate_str]);

varid_cast  = netcdf.defVar(ncid, 'cast', 'NC_INT', castDimId);
netcdf.putAtt(ncid, varid_cast, 'long_name', 'Cast number');
netcdf.putAtt(ncid, varid_cast, 'units', 'unitless');

varid1d = zeros(nvar1d,1);
varid2d = zeros(nvar2d,1);
for i = 1:nvar1d
    if strcmp(var1dType{i}, 'NC_INT')
        fillvalue = fv_int;
    elseif strcmp(var1dType{i}, 'NC_DOUBLE')
        fillvalue = fv_double;
    end
    varid1d(i) = netcdf.defVar(ncid, var1dListOut{i},...
        var1dType{i}, castDimId);
    netcdf.defVarFill(ncid, varid1d(i), false, fillvalue)
    netcdf.putAtt(ncid, varid1d(i), 'long_name', var1dLongName{i});
    netcdf.putAtt(ncid, varid1d(i), 'units', var1dUnits{i});
end
for i = 1:nvar2d
    if strcmp(var2dType{i}, 'NC_INT')
        fillvalue = fv_int;
    elseif strcmp(var2dType{i}, 'NC_DOUBLE')
        fillvalue = fv_double;
    end
    varid2d(i) = netcdf.defVar(ncid, var2dListOut{i},...
        var2dType{i}, [zDimId, castDimId]);
    netcdf.defVarFill(ncid, varid2d(i), false, fillvalue)
    netcdf.putAtt(ncid, varid2d(i), 'long_name', var2dLongName{i});
    netcdf.putAtt(ncid, varid2d(i), 'units', var2dUnits{i});
end
% global attributes
varid_global = netcdf.getConstant('NC_GLOBAL');
title = 'Argo profile data from 2004 to 2013 (WOD13_PFL_2000)';
netcdf.putAtt(ncid, varid_global, 'title', title);
netcdf.putAtt(ncid, varid_global, 'date', char(datetime));
history = 'Converted from WOD13_PFL_2000.mat';
netcdf.putAtt(ncid, varid_global, 'history', history);
netcdf.putAtt(ncid, varid_global, 'script', 'argo_mat2nc.m');
netcdf.putAtt(ncid, varid_global, '_FillValue', fv_double);

% end define mode
netcdf.endDef(ncid);

% write to netCDF file
netcdf.putVar(ncid, varid_z, int32(dat.Depth));
netcdf.putVar(ncid, varid_cast, 1:nnt);
netcdf.putVar(ncid, varid_time, time);
for i = 1:nvar1d
    vardata = squeeze(dat.(var1dListOut{i})(idx_sort));
    netcdf.putVar(ncid, varid1d(i), vardata);
end
for i = 1:nvar2d
    vardata = squeeze(dat.(var2dListOut{i})(idx_sort,:));
    netcdf.putVar(ncid, varid2d(i), vardata');
end

% close netCDF file
netcdf.close(ncid)
