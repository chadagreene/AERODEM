% This script converts Korsgaard's G150 AERODEM to Bedmachine's grid. 
% Download all the files from here: https://www.nodc.noaa.gov/archive/arc0088/0145405/1.1/data/0-data/G150AERODEM/
%
% This script takes a while to run! Give it 40 minutes on my laptop, give
% or take. 
%
% Chad A. Greene, NASA/JPL, December 2022. 

%% Load bedmachine data 

[geoid,x,y] = bedmachine_data('geoid','greenland'); 
sfz = bedmachine_data('surface','greenland'); 
th = bedmachine_data('thickness','greenland'); 
bed = bedmachine_data('bed','greenland'); 
mask = bedmachine_data('mask','greenland'); 

[X,Y] = meshgrid(x,y); 

[Lat,Lon] = psn2ll(X,Y); 

utmz = 19:27; % utm zones
X_utmz = nan(size(X,1),size(X,2),length(utmz)); 
Y_utmz = X_utmz; 

for k = 1:length(utmz)
    proj = projcrs(32600+utmz(k)); 
    [X_utmz(:,:,k),Y_utmz(:,:,k)] = projfwd(proj,Lat,Lon); 
    k
end

%%

yrz = [1978 1981 1985 1987]; 
sfz_korsgaard = NaN(size(X,1),size(X,2),length(yrz)); 

for kyr = 1:length(yrz)

    % This year's filenames: 
    d = dir(['/Users/cgreene/Documents/data/DEMs/korsgaard_g150_aerodem/aerodem_',num2str(yrz(kyr)),'_*.tif']);
    tmp_sfz = nan(size(X)); 

    % For each region (some utm zones occur more than once) 
    for k = 1:length(d)
    
        [tmp,xtmp,ytmp] = geoimread(d(k).name);
        tmp(tmp<-1000) = nan; 
    
        switch d(k).name(14:18)
            case 'utm19'
                Xtmp = X_utmz(:,:,1); 
                Ytmp = Y_utmz(:,:,1); 
            case 'utm20'
                Xtmp = X_utmz(:,:,2); 
                Ytmp = Y_utmz(:,:,2); 
            case 'utm21'
                Xtmp = X_utmz(:,:,3); 
                Ytmp = Y_utmz(:,:,3); 
            case 'utm22'
                Xtmp = X_utmz(:,:,4); 
                Ytmp = Y_utmz(:,:,4); 
            case 'utm23'
                Xtmp = X_utmz(:,:,5); 
                Ytmp = Y_utmz(:,:,5); 
            case 'utm24'
                Xtmp = X_utmz(:,:,6); 
                Ytmp = Y_utmz(:,:,6); 
            case 'utm25'
                Xtmp = X_utmz(:,:,7); 
                Ytmp = Y_utmz(:,:,7);  
            case 'utm26'
                Xtmp = X_utmz(:,:,8); 
                Ytmp = Y_utmz(:,:,8); 
            case 'utm27'
                Xtmp = X_utmz(:,:,9); 
                Ytmp = Y_utmz(:,:,9); 
            otherwise
                error('unrecognized utm zone') 
        end

        isn = isnan(tmp_sfz); 
        tmp_sfz(isn) = interp2(xtmp,ytmp,tmp,Xtmp(isn),Ytmp(isn)); 

        tmp_sfz((tmp_sfz-geoid)<0 | abs((sfz_korsgaard(:,:,1)- geoid) - sfz)>1000) = NaN; 

    k
    end
    sfz_korsgaard(:,:,kyr) = tmp_sfz; 
    
    disp(['Year ',num2str(kyr)])
end

%% Reliability mask

reliability = zeros(size(X,1),size(X,2),length(yrz),'uint8'); 

for kyr = 1:length(yrz)

    % This year's filenames: 
    d = dir(['/Users/cgreene/Documents/data/DEMs/korsgaard_g150_aerodem/rm_',num2str(yrz(kyr)),'_*.tif']);
    tmp_rm = reliability(:,:,kyr); 

    % For each region (some utm zones occur more than once) 
    for k = 1:length(d)
    
        [tmp,xtmp,ytmp] = geoimread(d(k).name);
    
        switch d(k).name(9:13)
            case 'utm19'
                Xtmp = X_utmz(:,:,1); 
                Ytmp = Y_utmz(:,:,1); 
            case 'utm20'
                Xtmp = X_utmz(:,:,2); 
                Ytmp = Y_utmz(:,:,2); 
            case 'utm21'
                Xtmp = X_utmz(:,:,3); 
                Ytmp = Y_utmz(:,:,3); 
            case 'utm22'
                Xtmp = X_utmz(:,:,4); 
                Ytmp = Y_utmz(:,:,4); 
            case 'utm23'
                Xtmp = X_utmz(:,:,5); 
                Ytmp = Y_utmz(:,:,5); 
            case 'utm24'
                Xtmp = X_utmz(:,:,6); 
                Ytmp = Y_utmz(:,:,6); 
            case 'utm25'
                Xtmp = X_utmz(:,:,7); 
                Ytmp = Y_utmz(:,:,7);  
            case 'utm26'
                Xtmp = X_utmz(:,:,8); 
                Ytmp = Y_utmz(:,:,8); 
            case 'utm27'
                Xtmp = X_utmz(:,:,9); 
                Ytmp = Y_utmz(:,:,9); 
            otherwise
                error('unrecognized utm zone') 
        end

        isn = tmp_rm==0; 
        tmp_rm(isn) = interp2(xtmp,ytmp,tmp,Xtmp(isn),Ytmp(isn),'nearest',0); 

    k
    end
    reliability(:,:,kyr) = tmp_rm; 
    
    disp(['Year ',num2str(kyr)])
end

%% Combine years and convert to geoid reference

year_mask = uint8(zeros(size(X))); 
reliability_mask = year_mask; 
sfz = NaN(size(X),'single'); 

for k = 1:4
   tmp_sfz = sfz_korsgaard(:,:,k) - geoid;
   tmp_rel = reliability(:,:,k); 

   isf = isfinite(tmp_sfz) & isnan(sfz) ;
   sfz(isf) = tmp_sfz(isf); 
   reliability_mask(isf) = tmp_rel(isf); 
   
   year_mask(isf) = yrz(k)-1900; 
end

%% Write the netcdf 

nanz = isnan(sfz); 
sfz_16 = int16(sfz); 
sfz_16(nanz) = -32767;

proj4 = '+proj=stere +lat_0=90 +lat_ts=70 +lon_0=-45 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +type=crs';

% 1. Create netCDF file handle and Attributes
mode = netcdf.getConstant('NETCDF4');
mode = bitor(mode,netcdf.getConstant('CLASSIC_MODEL'));
ncid=netcdf.create('AERODEM_Korsgaard_mosaic.nc',mode);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Conventions','CF-1.7');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Title','AERODEM');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Description','AERODEM digital elevation model of the Greenland Ice Sheet for 1978-1987, compiled, referenced to the geoid, and placed on the BedMachine grid for convenience.');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Author','Korsgaard, N., Nuth, C., Khan, S. et al.');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'creation_date',datestr(now,'yyyy-mm-dd'));
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'NetCDF_conversion','Chad A. Greene');
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'tmd_version',3.0);
netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'),'Data_citation','Korsgaard, N., Nuth, C., Khan, S. et al. Digital elevation model and orthophotographs of Greenland based on aerial photographs from 1978–1987. Sci Data 3, 160032 (2016). https://doi.org/10.1038/sdata.2016.32')

% 2. Define dimensions
% Define mapping variable
mapping_var_id= netcdf.defVar(ncid,'mapping','NC_CHAR',[]);
netcdf.putAtt(ncid,mapping_var_id,'grid_mapping_name', 'polar_stereographic');
netcdf.putAtt(ncid,mapping_var_id,'latitude_of_projection_origin',90);
netcdf.putAtt(ncid,mapping_var_id,'standard_parallel',70);
netcdf.putAtt(ncid,mapping_var_id,'straight_vertical_longitude_from_pole',-45);
netcdf.putAtt(ncid,mapping_var_id,'false_easting',0);
netcdf.putAtt(ncid,mapping_var_id,'false_northing',0); 
netcdf.putAtt(ncid,mapping_var_id,'spatial_proj4',proj4); 
netcdf.putAtt(ncid,mapping_var_id,'spatial_epsg',3413);
netcdf.putAtt(ncid,mapping_var_id,'datum','eigen-6c4 geoid (Forste et al 2014)'); 

% Define x: 
x_id     = netcdf.defDim(ncid,'x',length(x));
x_var_id = netcdf.defVar(ncid,'x','NC_FLOAT',x_id);
netcdf.putAtt(ncid,x_var_id,'standard_name','projection_x_coordinate');
netcdf.putAtt(ncid,x_var_id,'long_name',    'Cartesian x-coordinate, grid cell center');
netcdf.putAtt(ncid,x_var_id,'units',        'meters');

% Define y: 
y_id     = netcdf.defDim(ncid,'y',length(y));
y_var_id = netcdf.defVar(ncid,'y','NC_FLOAT',y_id);
netcdf.putAtt(ncid,y_var_id,'standard_name','projection_y_coordinate');
netcdf.putAtt(ncid,y_var_id,'long_name',    'Cartesian y-coordinate, grid cell center');
netcdf.putAtt(ncid,y_var_id,'units',        'meters');

% Define surface
sfz_var_id = netcdf.defVar(ncid,'surface','NC_SHORT',[x_id y_id]);
netcdf.putAtt(ncid,sfz_var_id,'standard_name','surface height');
netcdf.putAtt(ncid,sfz_var_id,'long_name',    'surface height relative to the EIGEN-6c4 geoid (Same as BedMachine v5).');
netcdf.putAtt(ncid,sfz_var_id,'grid_mapping', 'polar_stereographic');
netcdf.putAtt(ncid,sfz_var_id,'units',        'm');
netcdf.putAtt(ncid,sfz_var_id,'_FillValue',   int16(-32767));

% Define year_mask
yr_var_id = netcdf.defVar(ncid,'year_mask','NC_BYTE',[x_id y_id]);
netcdf.putAtt(ncid,yr_var_id,'standard_name','year mask');
netcdf.putAtt(ncid,yr_var_id,'long_name',    'Year of data collection.');
netcdf.putAtt(ncid,yr_var_id,'flag_values', [0 78 81 85 87]);
netcdf.putAtt(ncid,yr_var_id,'flag_meanings', 'no data, 1978, 1981, 1985, 1987');
netcdf.putAtt(ncid,yr_var_id,'grid_mapping', 'polar_stereographic');

% Define reliability_mask
rel_var_id = netcdf.defVar(ncid,'reliability_mask','NC_BYTE',[x_id y_id]);
netcdf.putAtt(ncid,rel_var_id,'standard_name','reliability mask');
netcdf.putAtt(ncid,rel_var_id,'long_name',    'Reliability estimates.');
netcdf.putAtt(ncid,rel_var_id,'flag_values',[0 100]);
netcdf.putAtt(ncid,rel_var_id,'flag_meanings', 'See the manuscript, but generally values >=40 are considered more reliable.');
netcdf.putAtt(ncid,rel_var_id,'grid_mapping', 'polar_stereographic');

% Compress and stop variable definition
netcdf.defVarDeflate(ncid,sfz_var_id,true,true,9);
netcdf.defVarDeflate(ncid,yr_var_id,true,true,9);
netcdf.defVarDeflate(ncid,rel_var_id,true,true,9);
netcdf.endDef(ncid);

%3. Place data
netcdf.putVar(ncid,x_var_id,single(x));
netcdf.putVar(ncid,y_var_id,single(y));
netcdf.putVar(ncid,sfz_var_id,ipermute(sfz_16,[2 1]));
netcdf.putVar(ncid,yr_var_id,ipermute(year_mask,[2 1]));
netcdf.putVar(ncid,rel_var_id,ipermute(reliability_mask,[2 1]));

%4. Close file 
netcdf.close(ncid)

disp done



