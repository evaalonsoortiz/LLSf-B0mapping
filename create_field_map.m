function create_field_map(mag_fname, phase_fname, mask_fname, rescale_flag, temp_unwrap_flag)

% *************************************************************************
% function create_field_map(mag_fname, phase_fname)
%
% DESCRIPTION: Function to create a B0 map, based on phase images obtained
%              with an mGRE sequence. Field map is calculated by linear 
%              fitting of the phase data over time. 
%              
%
% INPUTS: 
% mag_fname : magnitude volume (4D vol, from a mGRE acquisition)
% phase_fname : phase volume (4D vol, from a mGRE acquisition)
% mask_fname : masked used to avoid linearly fitting noise voxels
%
% OUTPUT:
% field_map.mnc : field map in unit Hz
% 
% AUTHOR: Eva Alonso Ortiz
% 
%*************************************************************************

% constants
gamma = 2*pi*42.58e6; % rad*Hz/T

% open data files
[phase_desc,phase] = niak_read_minc(phase_fname);
[mag_desc,mag] = niak_read_minc(mag_fname);
[mask_desc,mask] = niak_read_minc(mask_fname);

% check that mask and data_vol are the same dimensions
if size(mask,1)*size(mask,2) ~= size(mag,1)*size(mag,2)
    error(sprintf('\nError: Mask file dimensions do not match data image file.\n')); 
end

if max(reshape(mask,size(mask,1)*size(mask,2),size(mask,3))) > 1
    mask(find(mask > 0)) = 1;
end

% find echo times
num_echoes = mag_desc.info.dimensions(1,4);
echo_times = calc_echo_times(num_echoes);
echo_times = echo_times*1e-3; % convert to seconds

% even_echoes_flag = input('Do you want to use only even echoes? (y=1/n=0): ');
% if even_echoes_flag == 1
%     echo_times = echo_times(2:2:num_echoes);
% else
%     odd_echos_flag = input('Do you want to use only odd echoes? (y=1/n=0): ');
%     if odd_echos_flag == 1
%         echo_times = echo_times(1:2:num_echoes);
%     end
% end
% 
% num_echoes = input('How many echoes do you want to use?: ');

%% rescale phase images
if rescale_flag == 0
    rescaled_phase(:,:,:,:) = phase(:,:,:,:);
elseif rescale_flag == 1
    rescaled_phase(:,:,:,:) =  rescalePhaseImage(phase(:,:,:,:));
end

% create complex data volume
vol_compl(:,:,:,:) = mag(:,:,:,:).*exp(1i*rescaled_phase(:,:,:,:));

% temporal phase unwrapping
if temp_unwrap_flag == 0
    delPhaseNet(:,:,:,:) = rescaled_phase(:,:,:,:);
elseif temp_unwrap_flag == 1
    delPhaseNet = phaseUnwrap3D(vol_compl);
end

% frequency shift calculation (field map)
[delf, offset, STDX, MSE] = linfitFrequency3D_EAO(echo_times, delPhaseNet, num_echoes, mag, mask);
B0 = delf/gamma; % convert from rad*Hz to T

%% save rescaled images
[pathstr, name, ext] = fileparts(phase_fname);
% fname_rescaled_image = strcat(name,'_rescaled.mnc');
% phase_desc.file_name = fname_rescaled_image;
% niak_write_minc(phase_desc,rescaled_phase);

% save temporally unwrapped phase image
fname_temp_unwrapped_image = strcat(name,'_temp_unwarpped.mnc');;
phase_desc.file_name = fname_temp_unwrapped_image;
niak_write_minc(phase_desc,delPhaseNet);

% save field map
fname_field_map = strcat(name,'_field_map.mnc');
phase_desc.file_name = fname_field_map;
phase_desc.info.dimensions = [phase_desc.info.dimensions(1,1), phase_desc.info.dimensions(1,2), phase_desc.info.dimensions(1,3), 1];
niak_write_minc(phase_desc,B0);

% save offset map
% fname_offset_map = strcat(name,'_offset_map.mnc');
% phase_desc.file_name = fname_offset_map;
% phase_desc.info.dimensions = [phase_desc.info.dimensions(1,1), phase_desc.info.dimensions(1,2), phase_desc.info.dimensions(1,3), 1];
% niak_write_minc(phase_desc,offset);
