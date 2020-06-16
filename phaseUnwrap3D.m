function delPhaseNet = phaseUnwrap3D(vol_compl)
%  delPhaseNet = phaseUnwrap3D(vol_compl, num_echoes)
%  Temporally unwrap phase data
%
%  Input: vol_compl: complex dataset (mag.*exp(1i*phase))
%  Output: delPhaseNet: temporally unwrapped phase data
%
% Written by Avery Berman (ajberman@mgh.harvard.edu)
% Modified by Yuhan Ma (yuhanma1@gmail.com)
% Modifications: stopped setting the first echo phase data to zero
% Last modified: September 2013

%dim_mat = size(vol_compl(:,:,:,1));    % get just the image dimensions (w/o the time dimension)
num_echoes = size(vol_compl,4);

%delPhaseNet = zeros([dim_mat num_echoes]);
% delPhaseNet = zeros(size(vol_compl));

%delPhaseNet(:,:,:,1) = angle(squeeze(vol_compl(:,:,:,1)));
delPhaseNet(:,:,:,1) = angle(vol_compl(:,:,:,1));

for n=2:num_echoes
    
    vol_tmp0 = vol_compl(:,:,:,n-1);
    vol_tmp1 = vol_compl(:,:,:,n);

%               ^ Re
%               |
%               |  ^
%               | /
%               |/theta
%               --------> Im                       
    
    delPhase = angle( vol_tmp1.*conj(vol_tmp0) );
    delPhaseNet(:,:,:,n) = delPhaseNet(:,:,:,n-1) + delPhase;

end

end