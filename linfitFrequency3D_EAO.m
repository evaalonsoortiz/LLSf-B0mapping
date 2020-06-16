function [delf, offset, STDX, MSE] = linfitFrequency3D_EAO(delt, delPhaseNet, num_echoes, mag, mask)
%LINFITFREQUENCY3D linearly fits the change in phase vs. echo time
%   Detailed explanation goes here
% Input: delt: echo time 
%        delPhaseNet: temporarily unwrapped phase data
%        num_echoes: number of echoes used in the fitting
%        mag: magnitude dataset
% Ouptput: delf: frequency shift map (field map)'
%          offset: frequency shift at echo time 0
%          stdx: standard deviation of the fit
%          mse: mean square error of the fit
% Written by Avery Berman (ajberman@mgh.harvard.edu)
% Modified by Yuhan Ma (yuhanma1@gmail.com)
% Modification: Use my linearFit instead of Avery's linFit. Accept a 3D
% volume as an input
% Modification: use matlab built-in function lscov
% Last Modified July 2013


if max(reshape(mask,size(mask,1)*size(mask,2),size(mask,3))) > 1
    mask(find(mask > 0)) = 1;
end
 
delf = zeros(size(delPhaseNet,1),size(delPhaseNet,2),size(delPhaseNet,3));
offset = zeros(size(delPhaseNet,1),size(delPhaseNet,2),size(delPhaseNet,3));
STDX = zeros(size(delPhaseNet,1), size(delPhaseNet,2), size(delPhaseNet,3), 2);
MSE = zeros(size(delPhaseNet,1),size(delPhaseNet,2),size(delPhaseNet,3));

for i=1:size(delPhaseNet,1)
    for j=1:size(delPhaseNet,2)
        for k = 1:size(delPhaseNet,3)
            if (mask(i,j,k) == 0)
                delf(i,j,k) = 0;
            else
                A_design =  [delt(1:num_echoes)',ones(size(delt(1:num_echoes)'))];
                b = squeeze(delPhaseNet(i,j,k,1:num_echoes));
                weight = squeeze(mag(i,j,k,1:num_echoes)).^2;
                [fit_para_tmp, STDX(i,j,k,:), MSE(i,j,k)] = lscov(A_design,b, weight);
                delf(i,j,k) = fit_para_tmp(1);
                offset(i,j,k) = fit_para_tmp(2);
            end
        end
    end
end

end