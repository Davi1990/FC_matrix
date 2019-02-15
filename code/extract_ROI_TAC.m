function TAC = extract_ROI_TAC(EPI,parcels,ROI_id)

% Extract time activity curves of ROIs identified by Gordon parcellation
% 201902130 ES (Reduced Version)

% INPUT:
% EPI               fMRI 4D data
% parcels           Segmentation 3D data
% parcels_list      ROI IDs list

% OUTPUT:
% TAC               Region ID in Gordon Parcellation


nVols     = size(EPI,4);
EPI_2D    = reshape(EPI, [], nVols)';
Gordon1D  = reshape(parcels,1,[]);

tc        = [];
for jj = 1: length(ROI_id)
    tc      = [tc EPI_2D(:,find(Gordon1D == ROI_id(jj)))];
end
TAC       = mean(tc,2);

