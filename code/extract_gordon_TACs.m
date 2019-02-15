function gordonROI = extract_gordon_TACs(EPI,parcels,parcels_list)

% Extract time activity curves of ROIs identified by Gordon parcellation
% 201902130 ES (Reduced Version)

% INPUT:
% EPI                        fMRI data fullpath
% parcels                    Segmentation 
% parcels_list               Segmentation ROIs infos

% OUTPUT:
% gordonROI.ID               Region ID in Gordon Parcellation
% gordonROI.Hemi             Hemisphere
% gordonROI.Net              Network
% gordonROI.idx              Index of voxels in the EPI volume
% gordonROI.TAC              Time activity curve
% gordonROI.nCompExpl80      Number of principal components required to
%                            describe the 80% of the ROI signals variance

explained_thr = 80;

%% Reshape

nVols    = size(EPI,4);
EPI_2D        = reshape(EPI, [], nVols)';
Gordon1D      = reshape(parcels,1,[]);

%% Load Gordon parcels infos
[T, ~, ~, ~]    = make_Gordon_parcels_table(parcels_list);

%% Parcels infos and time course extraction:

disp('Extracting time courses of the parcels')

for jj = 1: length(T.ID_sorted)
    
    disp(['Working on parcel #' num2str(T.ID_sorted(jj))])
    
    gordonROI(jj).ID               = T.ID_sorted(jj);
    gordonROI(jj).Hemi             = T.HEM{jj};
    gordonROI(jj).Net              = T.NETID(jj);
    gordonROI(jj).idx              = find(Gordon1D == T.ID_sorted(jj));
    gordonROI(jj).epiROIdim        = length(gordonROI(jj).idx);

    tc                             = EPI_2D(:,gordonROI(jj).idx');
    gordonROI(jj).TAC              = mean(tc,2);
    
    if not(isempty(gordonROI(jj).idx))
        [~, score, ~, ~, explained] = pca(tc-mean(tc,1));
        nc = find(cumsum(explained)>explained_thr);
        eval(['gordonROI(jj).nCompExpl' num2str(explained_thr) ' = nc(1);'])
        gordonROI(jj).corr1PCA = corr(gordonROI(jj).TAC,score(:,1));
    end
    
    clear tc
    
end
