
function [D_compressed] = nbe_calc_PCA(D,nPCs,plot_f)

N        = size(D,1);
mean_D   = mean(D,1);
mean_sub = repmat(mean_D,N,1);
U        = D - mean_sub;

% Calculate the covariance matrix
covD = (transpose(U)*U) / (N-1);	  % Identical to MATLAB's  cov(U)

% Find eigenvalues, and eigenvectors (PCs)
[EMC_PC,covD_eigVals] = eig(covD);
covD_eigVals_arr = diag(covD_eigVals);

% Sort PCs according to magnitude
[~,eig_idcs] = sort(abs(covD_eigVals_arr),1,'descend');
EMC_PC = EMC_PC(:,eig_idcs);
EMC_PC_reduced = EMC_PC(:,1:nPCs);

% Calculate a compressed EMC mat
D_compressed = U * EMC_PC_reduced;

if (exist('plot_f','var')) && plot_f
	figure;
	subplot(121); imagesc(D);            title('Original   mat'); colorbar;
	subplot(122); imagesc(D_compressed); title('Compressed mat'); colorbar;
end;

return;

