% ------------------------------
% Create a mask
% ------------------------------
function mask = create_mask(mat,data_threshold,mask_radius)

if (~exist('mask_radius','var'))
	mask_radius = 0;
end;

mat = abs(mat);
mat = mat / max(mat(:));

mask = zeros(size(mat));
if exist('mask_radius','var') && (mask_radius ~= 0)
	Nx = size(mask,2);
	Ny = size(mask,1);
	x_mat = repmat( linspace(-Nx/2,+Nx/2,Nx)/Nx  ,Ny,1);
	y_mat = repmat((linspace(-Ny/2,+Ny/2,Ny)/Ny)',1,Nx);
	radius_mat = sqrt(x_mat.^2 + y_mat.^2);
	mask(radius_mat < mask_radius) = 1;
else
	mask(mat>data_threshold) = 1;
end;

return;
