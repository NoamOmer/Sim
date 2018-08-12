
function out_vol = freesurfer_erode_seg(in_vol,segs_info,nPx,erode_th,erode_level)

NON_EXIST_SEG_VAL = -111;

if (sum(in_vol(:))==0)
	return;
end;

% Create a list of segments that needs to be eroded
seg_vals_to_erode = [];
seg_rngs_to_erode = {};
for idx = 1:size(segs_info,1)
	switch segs_info{idx,2}
	case 'val'
		seg_vals_to_erode = [seg_vals_to_erode segs_info{idx,3}];
	case 'rng'
		seg_rngs_to_erode{end+1} = segs_info{idx,3};
	otherwise
		error('...');
	end;

end;
seg_vals_to_erode = unique(seg_vals_to_erode);
	

% Initialize an output volume
out_vol = zeros(size(in_vol));

% ------------------------------------
% Erode specific segments values
% ------------------------------------
for seg_val = seg_vals_to_erode
	tmp_vol = zeros(size(in_vol));

	% Set all currently calculated segment locations to a non-existant segment value
	tmp_vol(in_vol==seg_val) = NON_EXIST_SEG_VAL;
	% create the mask
	tmp_vol(tmp_vol~=NON_EXIST_SEG_VAL) = 0;
	tmp_vol(tmp_vol==NON_EXIST_SEG_VAL) = 1;

	cur_n_px = sum(tmp_vol(:));
	if (cur_n_px==0),  continue;  end;
	
	% Erode
	tmp_vol_orig = tmp_vol;
	for px_idx = 1:nPx
		tmp_vol1 = tmp_vol;  tmp_vol1(1:end-1,:,:) = tmp_vol1(2:end  ,:,:);
		tmp_vol2 = tmp_vol;  tmp_vol2(2:end  ,:,:) = tmp_vol2(1:end-1,:,:);
		tmp_vol3 = tmp_vol;  tmp_vol3(:,1:end-1,:) = tmp_vol3(:,2:end  ,:);
		tmp_vol4 = tmp_vol;  tmp_vol4(:,2:end  ,:) = tmp_vol4(:,1:end-1,:);
		tmp_vol5 = tmp_vol;  tmp_vol5(:,:,1:end-1) = tmp_vol5(:,:,2:end  );
		tmp_vol6 = tmp_vol;  tmp_vol6(:,:,2:end  ) = tmp_vol6(:,:,1:end-1);

		tmp_vol = tmp_vol1 + ...
		          tmp_vol2 + ...
		          tmp_vol3 + ...
		          tmp_vol4 + ...
		          tmp_vol5 + ...
		          tmp_vol6;
		tmp_vol(tmp_vol  < erode_level) = 0;
		tmp_vol(tmp_vol ~=           0) = 1;
	end;

	% Don't erode if erosion causes the number of pixels to reduce below a certain threhold
	new_n_px = sum(tmp_vol(:));
	if (new_n_px < erode_th*cur_n_px)
		out_vol = out_vol + seg_val*tmp_vol_orig;
	else
		out_vol = out_vol + seg_val*tmp_vol;
	end;
end;

% ------------------------------------
% Erode segments ranges
% ------------------------------------
for rng_idx = 1:length(seg_rngs_to_erode)
	tmp_vol = zeros(size(in_vol));
	cur_rng = seg_rngs_to_erode{rng_idx};
	rep_val = cur_rng(1);

	% Set all currently calculated segment locations to a non-existant segment value
	tmp_vol((in_vol>=cur_rng(1)) & ...
	        (in_vol<=cur_rng(2))) = NON_EXIST_SEG_VAL;  % important to use & and not && in this case

	% create the mask
	tmp_vol(tmp_vol~=NON_EXIST_SEG_VAL) = 0;
	tmp_vol(tmp_vol==NON_EXIST_SEG_VAL) = 1;

	cur_n_px = sum(tmp_vol(:));
	if (cur_n_px==0),  continue;  end;
	
	% Erode
	tmp_vol_orig = tmp_vol;
	for px_idx = 1:nPx
		tmp_vol1 = tmp_vol;  tmp_vol1(1:end-1,:,:) = tmp_vol1(2:end  ,:,:);
		tmp_vol2 = tmp_vol;  tmp_vol2(2:end  ,:,:) = tmp_vol2(1:end-1,:,:);
		tmp_vol3 = tmp_vol;  tmp_vol3(:,1:end-1,:) = tmp_vol3(:,2:end  ,:);
		tmp_vol4 = tmp_vol;  tmp_vol4(:,2:end  ,:) = tmp_vol4(:,1:end-1,:);
		tmp_vol5 = tmp_vol;  tmp_vol5(:,:,1:end-1) = tmp_vol5(:,:,2:end  );
		tmp_vol6 = tmp_vol;  tmp_vol6(:,:,2:end  ) = tmp_vol6(:,:,1:end-1);

		tmp_vol = tmp_vol1 + ...
		          tmp_vol2 + ...
		          tmp_vol3 + ...
		          tmp_vol4 + ...
		          tmp_vol5 + ...
		          tmp_vol6;
		tmp_vol(tmp_vol  < erode_level) = 0;
		tmp_vol(tmp_vol ~= 0          ) = 1;
	end;

	% Don't erode if erosion causes the number of pixels to reduce below a certain threhold
	new_n_px = sum(tmp_vol(:));
	if (new_n_px < erode_th*cur_n_px)
		out_vol = out_vol + rep_val*tmp_vol_orig;
	else
		out_vol = out_vol + rep_val*tmp_vol;
	end;
end;

return;


			
