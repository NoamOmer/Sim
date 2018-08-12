% Calculate 1D inhomogeneity polynomial
% Input:
%  map_mat   :  [Hz]   :  2D phase map (already converted from phase to freq.)
%  map_mask  :         :  2D mask map, marking the relevant pixels in the map
%  x_axis    :  [cm]   :  Map x-axis
%  y_axis    :  [cm]   :  Map y-axis
%  rot_phi   :  [deg]  :  Axis, along which the inhomogeneity polinomial need to be calculated.
%                         (rotation is counter-clockwise)
%  ax0       :         :  (optional) Handle to axes object, on which to plot the inhomogeneity fitting
%  context   :         :  (optional)
function [z_axis,P_DNu0,PolyOrder] = calc_oriented_1D_inhomo(map_mat,map_mask,x_axis,y_axis,rot_phi,ax0,ax1,PolyOrder)
set_globals;
% ----------------------------------------------------------------------------------------
% Rotate the map so that the wanted orientation will fall along the x-axis
%  ** It is important to use the 'crop' option in order for the spatial scaling
%     to stay the same (thus enabling the use of the un-rotated x_axis in the next stages)
%  ** 'imrotate' rotates the axes in a counter-clockwise direction where we
%      need a counter-clockwise rotation of the matix, hence the use of -rot_phi
% ----------------------------------------------------------------------------------------
% DEBUG:: % m = map_mat;% m(map_mask==0) = mean(mean(m(map_mask==1)));% figure;imagesc(m);% title('Phase map (pre-rotation)');
input_map_mat  = map_mat;
input_map_mask = map_mask;
crop_flag = 1;
if (crop_flag)
    map_mat  = imrotate(map_mat ,-rot_phi,'crop');
    map_mask = imrotate(map_mask,-rot_phi,'crop');
else
    map_mat  = imrotate(map_mat ,-rot_phi,'loose');
    map_mask = imrotate(map_mask,-rot_phi,'loose');
end;

% The map z_axis units should be a combination of x & y axes.
% The following method is therefore ok only for angles of 0 or 90!
% uiwait(msgbox('Check the assignment of z_axis (correct for 0 or 90 deg angles only!)'));
% REMOVED FOR COMFORT PURPOSES % uiwait(msgbox('Check the assignment of z_axis (correct for 0 or 90 deg angles only!)'));
map_z_axis = y_axis;

% DEBUG:: % m = map_mat;% m(map_mask==0) = mean(mean(m(map_mask==1)));% figure;imagesc(m); %set(gca,'ydir','normal');% title(sprintf('Phase map (rotation = rot phi = %f)',rot_phi));
% figure; 
% subplot(1,2,1); imagesc(map_mat);
% subplot(1,2,2); imagesc(map_mask);

% % The following line is redundant: in the next code we ignore any pixel where (map_mask == 1)
% map_mat(map_mask == 0) = mean(mean(map_mat(map_mask == 1))); 
% % The following line is redundant: already implicitly done in the previous line
% map_mat = map_mat - mean(mean(map_mat));

% -------------------------------------------------------------------------------------------
% (1) Reduce the map to a 1D vector by averaging the contribution of all y-pixels at each x-pixel
%     Do not use the data in the column if less than 10% of the pixels contribute.
% (2) Calculate the standard deviation in each averaged pixel of the 1D vector
% -------------------------------------------------------------------------------------------
for idx = 1:size(map_mat,1)
	ortho_vec = map_mat(idx,:);
	mask_vec  = map_mask(idx,:);
	if ((sum(mask_vec == 0) > 0.98*length(mask_vec)))
		map_1DAvgInhomo_vec(idx) = NaN;
        map_1DStdInhomo_vec(idx) = NaN;
	else
		map_1DAvgInhomo_vec(idx) = mean(ortho_vec(mask_vec == 1));
		map_1DStdInhomo_vec(idx) =  std(ortho_vec(mask_vec == 1));
        
%       % Remove points that are 1 standard deviations away from the mean
%       std_mask_vec = abs(ortho_vec - map_1DAvgInhomo_vec(idx)) > 1.0*map_1DStdInhomo_vec(idx);
%       std_mask_vec = ~std_mask_vec;
%         
% 		map_1DAvgInhomo_vec(idx) = mean(ortho_vec(mask_vec.*std_mask_vec == 1));
% 		map_1DStdInhomo_vec(idx) =  std(ortho_vec(mask_vec.*std_mask_vec == 1));
	end;
end;

% ------------------------------------------------------
% Create spatial and frequency axes in the relevant area
% ------------------------------------------------------
idxs = find(~isnan(map_1DAvgInhomo_vec));
map_1DAvgInhomo_vec = map_1DAvgInhomo_vec(idxs);
map_1DStdInhomo_vec = map_1DStdInhomo_vec(idxs);
map_vec_z_axis      = map_z_axis(idxs);

if (DEBUG_FLAG >= 3)
%     uiwait(msgbox('have you correctly rotated x & y axes of the map?'));
    figure;
    subplot(2,2,1); imagesc(input_map_mat.*input_map_mask); title('Input map');   set(gca,'ydir','normal');
    subplot(2,2,2); imagesc(map_mat      .*      map_mask); title('Rotated map'); set(gca,'ydir','normal');
    subplot(2,2,3); hold on; plot(x_axis,'b.-'); plot(y_axis,'k.-'); plot(map_z_axis,'r^');
    title('Map axes'); ylabel('[cm]'); legend({'x-axis','y-axis','max z-axis'});
    subplot(2,2,4); plot(map_vec_z_axis,'k.-'); title('Masked Map z-axis'); ylabel('[cm]');
end;

% --------------------------------------------
% Poly-fit of the frequency-shift versus space
% --------------------------------------------
DNu_of_z = map_1DAvgInhomo_vec;
z_axis   = map_vec_z_axis;

% NOT NEEDED ANY MORE SINCE WE SMOOTH THE INHOMOGENEITY VECTORS LATER ON
% Extrapolate DNu_of_z so that the edges will be smooth
nleft  = 0; % round(length(DNu_of_z)/50);
nright = 0;
if (nleft ~= 0 || nright ~= 0)
    DNu_of_z = [(linspace(1,1,nleft).^(0.5))*mean(DNu_of_z(1:2)) DNu_of_z];
    DNu_of_z = [DNu_of_z (linspace(1,1,nright).^(0.2))*mean(DNu_of_z(end:-1:end-1))];
%     DNu_of_z = window_1Dvec(DNu_of_z,2,5,1,'');
    map_1DStdInhomo_vec = [zeros(1,nleft) map_1DStdInhomo_vec zeros(1,nright)];
end;

dz_ = map_vec_z_axis(2) - map_vec_z_axis(1);
z_axis = linspace(map_vec_z_axis(1) - nleft*dz_, map_vec_z_axis(end) + nright*dz_, length(map_vec_z_axis) + nleft+nright);

P_DNu0 = polyfit(z_axis,DNu_of_z,PolyOrder);

% -----------------------------------------------------------------------------------------------
% Delete old plot from figure and plot
%  1. The rotated phase map with the inhomogeneity axis (x-axis)
%  2. The inhomogeneity poly-fit
% -----------------------------------------------------------------------------------------------
if ((nargin > 5) && (ax0 ~= 0))     % Plot the polyfit on an existing figure
    delete(get(ax0,'Children'));
    axes(ax0); hold on;
    plot(z_axis,DNu_of_z              ,'r.','markersize',1);
    plot(z_axis,polyval(P_DNu0,z_axis),'m--');
    xlabel('Spatial axis [cm]','FontSize',8);
    title(sprintf('Average 1D inhomo [Hz]\n(\\phi_{rot}=%2.0f [deg]   Poly Order=%1.0f)',rot_phi,PolyOrder),'FontSize',8);
    set(gca,'FontSize',8);
else                                % Plot the Phase-Map & polyfit on a new figure
    if (DEBUG_FLAG >= 3)
        figure;
        mean_map_mat = map_mat;
        mean_map_mat(map_mask == 0) = mean(mean(map_mat(map_mask == 1))); 
        subplot(2,2,1); imagesc(x_axis,y_axis,transpose(mean_map_mat)); set(gca,'ydir','normal'); hold on;
        title(sprintf('Phase Map [Hz] (Rotation angle = %f)',rot_phi)); xlabel('X-Axis [cm]'); ylabel('Y-Axis [cm]'); colorbar;
    
        line([x_axis(1) x_axis(end)],[0 0]);
        subplot(2,2,2); hold on;
        plot(z_axis,DNu_of_z              ,'r.');
        plot(z_axis,polyval(P_DNu0,z_axis),'m--');
        xlabel('Spatial axis [cm]','FontSize',8);
        title(sprintf('Average 1D inhomo [Hz]\n(\\phi_{rot}=%2.0f [deg]   Poly Order=%1.0f)',rot_phi,PolyOrder),'FontSize',8);
        set(gca,'FontSize',8);
    end;
end;

% -----------------------------------------------------------------------------------------------
% Delete old plot from figure and plot the inhomogeneity STD
% -----------------------------------------------------------------------------------------------
if ((nargin > 6) && (ax1 ~= 0))     % Plot the polyfit STD on an existing figure
    delete(get(ax1,'Children'));
    axes(ax1); hold on;
else                                % Plot the polyfit STD on a new figure
    if (DEBUG_FLAG < 3)
        return;
    end;
    subplot(2,2,3); hold on;
end;
plot(z_axis,map_1DStdInhomo_vec,'b-');
xlabel('Spatial axis [cm]','FontSize',8);
title(sprintf('1D inhomo STD [Hz]'),'FontSize',8);
set(gca,'FontSize',8);

return;

% - - - - -
%   JUNK
% - - - - -
% Convert the rotation angle to take into consideration the orientation
% direction of imrotate and the orientation of the phase-map matrix
% rot_phi = - ( 90 - abs(rot_phi))*sign(rot_phi)
% rot_phi = abs(rot_phi)*sign(rot_phi)

