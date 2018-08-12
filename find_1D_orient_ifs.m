function varargout = find_1D_orient_ifs(varargin)
% FIND_1D_ORIENT_IFS M-file for find_1D_orient_ifs.fig
%      FIND_1D_ORIENT_IFS, by itself, creates a new FIND_1D_ORIENT_IFS or raises the existing
%      singleton*.
%
%      H = FIND_1D_ORIENT_IFS returns the handle to a new FIND_1D_ORIENT_IFS or the handle to
%      the existing singleton*.
%
%      FIND_1D_ORIENT_IFS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FIND_1D_ORIENT_IFS.M with the given input arguments.
%
%      FIND_1D_ORIENT_IFS('Property','Value',...) creates a new FIND_1D_ORIENT_IFS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before find_1D_orient_ifs_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to find_1D_orient_ifs_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help find_1D_orient_ifs

% Last Modified by GUIDE v2.5 03-Feb-2008 14:27:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @find_1D_orient_ifs_OpeningFcn, ...
                   'gui_OutputFcn',  @find_1D_orient_ifs_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% ==============================================================================

% --- Executes just before find_1D_orient_ifs is made visible.
function find_1D_orient_ifs_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for find_1D_orient_ifs
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes find_1D_orient_ifs wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Delete all current figures
delete_all_axes(handles);

UserData.ge_fn         = [];
UserData.map_fn        = [];
UserData.ge_mat        = [];
UserData.map_mat       = [];
UserData.SmoothF       = 1;
UserData.dfSmoothF     = 1;
UserData.TrHld1        = 0;
UserData.TrHld2        = 100;
UserData.TrHld3        = 100;
UserData.ForcedMaxPhi  = -999;
UserData.PolyOrder     = 6;
UserData.OrientMatFlag = 0;
UserData.VOI_x1        = [];
UserData.VOI_x2        = [];
UserData.VOI_y1        = [];
UserData.VOI_y2        = [];
UserData.VOI_phi       = [];
set(handles.SmoothF    ,'string',num2str(   0));
set(handles.dfSmoothF  ,'string',num2str(   0));
set(handles.PolyOrder  ,'string',num2str(   0));
set(handles.TrHld1     ,'string',num2str(   0));
set(handles.TrHld2     ,'string',num2str( 100));
set(handles.TrHld3     ,'string',num2str( 100));
set(handles.MaxPhi     ,'string',num2str(-999));
set(handles.MinPhi     ,'string',num2str(-999));
set(handles.ForceMaxPhi,'string',num2str(   0));

set(handles.figure1,'UserData',UserData);

function delete_all_axes(handles)
delete(get(handles.axes1 ,'Children')); axes(handles.axes1 ); xlabel(''); ylabel(''); title('');
delete(get(handles.axes2 ,'Children')); axes(handles.axes2 ); xlabel(''); ylabel(''); title('');
delete(get(handles.axes3 ,'Children')); axes(handles.axes3 ); xlabel(''); ylabel(''); title('');
delete(get(handles.axes4 ,'Children')); axes(handles.axes4 ); xlabel(''); ylabel(''); title('');
delete(get(handles.axes5 ,'Children')); axes(handles.axes5 ); xlabel(''); ylabel(''); title('');
delete(get(handles.axes6 ,'Children')); axes(handles.axes6 ); xlabel(''); ylabel(''); title('');
delete(get(handles.axes7 ,'Children')); axes(handles.axes7 ); xlabel(''); ylabel(''); title('');
delete(get(handles.axes8 ,'Children')); axes(handles.axes8 ); xlabel(''); ylabel(''); title('');
delete(get(handles.axes9 ,'Children')); axes(handles.axes9 ); xlabel(''); ylabel(''); title('');
delete(get(handles.axes11,'Children')); axes(handles.axes11); xlabel(''); ylabel(''); title('');

% --- Outputs from this function are returned to the command line.
function varargout = find_1D_orient_ifs_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;

% ==============================================================================

function FileName_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of FileName as text
%        str2double(get(hObject,'String')) returns contents of FileName as a double
map_fn = get(hObject,'String');
ge_fn  = [map_fn(1:end-7) 'GE.mat'];
if (exist(map_fn,'file') ~= 2)
    errordlg(sprintf('Map File (%s) not found',map_fn),'File Error');
    UD.map_fn = [];
    UD.ge_fn  = [];
    return;
end;
if (exist(ge_fn,'file') ~= 2)
    errordlg(sprintf('GE File (%s) not found',ge_fn),'File Error');
    UD.map_fn = [];
    UD.ge_fn  = [];
    return;
end;

UD = get(handles.figure1,'UserData');
UD.map_fn = map_fn;
UD.ge_fn  = ge_fn;
set(handles.figure1,'UserData',UD);

load_and_plot_maps(handles,1)

% --- Executes during object creation, after setting all properties.
function FileName_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% ==============================================================================

function FileBrowser_Callback(hObject, eventdata, handles)
UD = get(handles.figure1,'UserData');
cur_pwd = pwd;
cd('D:\PhD\06 Experiments\T2 and Phase Mapping');
[map_fn, map_path] = uigetfile('*.mat','Enter Phase Map Filename');
UD.map_fn = [map_path map_fn];
UD.ge_fn  = [map_fn(1:end-7) 'GE.mat'];
if (exist(UD.ge_fn,'file') ~= 2)
    errordlg(sprintf('GE File (%s) not found',UD.ge_fn),'File Error');
    UD.map_fn = [];
    UD.ge_fn  = [];
    set(handles.figure1,'UserData',UD);
    return;
end;
cd(cur_pwd);
set(handles.FileName,'string',UD.map_fn);
set(handles.figure1,'UserData',UD);

load_and_plot_maps(handles,1);

% ==============================================================================

function TrHld1_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of TrHld1 as text
%        str2double(get(hObject,'String')) returns contents of TrHld1 as a double
UD = get(handles.figure1,'UserData');
if (isempty(UD.map_fn) || isempty(UD.ge_fn))
	uiwait(msgbox('Please choose map file'));
	return;
end;
TrHld1 = str2num(get(hObject,'String'));
if ((TrHld1 < 0) || (TrHld1 > 100))
	uiwait(msgbox('Please enter a valid threshold value [0..100]'));
    set(handles.TrHld1,'value',0);
    set(handles.TrHld1,'string',num2str(0));
	return;
end;
UD.TrHld1 = TrHld1;
UD = apply_filter_mat1(UD);
set(handles.figure1,'UserData',UD);
load_and_plot_maps(handles);

% ---------------------------------
function UD = apply_filter_mat1(UD)
% Set the conversion / filter matrix
M = abs(UD.ge_mat);
mn = min(min(M));
mx = max(max(M));
cm1 = ones(size(UD.ge_mat ,1),size(UD.ge_mat,2));    % just as a dummy starting value
cm1(M <  (mn + (mx-mn)*UD.TrHld1/100)) = 0;
cm1(M >= (mn + (mx-mn)*UD.TrHld1/100)) = 1;
UD.filter_mat1 = cm1;
return;

% ----------------------------------------------------
function TrHld1_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end;

% ==============================================================================

function TrHld2_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of TrHld2 as text
%        str2double(get(hObject,'String')) returns contents of TrHld2 as a double
UD = get(handles.figure1,'UserData');
if (isempty(UD.map_fn) || isempty(UD.ge_fn))
	uiwait(msgbox('Please choose map file'));
	return;
end;
TrHld2 = str2num(get(hObject,'String'));
if (TrHld2 < 0) || (TrHld2 > 100)
	uiwait(msgbox('Please enter a valid threshold value [0..100]'));
    set(handles.TrHld2,'value',100);
    set(handles.TrHld2,'string',num2str(100));
	return;
end;
UD.TrHld2 = TrHld2;
UD = apply_filter_mat2(UD);
set(handles.figure1,'UserData',UD);
load_and_plot_maps(handles);

% ---------------------------------
function UD = apply_filter_mat2(UD)
M = abs(UD.ge_mat);

% Set the conversion / filter matrix -- derivative along rows
df_row = diff(M,1,1);   df_row(end+1,:) = df_row(end,:);
mn = min(min(df_row));
mx = max(max(df_row));
cm2 = ones(size(UD.ge_mat ,1),size(UD.ge_mat,2));    % just as a dummy starting value
cm2(df_row >  (mn + (mx-mn)*UD.TrHld2/100)) = 0;
cm2(df_row <= (mn + (mx-mn)*UD.TrHld2/100)) = 1;
UD.filter_mat2 = cm2;

% Set the conversion / filter matrix -- derivative along columns
df_col = diff(M,1,2);   df_col(:,end+1) = df_col(:,end);
mn = min(min(df_col));
mx = max(max(df_col));
cm2 = ones(size(UD.ge_mat ,1),size(UD.ge_mat,2));    % just as a dummy starting value
cm2(df_col >  (mn + (mx-mn)*UD.TrHld2/100)) = 0;
cm2(df_col <= (mn + (mx-mn)*UD.TrHld2/100)) = 1;
UD.filter_mat2 = UD.filter_mat2 .* cm2;  % retain the previous threshold information

% ----------------------------------------------------
function TrHld2_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end;

% ==============================================================================

function TrHld3_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of TrHld3 as text
%        str2double(get(hObject,'String')) returns contents of TrHld3 as a double
UD = get(handles.figure1,'UserData');
if (isempty(UD.map_fn) || isempty(UD.ge_fn))
	uiwait(msgbox('Please choose map file'));
	return;
end;
TrHld3 = str2num(get(hObject,'String'));
if ((TrHld3 < 0) || (TrHld3 > 100))
	uiwait(msgbox('Please enter a valid threshold value [0..100]'));
    set(handles.TrHld3,'value',100);
    set(handles.TrHld3,'string',num2str(100));
	return;
end;
UD.TrHld3 = TrHld3;
UD = apply_filter_mat3(UD,1);
set(handles.figure1,'UserData',UD);
load_and_plot_maps(handles);

% ----------------------------------------
% Calculate the conversion / filter matrix
function UD = apply_filter_mat3(UD,reset)
[dfdx, dfdy] = calc_map_1D_derivatives(UD.map_mat);

if ((nargin > 1) && reset)
    UD.filter_mat3 = ones(size(UD.map_mat,1),size(UD.map_mat,2));    % just as a dummy starting value
end;

% We are only interested in the absolute value of the phase change, thus:
dfdx = abs(dfdx);
dfdy = abs(dfdy);

% Set the conversion / filter matrix -- derivative along x-axis
mn = min(min(dfdx));
mx = max(max(dfdx));
cm3 = ones(size(UD.map_mat,1),size(UD.map_mat,2));    % just as a dummy starting value
cm3(dfdx >  (mn + (mx-mn)*UD.TrHld3/100)) = 0;
cm3(dfdx <= (mn + (mx-mn)*UD.TrHld3/100)) = 1;
UD.filter_mat3 = UD.filter_mat3 .* cm3;

% Set the conversion / filter matrix -- derivative along y-axis
mn = min(min(dfdy));
mx = max(max(dfdy));
cm3 = ones(size(UD.map_mat,1),size(UD.map_mat,2));    % just as a dummy starting value
cm3(dfdy >  (mn + (mx-mn)*UD.TrHld3/100)) = 0;
cm3(dfdy <= (mn + (mx-mn)*UD.TrHld3/100)) = 1;
UD.filter_mat3 = UD.filter_mat3 .* cm3;  % retain the previous threshold information

% ----------------------------------------------------
function TrHld3_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end;

% ==============================================================================

function SmoothF_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of SmoothF as text
%        str2double(get(hObject,'String')) returns contents of SmoothF as a double
UD = get(handles.figure1,'UserData');
if (isempty(UD.map_fn) || isempty(UD.ge_fn))
	uiwait(msgbox('Please choose map file'));
    set(handles.SmoothF,'string','1');
	return;
end;
SmoothF = str2num(get(hObject,'String'));
if (SmoothF ~= round(SmoothF))
	uiwait(msgbox('Please enter integer value'));
	return;
end;
UD.SmoothF = SmoothF;

% Apply the phase map derivative filter before and after the smoothing
UD.map_mat = UD.Orig_map_mat;
UD = apply_filter_mat3(UD,1);

UD.ge_mat  = smooth_mat(UD.Orig_ge_mat , SmoothF);
UD.map_mat = smooth_mat(UD.Orig_map_mat, SmoothF);

UD = apply_filter_mat1(UD);
UD = apply_filter_mat2(UD);
UD = apply_filter_mat3(UD);
set(handles.figure1,'UserData',UD);
load_and_plot_maps(handles);

function SmoothF_CreateFcn(hObject, eventdata, handles)
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end;

% --------------------------------------------
function out_mat = smooth_mat(in_mat, SmoothF)
CONV_MAT = ones(SmoothF,SmoothF) / sum(sum(ones(SmoothF,SmoothF)));
out_mat = conv2(in_mat,CONV_MAT,'same');
return;

% ==============================================================================

function dfSmoothF_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'String') returns contents of dfSmoothF as text
%        str2double(get(hObject,'String')) returns contents of dfSmoothF as a double
UD = get(handles.figure1,'UserData');
if (isempty(UD.map_fn) || isempty(UD.ge_fn))
	uiwait(msgbox('Please choose map file'));
    set(handles.dfSmoothF,'string','1');
	return;
end;
dfSmoothF = str2num(get(hObject,'String'));
if (dfSmoothF ~= round(dfSmoothF))
	uiwait(msgbox('Please enter integer value'));
	return;
end;
UD.dfSmoothF = dfSmoothF;
set(handles.figure1,'UserData',UD);
load_and_plot_maps(handles);

function dfSmoothF_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end;

% ==============================================================================

function ForceMaxPhi_Callback(hObject, eventdata, handles)
UD = get(handles.figure1,'UserData');
if (isempty(UD.map_fn) || isempty(UD.ge_fn))
	uiwait(msgbox('Please choose map file'));
    set(handles.ForceMaxPhi,'string','--');
	return;
end;
UD.ForcedMaxPhi = str2num(get(hObject,'String'));
set(handles.figure1,'UserData',UD);
CalcOrient_Callback(0, 0, handles);
return;

function ForceMaxPhi_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end;
return;

% ==============================================================================

function CalcOrient_Callback(hObject, eventdata, handles)
UD = get(handles.figure1,'UserData');
if (isempty(UD.map_fn) || isempty(UD.ge_fn))
	uiwait(msgbox('Please choose map file'));
	return;
end;
M  = UD.map_mat;
fm1 = UD.filter_mat1; fm2 = UD.filter_mat2; fm3 = UD.filter_mat3; fm4 = UD.filter_mat_VOI;
fm = fm1 .* fm2 .* fm3 .* fm4;

if (UD.ForcedMaxPhi == -999)
    % Calculate 'P': the phase derivatives tensor -
    % 1. perform the calculations on the original (yet, smoothed) phase map
    % 2. SKIP, we need the signal sign: Take only the absolute value
    %    (we are not interested in the sign of change, only in it's size)
    % 3. Smooth the derivative matrix
    % 4. Do the trunction
    % 5. Remove the matrix edges and relating edge effects
    [dfdx, dfdy] = calc_map_1D_derivatives(M);
    % dfdx = abs(dfdx); NO NO
    dfdx = smooth_mat(dfdx,UD.dfSmoothF);   dfdx = dfdx.*fm;   dfdx = dfdx(3:end-3,3:end-3);
    % dfdy = abs(dfdy); NO NO
	dfdy = smooth_mat(dfdy,UD.dfSmoothF);   dfdy = dfdy.*fm;   dfdy = dfdy(3:end-3,3:end-3);
	dxdy = dfdy + dfdx;
	df = sqrt((dfdy.^2) + (dfdx.^2));% .* sign(dxdy);

    imagesc_mat(handles.axes6,[],[],dfdx - mean(mean(dfdx)),0);   set(gca,'xTickLabel',[],'yTickLabel',[]);
    imagesc_mat(handles.axes7,[],[],dfdy - mean(mean(dfdy)),0);   set(gca,'xTickLabel',[],'yTickLabel',[]);
	
	if (0)
	figure; caxis_F = 4;
	subplot(2,2,1); imagesc(df  );      set(gca,'ydir','normal'); title('sqrt(dfdx^2 + dfdy^2)'); colorbar; a = caxis; caxis(a/caxis_F); caxis([0 15]);
	subplot(2,2,2); imagesc(abs(dfdx)); set(gca,'ydir','normal'); title('|dfdx|');                colorbar; a = caxis; caxis(a/caxis_F); caxis([0  4]);
	subplot(2,2,3); imagesc(abs(dfdy)); set(gca,'ydir','normal'); title('|dfdy|');                colorbar; a = caxis; caxis(a/caxis_F); caxis([0 15]);
	subplot(2,2,4); imagesc(abs(dxdy)); set(gca,'ydir','normal'); title('|dfdx + dfdy|');         colorbar; a = caxis; caxis(a/caxis_F); caxis([0 16]);
		
	interpF = 1;
	[xi,yi] = meshgrid(linspace(1,size(M,2),size(M,2)*interpF),linspace(1,size(M,1),size(M,1)*interpF));
	M_interp = interp2(M.*fm,xi,yi);
	
	reductionF = 6;
	for idx1 = 1:size(M,1)/reductionF
		row1 = reductionF*interpF*(idx1-1) + 1;
		row2 = reductionF*interpF*(idx1);
		for idx2 = 1:size(M,2)/reductionF
			col1 = reductionF*interpF*(idx2-1) + 1;
			col2 = reductionF*interpF*(idx2);
			tmp  = M_interp(row1:row2,col1:col2);
			M_dephase1(idx1,idx2) = abs(sum(sum(exp(2*pi*i*M_interp(row1:row2,col1:col2)))));
			M_dephase2(idx1,idx2) = exp(-1*(std(tmp(:))));
		end;
	end;
	M_dephase1 = M_dephase1 / max(max(M_dephase1));
% 	M_dephase = 1 - M_dephase;
	figure;
	subplot(2,2,1); imagesc(M.*fm);      set(gca,'ydir','normal'); title('M');            colorbar;
	subplot(2,2,2); imagesc(M_interp);   set(gca,'ydir','normal'); title('M_{interp}');   colorbar;
	subplot(2,2,3); imagesc(M_dephase1); set(gca,'ydir','normal'); title('M_{dephase}1'); colorbar;
	subplot(2,2,4); imagesc(M_dephase2); set(gca,'ydir','normal'); title('M_{dephase}2'); colorbar;
	end;
	
    P(1,1) = sum(sum(dfdx.*dfdx));
    P(1,2) = sum(sum(dfdx.*dfdy));
    P(2,1) = P(1,2);
    P(2,2) = sum(sum(dfdy.*dfdy));

    % Calculate the eigenvalues and eigenvectors of P
    %    Diagonal of D are the eigenvalues
    %    Columns  of V are the eigenvectors
    [V,D] = eigs(P);
    [maxD,maxIdx] = max([D(1,1), D(2,2)]);   maxV = V(:,maxIdx);
    [minD,minIdx] = min([D(1,1), D(2,2)]);   minV = V(:,minIdx);
    max_phi1 = atan2(maxV(2),maxV(1)) * 180 / pi;     % Find the angle, and it's relevant quadrant
    max_phi1 = mod(270+max_phi1,180) - 90;            % Convert the angle to [-90..+90] region
    min_phi1 = atan2(minV(2),minV(1)) * 180 / pi;
    min_phi1 = mod(270+min_phi1,180) - 90;
else
    max_phi1 = UD.ForcedMaxPhi;
    min_phi1 = -999;
end;
set(handles.MaxPhi,'string',num2str(max_phi1));
set(handles.MinPhi,'string',num2str(min_phi1));

% Plot filtered Phase map
delete(get(handles.axes2 ,'Children'));
new_map = UD.map_mat;
mean_map_val = mean(mean(new_map(fm == 1)));
new_map(fm == 0) = mean_map_val;
% new_map = new_map - mean(mean(new_map)); % redundant -- already implicitly done in the previous line
imagesc_mat(handles.axes2,UD.map_x,UD.map_y,new_map,0); % caxis([-80 0]);
xlabel('X-axis [cm]','Fontsize',8);
% ylabel('Y-axis [cm]','Fontsize',8);

% Debugging:
figure;
imagesc(UD.map_x,UD.map_y,new_map);
xlabel('X-axis [cm]','Fontsize',8);

GEmat = UD.ge_mat;
figure;
imagesc(UD.map_x,UD.map_y,abs(GEmat));
xlabel('X-axis [cm]','Fontsize',8);
mean_GEmat_val = mean(mean(GEmat(fm == 1)));
GEmat(fm == 0) = mean_GEmat_val;
uiwait(msgbox('interp GEmap figure?')); % GEmat = interp_2D_mat(GEmat,2);
figure;
imagesc(UD.map_y,UD.map_x,transpose(abs(GEmat)));
% xlabel('X-axis [cm]','Fontsize',8);

% Plot Inhomogeneity orientation on Phase map
%  max_phi denotes counter-clockwise rotation from the x-axis [-x --> +x]
a_max = tan(max_phi1*pi/180);
a_min = tan(min_phi1*pi/180);
if (a_max == 0)
    x_max = UD.map_x;
    y_max = zeros(1,length(x_max));
else
    y_max = UD.map_y;
    x_max = y_max / a_max;
end;
if (a_min == 0)
    x_min = UD.map_x;
    y_min = zeros(1,length(x_min));
else
    y_min = UD.map_y;
    x_min = y_min / a_min;
end;
axes(handles.axes2); hold on; plot(x_max,y_max,'w-'); if (min_phi1 ~= -999) plot(x_min,y_min,'w--'); end;

% Save the Phase Map and calculated orientation in a MAT file for later use
UD.mat_fn = sprintf('%sOrient_RP%s_SF%s_DSF%s_PO%s',UD.map_fn(1:end-4)      ,...
                                                    num2str(round(max_phi1)),...
                                                    num2str(UD.SmoothF)     ,...
                                                    num2str(UD.dfSmoothF)   ,...
                                                    num2str(UD.PolyOrder));
map_mat  = UD.map_mat;
map_mask = fm;
rot_phi  = max_phi1;
x_axis   = UD.map_x;
y_axis   = UD.map_y;

if (UD.OrientMatFlag)
    %                 [Hz]     [none]     [cm]     [cm]     [deg]
    save(UD.mat_fn,'map_mat','map_mask','x_axis','y_axis','rot_phi');
end;

% Calculate the average inhomogeneity along the maximal-change orientation.
% DEBUG % m = map_mat;
% DEBUG % m(map_mask==0) = mean(mean(m(map_mask==1)));
% DEBUG % figure; 
% DEBUG % subplot(2,2,1); imagesc(m);                            title(sprintf('Phase map (rot = %f)',0));
% DEBUG % subplot(2,2,2); imagesc(imrotate(m,-max_phi1,'crop')); title(sprintf('Phase map (rot = -max_phi1 = %f)',-max_phi1));
[z_axis,P_DNu0,Porder] = calc_oriented_1D_inhomo(UD.map_mat,map_mask,x_axis,y_axis,max_phi1,handles.axes9,handles.axes11,UD.PolyOrder);
axes(handles.axes10); plot((0:length(P_DNu0)-1),P_DNu0(end:-1:1),'.-'); title('Poly-fit Coeff. Vs. n');
set(gca,'fontsize',8,'YTick',(0));

return;

% ------------------------------------------------------
function [dfdx, dfdy] = calc_map_1D_derivatives(M)
% We note that since M is 
dfdx = diff(M,1,1);   dfdx(end+1,:) = mean(mean(dfdx));   % dfdx(end,:);    % Diff along rows
dfdy = diff(M,1,2);   dfdy(:,end+1) = mean(mean(dfdy));   % dfdy(:,end);    % Diff along cols

% ==============================================================================

function MaxPhi_Callback(hObject, eventdata, handles)

function MaxPhi_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end;

% ==============================================================================

function MinPhi_Callback(hObject, eventdata, handles)

function MinPhi_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end;

% ==============================================================================

function PolyOrder_Callback(hObject, eventdata, handles)
UD = get(handles.figure1,'UserData');
UD.PolyOrder = str2num(get(hObject,'String'));
set(handles.figure1,'UserData',UD);

function PolyOrder_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end;

% ==============================================================================

function OrientMatFlag_Callback(hObject, eventdata, handles)
UD = get(handles.figure1,'UserData');
UD.OrientMatFlag = get(handles.OrientMatFlag,'value');
set(handles.figure1,'UserData',UD);

% ==============================================================================

function VOI_x1_Callback(hObject, eventdata, handles)
UD = get(handles.figure1,'UserData');
UD.VOI_x1 = str2num(get(hObject,'String'));
% UD = apply_filter_VOI(UD);
set(handles.figure1,'UserData',UD);
% load_and_plot_maps(handles);
return;

function VOI_x1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end;
return;

function VOI_x2_Callback(hObject, eventdata, handles)
UD = get(handles.figure1,'UserData');
UD.VOI_x2 = str2num(get(hObject,'String'));
% UD = apply_filter_VOI(UD);
set(handles.figure1,'UserData',UD);
% load_and_plot_maps(handles);
return;

function VOI_x2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end;
return;

function VOI_y1_Callback(hObject, eventdata, handles)
UD = get(handles.figure1,'UserData');
UD.VOI_y1 = str2num(get(hObject,'String'));
% UD = apply_filter_VOI(UD);
set(handles.figure1,'UserData',UD);
% load_and_plot_maps(handles);
return;

function VOI_y1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end;
return;

function VOI_y2_Callback(hObject, eventdata, handles)
UD = get(handles.figure1,'UserData');
UD.VOI_y2 = str2num(get(hObject,'String'));
% UD = apply_filter_VOI(UD);
set(handles.figure1,'UserData',UD);
% load_and_plot_maps(handles);
return;

function VOI_y2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end;
return;

% ==============================================================================

function VOI_phi_Callback(hObject, eventdata, handles)
UD = get(handles.figure1,'UserData');
UD.VOI_phi = str2num(get(hObject,'String'));
UD = apply_filter_VOI(UD);
set(handles.figure1,'UserData',UD);
load_and_plot_maps(handles);
return;

function VOI_phi_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end;
return;

% ==============================================================================

function UD = apply_filter_VOI(UD)
% Set the conversion / filter matrix
M = abs(UD.ge_mat);
% [dummy  left_col] = min(abs(UD.ge_x - UD.VOI_x1));
% [dummy right_col] = min(abs(UD.ge_x - UD.VOI_x2));
% [dummy   low_row] = min(abs(UD.ge_y - UD.VOI_y1));
% [dummy  high_row] = min(abs(UD.ge_y - UD.VOI_y2));

if (isempty(UD.VOI_x1))
    UD.VOI_x1  = min(UD.ge_x);
    UD.VOI_x2  = max(UD.ge_x);
    UD.VOI_y1  = min(UD.ge_y);
    UD.VOI_y2  = max(UD.ge_y);
    UD.VOI_phi = 0;
end;

x_axis = UD.ge_x;
y_axis = UD.ge_y;
x1 = UD.VOI_x1;
x2 = UD.VOI_x2;
y1 = UD.VOI_y1;
y2 = UD.VOI_y2;
phi = UD.VOI_phi*pi/180;
for x_idx = 1:length(x_axis)
    for y_idx = 1:length(y_axis)
        x = x_axis(x_idx);
        y = y_axis(y_idx);
        [new_x, new_y] = rot_xy_Point(x,y,phi);
        VOI_mask_mat(x_idx,y_idx) = ((new_x > x1) && (new_x < x2) && (new_y > y1) && (new_y < y2));
    end;
end;

% figure;
% subplot(1,2,1); imagesc(VOI_mask_mat); title('original map mask');
% VOI_mask_mat(425:615,:) = 0;
% subplot(1,2,2); imagesc(VOI_mask_mat); title('manually filtered map mask');

UD.filter_mat_VOI = VOI_mask_mat;
return;

% ==============================================================================

function load_and_plot_maps(handles,reset)
UD = get(handles.figure1,'UserData');
if (isempty(UD.map_fn) || isempty(UD.ge_fn))
    error('Invoked plot function on empty User Data');
end;

if (isempty(UD.map_mat) || ((nargin > 1) && reset)) % If it's the first time -- load map and reset params
    
    % Shape3D -- rows(1:end) are x-axis (1st coordiante)
    %            cols(1:end) are y-axis (2nd coordiante)
    % Shape3D(x(1),y(1)) gives the bottom left element in the image when presented as:
    % >> imagesc(x,y,transpose(Shape3D))
    load(UD.ge_fn);
    UD.Orig_ge_mat = Shape3D;
    UD.ge_mat      = Shape3D;
    UD.ge_x        = X;
    UD.ge_y        = Y;
    
    load(UD.map_fn);

    % Test #1 >>>>
%     [nrows,ncols] = size(Shape3D);
%     Shape3D = ((1:nrows)') * ones(1,ncols);
%     amp = max(max(Shape3D)) - min(min(Shape3D));
%     Shape3D = Shape3D + (amp * rand(nrows) / 2);   % add noise. assume matrix is square.
%     Shape3D = imrotate(Shape3D,40,'nearest','crop');
%     UD.ge_mat      = imrotate(ones(nrows,ncols),40,'nearest','crop');
%     UD.Orig_ge_mat = imrotate(ones(nrows,ncols),40,'nearest','crop');
    %      <<<<

    % Test #2 >>>>
%     [nrows,ncols] = size(Shape3D);
%     gwin = gausswin(ncols,1.5);
%     Shape3D = ((1:nrows)') * (gwin');
%     amp = max(max(Shape3D)) - min(min(Shape3D));
% %     Shape3D = Shape3D + (amp * rand(nrows) / 2);   % add noise. assume matrix is square.
%     Shape3D = imrotate(Shape3D,45,'nearest','crop');
%     UD.ge_mat      = imrotate(ones(nrows,ncols),45,'nearest','crop');
%     UD.Orig_ge_mat = imrotate(ones(nrows,ncols),45,'nearest','crop');
    %        <<<<

    % Test #3 >>>>
%     [nrows,ncols] = size(Shape3D);
%     [xm,ym] = ndgrid(X,Y);
%     Shape3D = xm.^2 - sin(5*xm) - 3*ym.^2 + ym;
%     amp = max(max(Shape3D)) - min(min(Shape3D));
% %     Shape3D = Shape3D + (amp * rand(nrows) / 2);   % add noise. assume matrix is square.
%     Shape3D = imrotate(Shape3D,45,'nearest','crop');
%     UD.ge_mat      = imrotate(ones(nrows,ncols),45,'nearest','crop');
%     UD.Orig_ge_mat = imrotate(ones(nrows,ncols),45,'nearest','crop');
    %        <<<<

    UD.Orig_map_mat = Shape3D*1E+3;  % [kHz --> Hz]
    UD.map_mat      = Shape3D*1E+3;
    UD.map_x        = X;
    UD.map_y        = Y;

    % Apply the phase map derivatives filter before and after the smoothing
    UD = apply_filter_mat3(UD,1);
    
    % Smooth the newly loaded map, according to current smoothing parameters
    SmoothF = UD.SmoothF;
    UD.ge_mat  = smooth_mat(UD.Orig_ge_mat , SmoothF);
    UD.map_mat = smooth_mat(UD.Orig_map_mat, SmoothF);

    % Re-calculate the truncation filters for the new maps
    UD = apply_filter_mat1(UD);
    UD = apply_filter_mat2(UD);
    UD = apply_filter_mat3(UD);
    UD = apply_filter_VOI(UD);
end;
imagesc_mat(handles.axes1,UD.ge_x ,UD.ge_y ,abs(UD.ge_mat),0);  xlabel('X-axis [cm]','FontSize',8);  ylabel('Y-axis [cm]','FontSize',8);
imagesc_mat(handles.axes2,UD.map_x,UD.map_y,    UD.map_mat,0);  xlabel('X-axis [cm]','FontSize',8);  ylabel(''); set(gca,'yTickLabel',[]);

% Reset thresholds and plot conversion / filter maps
fm1 = UD.filter_mat1;
fm2 = UD.filter_mat2;
fm3 = UD.filter_mat3;
fm4 = UD.filter_mat_VOI;
fm  = fm1 .* fm2 .* fm3 .* fm4;
imagesc_mat(handles.axes3,[],[],fm1);   set(gca,'xTickLabel',[],'yTickLabel',[]);
imagesc_mat(handles.axes4,[],[],fm2);   set(gca,'xTickLabel',[],'yTickLabel',[]);
imagesc_mat(handles.axes8,[],[],fm3);   set(gca,'xTickLabel',[],'yTickLabel',[]);
imagesc_mat(handles.axes5,[],[],fm);    set(gca,'xTickLabel',[],'yTickLabel',[]);
delete(get(handles.axes6 ,'Children'));
delete(get(handles.axes7 ,'Children'));
delete(get(handles.axes9 ,'Children'));
delete(get(handles.axes11,'Children'));

set(handles.figure1,'UserData',UD);

function imagesc_mat(ax,x,y,z,hold_flag)
if (nargin > 4)
    if (hold_flag)
        hold on;
    else
        hold off;
    end;
end;

delete(get(ax,'Children'));
axes(ax);   imagesc(x,y,transpose(z));  set(ax,'ydir','normal');


% ================
%      JUNK
% ================

%   UD.SmoothF   =   1;   set(handles.SmoothF  ,'value',  1);   set(handles.SmoothF  ,'string',num2str(  1));
%   UD.dfSmoothF =   1;   set(handles.dfSmoothF,'value',  1);   set(handles.dfSmoothF,'string',num2str(  1));
%   UD.filter_mat1 = ones(size(UD.ge_mat ,1),size(UD.ge_mat ,2));
%   UD.filter_mat2 = ones(size(UD.ge_mat ,1),size(UD.ge_mat ,2));
%   UD.filter_mat3 = ones(size(UD.map_mat,1),size(UD.map_mat,2));    
%   UD.TrHld1    =   0;   set(handles.TrHld1   ,'value',  0);   set(handles.TrHld1   ,'string',num2str(  0));
%   UD.TrHld2    = 100;   set(handles.TrHld2   ,'value',100);   set(handles.TrHld2   ,'string',num2str(100));
%   UD.TrHld3    = 100;   set(handles.TrHld3   ,'value',100);   set(handles.TrHld3   ,'string',num2str(100));

