function varargout = DICOMViewer(varargin)
% DICOMViewer M-file for DICOMViewer.fig
%      DICOMViewer, by itself, creates a new DICOMViewer or raises the existing
%      singleton*.
%
%      H = DICOMViewer returns the handle to a new DICOMViewer or the handle to
%      the existing singleton*.
%
%      DICOMViewer('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DICOMViewer.M with the given input arguments.
%
%      DICOMViewer('Property','Value',...) creates a new DICOMViewer or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before DICOMViewer_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to DICOMViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help DICOMFILES

% Last Modified by GUIDE v2.5 07-Apr-2014 11:47:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DICOMViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @DICOMViewer_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end

% End initialization code - DO NOT EDIT


% --- Executes just before DICOMFiles is made visible.
function DICOMViewer_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for DICOMFiles
handles.output = hObject;
handles.hidden = [];
guidata(hObject, handles);
colormap gray(256)

% use push-button callback
if length(varargin) && ischar (varargin{1})
   handles.dfolder = varargin{1};
   SetFolder(handles);
   ListBox_Callback(hObject, eventdata, handles);
else
%  newFolder_Callback(hObject, eventdata, handles)
   handles.dfolder = pwd;
   SetFolder(handles);
   ListBox_Callback(hObject, eventdata, handles);
end

% UIWAIT makes DICOMFiles wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = DICOMViewer_OutputFcn(hObject, eventdata, handles)
% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in newFolder.
function newFolder_Callback(hObject, eventdata, handles)

P = fileparts(mfilename('fullpath'));
% P = 'E:\01 Post\06 Projects\Diffusion';                                                         % [nbe]
% P = 'E:\01 Post\06 Projects\Prostate\900 Experiments\2013_06_19_MSK_phantom_CMC\DICOMs';        % [nbe]
P = 'E:\01 Post\06 Projects\Prostate'                                                    ;        % [nbe]
nfolder=uigetdir(P,'Select DICOM Directory');
if ~ischar(nfolder)
    disp('no valid Directory selected.')
    return;
end
handles.dfolder=nfolder;
guidata(hObject, handles);
SetFolder(handles);
ListBox_Callback(hObject, eventdata, handles)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SetFolder (handles)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dfiles=dir(handles.dfolder);
dfiles=dfiles(3:end);                   % avoid . and ..
nfiles=length(dfiles);
idx2=0;
for idx = 1:nfiles
	nm = dfiles(idx).name;
	if (length(nm) < 3) % || ((strcmp(nm(end-2:end),'IMA') == 0) && (strcmp(nm(end-2:end),'dcm') == 0))
		continue;
	end;
	idx2 = idx2 + 1;
	new_dfiles(idx2) = dfiles(idx);
end;
dfiles = new_dfiles;
nfiles = length(dfiles);

if nfiles<1
    disp('no files availabel.')
    return;
end
set(handles.ListBox,'String',char(dfiles.name),'value',1);
s = [num2str(nfiles) ' files in: ' handles.dfolder];
set(handles.NofFiles,'String', s);
guidata(handles.figure1, handles);

% --- Executes during object creation, after setting all properties.
function ListBox_CreateFcn(hObject, eventdata, handles)

if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on selection change in ListBox.
function ListBox_Callback(hObject, eventdata, handles)

fname = get(handles.ListBox,'String');
fname = fname(get(handles.ListBox,'value'),:);
try
   metadata = dicominfo([handles.dfolder '\' fname]);
catch
   disp ('apparently not a DICOM file');
   return
end
img      = dicomread([handles.dfolder '\' fname]);

zoomF   = get(handles.zoomF  , 'string'); zoomF=str2num(zoomF);  if (isempty(zoomF) || zoomF==0), zoomF=1; end;
caF     = get(handles.caF    , 'string'); caF  =str2num(  caF);  if (isempty(caF  ) || caF  ==0), caF  =1; end;
% imagesc(img);  % orig
interpF =3;
imagesc(interp_2D_mat(double(img),interpF,interpF));  % nbe
ai;
zoom reset;
if (zoomF ~= 1.0), zoom(double(zoomF)); end;
ca=caxis; caxis(ca/caF);
axis off;
% axis([200  432.5000    0.5000  432.5000]); %ai;

plotFig = get(handles.plotFig, 'value');
if plotFig
% [nbe] 2011_10_26: plot image in a separate figure
figure;
% smooth_f = 1;
% if    (smooth_f == 1), imagesc(img);
% else  imagesc(smooth_mat(img,smooth_f,smooth_f)/10);
% end;
imagesc(interp_2D_mat(double(img),interpF,interpF));
cg; ai; colorbar; imgui;  title(sprintf('%s',fname),'interpreter','none'); 
ca=caxis; caxis(ca/caF);
zoom(zoomF);
end;

ch = get(handles.HeaderList, 'value');
fields=char(fieldnames(metadata));
len = setdiff (1:size(fields,1), handles.hidden);
id=0;
for k=len,
    estr=eval(['metadata.' fields(k,:)]);
    if ischar(estr)
        str=[fields(k,:) ' : ' estr];
    elseif isnumeric(estr)
        str=[fields(k,:) ' : ' num2str(estr(1:min(3,end))')];
    else
        str=[fields(k,:) ' : ...'];
    end
    id = id+1;
    cstr{id}=sprintf('%3d %s',k,str);
end
set(handles.HeaderList,'Value',ch);
set(handles.HeaderList,'String',cstr);
guidata(hObject, handles);
return;

% --- Executes during object creation, after setting all properties.
function HeaderList_CreateFcn(hObject, eventdata, handles)
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on button press in cmdHide.
function cmdHide_Callback(hObject, eventdata, handles)
st = get (handles.HeaderList, 'string');
hide = get(handles.HeaderList, 'value');
if length(hide)==length(st)
   disp ('WARNING: at least one field must be shown');
   return
end
hidev=[];
for id=hide
   hidev = [hidev str2num(st{id}(1:3))];
end
handles.hidden = union (handles.hidden, hidev);
set(handles.HeaderList,'Value',1);
guidata (hObject, handles);
ListBox_Callback(hObject, eventdata, handles);

% --- Executes on button press in cmdShowAll.
function cmdShowAll_Callback(hObject, eventdata, handles)

handles.hidden = [];
guidata (hObject, handles);
ListBox_Callback(hObject, eventdata, handles);


% --- Executes on button press in plotFig.
function plotFig_Callback(hObject, eventdata, handles)
% hObject    handle to plotFig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotFig



function zoomF_Callback(hObject, eventdata, handles)
% hObject    handle to zoomF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zoomF as text
%        str2double(get(hObject,'String')) returns contents of zoomF as a double


% --- Executes during object creation, after setting all properties.
function zoomF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zoomF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function caF_Callback(hObject, eventdata, handles)
% hObject    handle to caF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of caF as text
%        str2double(get(hObject,'String')) returns contents of caF as a double


% --- Executes during object creation, after setting all properties.
function caF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to caF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
