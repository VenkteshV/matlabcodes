function varargout = cellcountui(varargin)
% CELLCOUNTUI MATLAB code for cellcountui.fig
%      CELLCOUNTUI, by itself, creates a new CELLCOUNTUI or raises the existing
%      singleton*.
%
%      H = CELLCOUNTUI returns the handle to a new CELLCOUNTUI or the handle to
%      the existing singleton*.
%
%      CELLCOUNTUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CELLCOUNTUI.M with the given input arguments.
%
%      CELLCOUNTUI('Property','Value',...) creates a new CELLCOUNTUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cellcountui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cellcountui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cellcountui

% Last Modified by GUIDE v2.5 09-Mar-2017 14:41:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cellcountui_OpeningFcn, ...
                   'gui_OutputFcn',  @cellcountui_OutputFcn, ...
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


% --- Executes just before cellcountui is made visible.
function cellcountui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cellcountui (see VARARGIN)

% Choose default command line output for cellcountui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes cellcountui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = cellcountui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename pathname] = uigetfile({'*.png'},'File Selector');
 handles.myImage = strcat(pathname, filename);
 axes(handles.axes1);
 imshow(handles.myImage)
 handles.MyImage = imread(handles.myImage);
 % save the updated handles object
 guidata(hObject,handles);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
 if isfield(handles,'myImage')
     A= handles.MyImage;
I= rgb2gray(A);
I = adapthisteq(I);
I = imclearborder(I);

I = wiener2(I, [5 5]);
bw = im2bw(I, graythresh(I));
bw2 = imfill(bw,'holes');
bw3 = imopen(bw2,strel('disk',2));
bw4 = bwareaopen(bw3, 100);
bw4_perim = bwperim(bw4);
overlay1 = imoverlay(I, bw4_perim, [1 .3 .3]);
maxs = imextendedmax(I, 5);
maxs = imclose(maxs ,strel('disk',3));
maxs = imfill(maxs, 'holes');
maxs = bwareaopen(maxs, 2);
overlay = imoverlay(I, bw4_perim | maxs, [1 .3 .3]);
Jc = imcomplement(I);
I_mod = imimposemin(Jc, ~bw4 | maxs);
L = watershed(I_mod);
imshow(overlay1);
labeledImage = label2rgb(L);
[L, num] = bwlabel(L);
disp(numel(L));
mask = im2bw(L, 1);
overlay3 = imoverlay(I,mask, [1 .3 .3]);
 end
