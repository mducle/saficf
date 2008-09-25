function varargout = fourpanels(varargin)
% FOURPANELS M-file for fourpanels.fig
%      FOURPANELS, by itself, creates a new FOURPANELS or raises the existing
%      singleton*.
%
%      H = FOURPANELS returns the handle to a new FOURPANELS or the handle to
%      the existing singleton*.
%
%      FOURPANELS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FOURPANELS.M with the given input arguments.
%
%      FOURPANELS('Property','Value',...) creates a new FOURPANELS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before fourpanels_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to fourpanels_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help fourpanels

% Last Modified by GUIDE v2.5 22-Sep-2008 19:47:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @fourpanels_OpeningFcn, ...
                   'gui_OutputFcn',  @fourpanels_OutputFcn, ...
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


% --- Executes just before fourpanels is made visible.
function fourpanels_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles (field names==tags) and user data (see GUIDATA)
% varargin   command line arguments to fourpanels (see VARARGIN)

% Choose default command line output for fourpanels
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes fourpanels wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = fourpanels_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in tbSinglet.
function tbSinglet_Callback(hObject, eventdata, handles)
% hObject    handle to tbSinglet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in tbDoublet.
function tbDoublet_Callback(hObject, eventdata, handles)
% hObject    handle to tbDoublet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in tbQuartet.
function tbQuartet_Callback(hObject, eventdata, handles)
% hObject    handle to tbQuartet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in tbTriplet.
function tbTriplet_Callback(hObject, eventdata, handles)
% hObject    handle to tbTriplet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in savepars.
function savepars_Callback(hObject, eventdata, handles)
% hObject    handle to savepars (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in loadpars.
function loadpars_Callback(hObject, eventdata, handles)
% hObject    handle to loadpars (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in savespec.
function savespec_Callback(hObject, eventdata, handles)
% hObject    handle to savespec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in savelevels.
function savelevels_Callback(hObject, eventdata, handles)
% hObject    handle to savelevels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in savephys.
function savephys_Callback(hObject, eventdata, handles)
% hObject    handle to savephys (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in susc.
function susc_Callback(hObject, eventdata, handles)
% hObject    handle to susc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function Tcalc_Callback(hObject, eventdata, handles)
% hObject    handle to Tcalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tcalc as text
%        str2double(get(hObject,'String')) returns contents of Tcalc as a double


% --- Executes during object creation, after setting all properties.
function Tcalc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tcalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in invsusc.
function invsusc_Callback(hObject, eventdata, handles)
% hObject    handle to invsusc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in heatcap.
function heatcap_Callback(hObject, eventdata, handles)
% hObject    handle to heatcap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in mag.
function mag_Callback(hObject, eventdata, handles)
% hObject    handle to mag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function Bcalc_Callback(hObject, eventdata, handles)
% hObject    handle to Bcalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Bcalc as text
%        str2double(get(hObject,'String')) returns contents of Bcalc as a double


% --- Executes during object creation, after setting all properties.
function Bcalc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Bcalc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function Emin_Callback(hObject, eventdata, handles)
% hObject    handle to Emin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Emin as text
%        str2double(get(hObject,'String')) returns contents of Emin as a double


% --- Executes during object creation, after setting all properties.
function Emin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Emin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function Emax_Callback(hObject, eventdata, handles)
% hObject    handle to Emax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Emax as text
%        str2double(get(hObject,'String')) returns contents of Emax as a double


% --- Executes during object creation, after setting all properties.
function Emax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Emax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in estimatediff.
function estimatediff_Callback(hObject, eventdata, handles)
% hObject    handle to estimatediff (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function Jval_Callback(hObject, eventdata, handles)
% hObject    handle to Jval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Jval as text
%        str2double(get(hObject,'String')) returns contents of Jval as a double


% --- Executes during object creation, after setting all properties.
function Jval_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Jval (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on button press in fitengy.
function fitengy_Callback(hObject, eventdata, handles)
% hObject    handle to fitengy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in fitspec.
function fitspec_Callback(hObject, eventdata, handles)
% hObject    handle to fitspec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in tSymm.
function tSymm_Callback(hObject, eventdata, handles)
% hObject    handle to tSymm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns tSymm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tSymm


% --- Executes during object creation, after setting all properties.
function tSymm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tSymm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


