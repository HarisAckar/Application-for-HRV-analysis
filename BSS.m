function varargout = BSS(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BSS_OpeningFcn, ...
                   'gui_OutputFcn',  @BSS_OutputFcn, ...
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


% --- Executes just before BSS is made visible.
function BSS_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BSS (see VARARGIN)

% Choose default command line output for BSS
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BSS wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = BSS_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function ime_datoteke_Callback(hObject, eventdata, handles)
% hObject    handle to ime_datoteke (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ime_datoteke as text
%        str2double(get(hObject,'String')) returns contents of ime_datoteke as a double


% --- Executes during object creation, after setting all properties.
function ime_datoteke_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ime_datoteke (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ucitaj.
function ucitaj_Callback(hObject, eventdata, handles)
% hObject    handle to ucitaj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
naziv = get(handles.ime_datoteke, 'String');
naziv = [naziv '.txt']
global X
X = dlmread(naziv,'',2)
[duzina_trajanja,~] = max(X(:,1));
set(handles.duzina, 'String', duzina_trajanja);
[mini,~] = min(X(:,2));
[maxi,~] = max(X(:,2));
set(handles.minimalni, 'String', mini);
set(handles.maksimalni, 'String', maxi);

rr = [];
for i = 1:length(X(:,2))
    rr = [rr 60/X(i,2)];
end
x = rr; x(end) = [];
y = rr; y(1) = [];
L = length(x);

SD1C = sqrt((1/L) * sum(((x - y) - mean(x - y)).^2)/2);
SD2C = sqrt((1/L)*sum(((x+y) - mean(x+y)).^2)/2);
set(handles.s1, 'String', SD1C);
set(handles.s2, 'String', SD2C);
set(handles.s12, 'String', SD1C/SD2C);



% --- Executes on button press in vr_hr.
function vr_hr_Callback(hObject, eventdata, handles)
% hObject    handle to vr_hr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global X
[mn,~] = min(X(:,2));
[mx,~] = max(X(:,2));
figure;
plot(X(:,1),X(:,2), '-o');
title('Time representation of HR');
xlabel('Time [s]');
ylabel('Heart rate [BPM]');
hold on
grid on
axis([X(1,1) X(end,1) mn mx]);
% Hint: get(hObject,'Value') returns toggle state of vr_hr


% --- Executes on button press in bu_hr.
function bu_hr_Callback(hObject, eventdata, handles)
% hObject    handle to bu_hr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bu_hr


% --- Executes on button press in his_hr.
function his_hr_Callback(hObject, eventdata, handles)
% hObject    handle to his_hr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of his_hr


% --- Executes on button press in sc_hr.
function sc_hr_Callback(hObject, eventdata, handles)
% hObject    handle to sc_hr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sc_hr


% --- Executes on button press in delta_hr.
function delta_hr_Callback(hObject, eventdata, handles)
% hObject    handle to delta_hr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global X
y = [0]
for i = 1:length(X(:,2)) - 1
    y = [y X(i+1,2)-X(i,2)];
end
[mn,~] = min(y);
[mx,~] = max(y);
figure;
plot(X(:,1), y, '-o');
title('Delta HR');
xlabel('Time [s]');
ylabel('Delta HR');
hold on
grid on
axis([X(1,1) X(end,1) mn mx]);
% Hint: get(hObject,'Value') returns toggle state of delta_hr


% --- Executes during object creation, after setting all properties.
function duzina_CreateFcn(hObject, eventdata, handles)
% hObject    handle to duzina (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function minimalni_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minimalni (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function maksimalni_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maksimalni (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in fft.
function fft_Callback(hObject, eventdata, handles)
% hObject    handle to fft (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global X
Fs = 0.5;
T = 1/Fs;
L = length(X(:,2));
Y = fft(X(:,2));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
f = Fs * (0:(L/2))/L;
figure;
plot(f, P1);
title('Amplitude spectrum HR-a');
xlabel('f [Hz]');
ylabel('|HR(f)|');

% --- Executes on button press in spektar_snage.
function spektar_snage_Callback(hObject, eventdata, handles)
global X
Fs = 0.5;
N = length(X(:,2));
xfft = fft(X(:,2));
xfft = xfft(1:N/2+1);
psd = (1/(Fs*N)) * abs(xfft).^2;
psd(2:end-1) = 2*psd(2:end-1);
f = 0:Fs/length(X(:,2)): Fs/2;

figure;
plot(f, 10*log10(psd));
grid on

title('Periodogram using FFT')
xlabel('Frequency [Hz]')
ylabel('Power/Frequency [dB/Hz]')


% --- Executes on button press in rr_hr.
function rr_hr_Callback(hObject, eventdata, handles)
% hObject    handle to rr_hr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global X
[mn,~] = min(X(:,2));
[mx,~] = max(X(:,2));
rr = [];
for i = 1:length(X(:,2))
    rr = [rr 60/X(i,2)];
end
figure;
plot(X(:,1),rr, '-o');
title('Time representation of RR intervals');
xlabel('Time [s]');
ylabel('RR interval [s]');
hold on
grid on
%axis([X(1,1) X(end,1) 60/mn 60/mx]);

% --- Executes on button press in rr_beat.
function rr_beat_Callback(hObject, eventdata, handles)
% hObject    handle to rr_beat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global X
[mn,~] = min(X(:,2));
[mx,~] = max(X(:,2));
rr = [];
for i = 1:length(X(:,2))
    rr = [rr 60/X(i,2)];
end
figure;
plot(rr, '-o');
title('RR intervals by heart beats');
xlabel('Beat sample');
ylabel('RR interval [s]');
hold on
grid on
%axis([1 length(X(:,2)) 60/mn 60/mx]);

% --- Executes on button press in rr_hist.
function rr_hist_Callback(hObject, eventdata, handles)
% hObject    handle to rr_hist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global X
rr = [];
for i = 1:length(X(:,2))
    rr = [rr 60/X(i,2)];
end
figure;
hist(rr);
xlabel('Values of RR intervals');
ylabel('Frequency');
title('RR Histogram');

% --- Executes on button press in poincare_rr.
function poincare_rr_Callback(hObject, eventdata, handles)
% hObject    handle to poincare_rr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global X
x1 = [];
y1 = [];
rr = [];
for i = 1:length(X(:,2))
    rr = [rr 60/X(i,2)];
end
for i = 1:length(rr) - 1
    x1 = [x1 rr(i+1)];
    y1 = [y1 rr(i)];
end
figure;
scatter(x1,y1)
hold on
title('Poincare plot');
xlabel('R - R_{n} [s]');
ylabel('R - R_{n+1} [s]');
x = 0.4:0.01:1;
y = x;
plot(x,y,'r');
axis([0.4 1 0.4 1]);

% --- Executes on button press in o_projektu.
function o_projektu_Callback(hObject, eventdata, handles)
% hObject    handle to o_projektu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
h = msgbox({'University of Sarajevo' 'Faculty of Electrical Engineering' 'Department for automation control and electronics'...
    'Team leaders: Doc. dr Dušanka Boškoviæ, Doc. dr Almir Badnjeviæ'...
    'Students: Aèkar Haris, Kico Iris, Tahirbegoviæ Anel'}','About project!');


% --- Executes on button press in deskriptor.
function deskriptor_Callback(hObject, eventdata, handles)
% hObject    handle to deskriptor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global X;
rr = [];
for i = 1:length(X(:,2))
    rr = [rr 60/X(i,2)];
end
x = rr; x(end) = [];
y = rr; y(1) = [];
L = length(x);

SD1C = sqrt((1/L) * sum(((x - y) - mean(x - y)).^2)/2)
SD2C = sqrt((1/L)*sum(((x+y) - mean(x+y)).^2)/2)
