function varargout = Majda_Zavrsni(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Majda_Zavrsni_OpeningFcn, ...
                   'gui_OutputFcn',  @Majda_Zavrsni_OutputFcn, ...
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


% --- Executes just before Majda_Zavrsni is made visible.
function Majda_Zavrsni_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Majda_Zavrsni (see VARARGIN)

% Choose default command line output for Majda_Zavrsni
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
global data;
global time_cut;
global signali;
global rezani_signali;
global studija;
% UIWAIT makes Majda_Zavrsni wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Majda_Zavrsni_OutputFcn(hObject, eventdata, handles) 
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


% --- Executes on button press in Load.
function Load_Callback(hObject, eventdata, handles)
% hObject    handle to Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
d = uigetdir;
%c=get(handles.popup,'String');
c=cell(0);
global data;
global odr;
global signali;
global signali_rr;
global time_cut;
global typerr;

typerr= get(handles.sigrr,'Value');
info = dir(fullfile(d, '*.txt'));
data = cell(max(size(info)),8);
odr = cell(max(size(info)),1);
signali=cell(max(size(info)),2);
signali_rr=cell(max(size(info)),1);
time_cut=cell(max(size(info)),2);
%leg = cell(max(size(info)),1);
%size(leg)
%info=dir(filename);
%info.name
%size(info.name)
%x=ewwwwwww
max_x= -1;
max_y= -1;
min_y= 100000;
min_x = 100000;
for qw=1:max(size(info))
    c{max(size(c))+1}=strcat('Signal:',num2str(qw));
    X = dlmread(fullfile(d,info(qw).name),'',2);
    if(typerr ~= 1)
        signali{qw,1}=X(:,1);
        signali{qw,2}=X(:,2);
        if(min(X(:,1)) < min_x)
            min_x=min(X(:,1));
        end
        if(max(X(:,1)) > max_x)
            max_x=max(X(:,1));
        end
        if(min_y > min(X(:,2)))
            min_y = min(X(:,2));
        end
        if(max(X(:,2)) > max_y)
            max_y = max(X(:,2));
        end
        rr = [];
        for i = 1:length(X(:,2))%pretvara HR u RR
            rr = [rr 60/X(i,2)];
        end
        x = rr; x(end) = [];
        y = rr; y(1) = [];
        L = length(x);
        signali_rr{qw} = rr;
        SD1C = sqrt((1/L) * sum(((x - y) - mean(x - y)).^2)/2);
        SD2C = sqrt((1/L)*sum(((x+y) - mean(x+y)).^2)/2);
        data{qw,1} = max(X(:,1)); %% trajanje signala
        data{qw,2} = max(X(:,2)); %% max amplituda
        data{qw,3} = min(X(:,2)); %% min amplituda
        time_cut{qw,1}=min(X(:,1));
        time_cut{qw,2}=max(X(:,1));
        data{qw,4} = SD1C; %% sd1c
        data{qw,5} = SD2C; %% sd2c
        data{qw,6} = SD1C/SD2C; %% sd1c/sd2c
        data{qw,7} = mean(X(:,2)); %% average
        data{qw,8} = min(X(:,1)); %% average
    else
        signali{qw,1}=X(:,1);
        signali_rr{qw}=X(:,2);
        hr = [];
        for i = 1:length(X(:,2))
            hr = [hr 60/X(i,2)];
        end
        if(min(X(:,1)) < min_x)
            min_x=min(X(:,1));
        end
        if(max(X(:,1)) > max_x)
            max_x=max(X(:,1));
        end
        if(min_y > min(X(:,2)))
            min_y = min(X(:,2));
        end
        if(max(X(:,2)) > max_y)
            max_y = max(X(:,2));
        end


        x = X(:,2); x(end) = [];
        y = X(:,2); y(1) = [];
        L = length(x);
        signali{qw,2} = hr;
        SD1C = sqrt((1/L) * sum(((x - y) - mean(x - y)).^2)/2);
        SD2C = sqrt((1/L)*sum(((x+y) - mean(x+y)).^2)/2);
        data{qw,1} = max(X(:,1)); %% trajanje signala
        data{qw,2} = max(X(:,2)); %% max amplituda
        data{qw,3} = min(X(:,2)); %% min amplituda
        time_cut{qw,1}=min(X(:,1));
        time_cut{qw,2}=max(X(:,1));
        data{qw,4} = SD1C; %% sd1c
        data{qw,5} = SD2C; %% sd2c
        data{qw,6} = SD1C/SD2C; %% sd1c/sd2c
        data{qw,7} = mean(X(:,2)); %% average
        data{qw,8} = min(X(:,1)); %% average
    end
    strcat('Signal:',num2str(qw))
    %leg{qw}=strcat('Signal:',num2str(qw));
end
set(handles.popup,'String',c);
set(handles.sl1,'Min',min_x);
set(handles.sl2,'Max',max_x);
set(handles.sl1,'Min',min_x);
set(handles.sl2,'Max',max_x);
set(handles.minsl,'String',num2str(min_x));
set(handles.maxsl,'String',num2str(max_x));
axis(handles.fig1,[min_x-1 max_x+1 min(X(:,2)) - 1 max(X(:,2)) + 1]);
%legend(handles.fig1,leg);
% --- Executes on button press in vr_hr.
function vr_hr_Callback(hObject, eventdata, handles)
% hObject    handle to vr_hr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global signali;
global signali_rr;
global time_cut;
global typerr;
if(typerr ~= 1)
    c=get(handles.popup,'Value');
    X1=signali{c,1};
    X2=signali{c,2};
    k1 = find(X1 <= time_cut{c,1});k1=k1(end);
    k2 = find(X1 >= time_cut{c,2});k2=k2(1);
    X1=X1(k1:k2);
    X2=X2(k1:k2);
    X=[X1 X2];
    [mn,~] = min(X(:,2));
    [mx,~] = max(X(:,2));
    figure;
    plot(X(:,1),X(:,2), 'Marker','o','MarkerFaceColor','r','Color','b','LineWidth',1.7);
    title('Time representation of HR');
    xlabel('Time [s]');
    ylabel('Heart rate [BPM]');
    hold on
    grid on
    axis([X(1,1) X(end,1) mn mx]);
else
    c=get(handles.popup,'Value');
    X1=signali{c,1};
    X2=signali_rr{c};
    k1 = find(X1 <= time_cut{c,1});k1=k1(end);
    k2 = find(X1 >= time_cut{c,2});k2=k2(1);
    X1=X1(k1:k2);
    X2=X2(k1:k2);
    
    hr=[];
    for i = 1:length(X2)
        hr = [hr 60/X2(i)];
    end
    X2 = hr;
    k1 = find(X1 <= time_cut{c,1});k1=k1(end);
    k2 = find(X1 >= time_cut{c,2});k2=k2(1);
    X1=X1(k1:k2);
    X2=X2(k1:k2);
    X=[X1 X2'];
    [mn,~] = min(X(:,2));
    [mx,~] = max(X(:,2));
    figure;
    plot(X(:,1),X(:,2), 'Marker','o','MarkerFaceColor','r','Color','b','LineWidth',1.7);
    title('Time representation of HR');
    xlabel('Time [s]');
    ylabel('Heart rate [BPM]');
    hold on
    grid on
    axis([X(1,1) X(end,1) mn mx]);
end
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
global time_cut;
global signali;
global signali_rr;
global typerr;
c=get(handles.popup,'Value');

if(typerr ~= 1)
    X1=signali{c,1};
    X2=signali{c,2};
    k1 = find(X1 <= time_cut{c,1});k1=k1(end);
    k2 = find(X1 >= time_cut{c,2});k2=k2(1);
    X1=X1(k1:k2);
    X2=X2(k1:k2);
    X=[X1 X2];
    y = [0];
    for i = 1:length(X(:,2)) - 1
        y = [y X(i+1,2)-X(i,2)];
    end
    [mn,~] = min(y);
    [mx,~] = max(y);
    figure;
    plot(X(:,1), y, 'Color','b', 'Marker','o','MarkerFaceColor','r','LineWidth',1.7);
    title('Delta HR');
    xlabel('Time [s]');
    ylabel('Delta HR [bpm]');
    hold on
    grid on
    axis([X(1,1) X(end,1) mn mx]);
else
    hr=[];
    X1=signali{c,1};
    X2=signali_rr{c};
    k1 = find(X1 <= time_cut{c,1});k1=k1(end);
    k2 = find(X1 >= time_cut{c,2});k2=k2(1);
    X1=X1(k1:k2);
    X2=X2(k1:k2);

    for i = 1:length(X2)
        hr = [hr 60/X2(i)];
    end
    X2 = hr;
    X=[X1 X2'];
    y = [0];
    for i = 1:length(X(:,2)) - 1
        y = [y X(i+1,2)-X(i,2)];
    end
    [mn,~] = min(y);
    [mx,~] = max(y);
    figure;
    plot(X(:,1), y, 'Color','b', 'Marker','o','MarkerFaceColor','r','LineWidth',1.7);
    title('Delta HR');
    xlabel('Time [s]');
    ylabel('Delta HR [bpm]');
    hold on
    grid on
    axis([X(1,1) X(end,1) mn mx]);
end
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
global time_cut;
global signali;
c=get(handles.popup,'Value');
X1=signali{c,1};
X2=signali{c,2};
k1 = find(X1 <= time_cut{c,1});k1=k1(end);
k2 = find(X1 >= time_cut{c,2});k2=k2(1);
X1=X1(k1:k2);
X2=X2(k1:k2);
X=[X1 X2];

Fs = 0.5;
T = 1/Fs;
L = length(X(:,2));
Y = fft(X(:,2));
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
f = Fs * (0:(L/2))/L;
figure;
plot(f, P1);
title('Amplituda spektra HR-a');
xlabel('f [Hz]');
ylabel('|HR(f)|');

% --- Executes on button press in spektar_snage.
function spektar_snage_Callback(hObject, eventdata, handles)
global time_cut;
global signali;
global signali_rr;
global typerr;
c=get(handles.popup,'Value');
if(typerr ~= 1)
    X1=signali{c,1};
    X2=signali{c,2};
    k1 = find(X1 <= time_cut{c,1});k1=k1(end);
    k2 = find(X1 >= time_cut{c,2});k2=k2(1);
    X1=X1(k1:k2);
    X2=X2(k1:k2);
    X=[X1 X2];

    Fs = 0.5;
    N = length(X(:,2));
    xfft = fft(X(:,2));
else
    X1=signali{c,1};
    X2=signali_rr{c};
    k1 = find(X1 <= time_cut{c,1});k1=k1(end);
    k2 = find(X1 >= time_cut{c,2});k2=k2(1);
    X1=X1(k1:k2);
    X2=X2(k1:k2);
    X=[X1 X2];
    Fs = 0.5;
    N = length(X(:,2));
    xfft = fft(X(:,2));
end
xfft = xfft(1:N/2+1);
psd = (1/(Fs*N)) * abs(xfft).^2;
psd(2:end-1) = 2*psd(2:end-1);
f = 0:Fs/length(X(:,2)): Fs/2;

figure;
area(f, 10*log10(psd));
grid on

title('Power spectrum')
xlabel('Frequency [Hz]')
ylabel('Power/Frequency [dB/Hz]')


% --- Executes on button press in rr_hr.
function rr_hr_Callback(hObject, eventdata, handles)
% hObject    handle to rr_hr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global time_cut;
global signali;
global signali_rr;
global typerr;
c=get(handles.popup,'Value');

if(typerr ~= 1)
    X1=signali{c,1};
    X2=signali{c,2};
    k1 = find(X1 <= time_cut{c,1});k1=k1(end);  
    k2 = find(X1 >= time_cut{c,2});k2=k2(1);
    X1=X1(k1:k2);
    X2=X2(k1:k2);
    X=[X1 X2];
    [mn,~] = min(X(:,2));
    [mx,~] = max(X(:,2));
rr = [];
    for i = 1:length(X(:,2))
        rr = [rr 60/X(i,2)];
    end
else
    rr=signali_rr{c};
     X1=signali{c,1};
    k1 = find(X1 <= time_cut{c,1});k1=k1(end);  
    k2 = find(X1 >= time_cut{c,2});k2=k2(1);
    X1=X1(k1:k2);
    rr=rr(k1:k2);
end
figure;
plot(X1,rr, 'Color','b', 'Marker','o','MarkerFaceColor','r','LineWidth',1.7);
title('Time representation of RR');
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
global time_cut;
global signali;
global signali_rr;
global typerr;
c=get(handles.popup,'Value');
if(typerr ~= 1)
    X1=signali{c,1};
    X2=signali{c,2};
    
    k1 = find(X1 <= time_cut{c,1});k1=k1(end);
    k2 = find(X1 >= time_cut{c,2});k2=k2(1);
    X1=X1(k1:k2);
    X2=X2(k1:k2);
  

    [mn,~] = min(X2);
    [mx,~] = max(X2);
    rr = [];
     for i = 1:length(X2)
         rr = [rr 60/X2(i)];
     end
else
    X1=signali{c,1};
    X2=signali_rr{c};
    k1 = find(X1 <= time_cut{c,1});k1=k1(end);
    k2 = find(X1 >= time_cut{c,2});k2=k2(1);
    X1=X1(k1:k2);
    X2=X2(k1:k2);
    rr = X2;
end
figure;
plot(rr, 'Color','b', 'Marker','o','MarkerFaceColor','r','LineWidth',1.7);
title('RR intervals by heart beats');
xlabel('Heart beat');
ylabel('RR interval [s]');
hold on
grid on
%axis([1 length(X(:,2)) 60/mn 60/mx]);

% --- Executes on button press in rr_hist.
function rr_hist_Callback(hObject, eventdata, handles)
% hObject    handle to rr_hist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global time_cut;
global signali;
global signali_rr;
global typerr;
c=get(handles.popup,'Value');
if(typerr  ~= 1)
        X1=signali{c,1};
        X2=signali{c,2};
        k1 = find(X1 <= time_cut{c,1});k1=k1(end);
        k2 = find(X1 >= time_cut{c,2});k2=k2(1);
        X1=X1(k1:k2);
        X2=X2(k1:k2);
        X=[X1 X2];

        rr = [];
        for i = 1:length(X(:,2))
            rr = [rr 60/X(i,2)];
        end
else
        X1=signali{c,1};
        X2=signali_rr{c};
        k1 = find(X1 <= time_cut{c,1});k1=k1(end);
        k2 = find(X1 >= time_cut{c,2});k2=k2(1);
        X1=X1(k1:k2);
        X2=X2(k1:k2);
        rr = X2;

end
figure;
hist(rr);
xlabel('RR [s]');
ylabel('Number of RR intervals');
title('RR Histogram');

% --- Executes on button press in poincare_rr.
function poincare_rr_Callback(hObject, eventdata, handles)
% hObject    handle to poincare_rr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global time_cut;
global signali;
global signali_rr;
global typerr;
global data;
c=get(handles.popup,'Value');
if(typerr ~= 1)
    X1=signali{c,1};
    X2=signali{c,2};
    k1 = find(X1 <= time_cut{c,1});k1=k1(end);
    k2 = find(X1 >= time_cut{c,2});k2=k2(1);
    X1=X1(k1:k2);
    X2=X2(k1:k2);
    X=[X1 X2];
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
else
    X1=signali{c,1};
    X2=signali_rr{c};
    k1 = find(X1 <= time_cut{c,1});k1=k1(end);
    k2 = find(X1 >= time_cut{c,2});k2=k2(1);
    X1=X1(k1:k2);
    X2=X2(k1:k2);
    X=[X1 X2];
    x1 = [];
    y1 = [];
    rr = X2;
    for i = 1:length(rr) - 1
        x1 = [x1 rr(i+1)];
        y1 = [y1 rr(i)];
    end
end
figure;
sd1 = data{c,4};
sd2 = data{c,5};
% xc = mean(x1);
% yc = mean(y1);
% fi = 45*180/pi;
% th = 0 : 0.01 : 2*pi; th = th + fi;
% xell = sd1*cos(th) + xc;
% yell = sd2*sin(th) + yc;
scatter(x1,y1)
hold on
title('Poincare plot');
xlabel('R - R_{n} [s]');
ylabel('R - R_{n+1} [s]');
x = 0.4:0.01:1;
y = x;
plot(x,y,'r','LineWidth',1.7);

axis([0.4 1 0.4 1]);
% hold on;
% plot(xell,yell,'k','LineWidth',2);
% hold on;
% line([xc, xc], [yc - sd2, yc + sd2], ...
% 	'LineWidth', 4, 'Color', [1,0,0]);
% line([xc - sd1], xc + sd1, [yc, yc], ...
% 	'LineWidth', 4, 'Color', [1,0,0]);
% --- Executes on button press in o_projektu.
function o_projektu_Callback(hObject, eventdata, handles)
% hObject    handle to o_projektu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



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


% --- Executes on slider movement.
function sl1_Callback(hObject, eventdata, handles)
% hObject    handle to sl1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global time_cut;
global signali;
global signali_rr;
global data;
global typerr;
if(typerr ~= 1)
    c=get(handles.popup,'Value');
    time_cut{c,1}= round(get(handles.sl1,'Value'));
    set(handles.sl1,'Value',time_cut{c,1});
    set(handles.sl2,'Value',time_cut{c,2});
    set(handles.minsl,'String',num2str(time_cut{c,1}));
    plot(handles.fig1,signali{c,1},signali{c,2}, 'Color','b','Marker','o','MarkerFaceColor','r','LineWidth',1.7); title(handles.fig1,'HR tachogram');
    xlabel(handles.fig1,'Time [s]');
    ylabel(handles.fig1,'Heart beat [BPM]');
    grid(handles.fig1,'on');
    legend(handles.fig1,strcat('Signal:',num2str(c)));
    xlim(handles.fig1,[time_cut{c,1} time_cut{c,2}]);

    rr = [];
    X1=signali{c,1};
    X2=signali{c,2};
    k1 = find(X1 <= time_cut{c,1});k1=k1(end);
    k2 = find(X1 >= time_cut{c,2});k2=k2(1);
    X1=X1(k1:k2);
    X2=X2(k1:k2);
    X=[X1 X2];
    for i = 1:length(X(:,2))
        rr = [rr 60/X(i,2)];
    end
    x = rr; x(end) = [];
    y = rr; y(1) = [];
    L = length(x);

    SD1C = sqrt((1/L) * sum(((x - y) - mean(x - y)).^2)/2);
    SD2C = sqrt((1/L)*sum(((x+y) - mean(x+y)).^2)/2);
    ylim(handles.fig1,[min(X1) max(X2)]);
    %data{c,1} = max(X(:,1)); %% trajanje signala
    %data{c,2} = max(X(:,2)); %% max amplituda
    %data{c,3} = min(X(:,2)); %% min amplituda
    data{c,4} = SD1C; %% sd1c
    data{c,5} = SD2C; %% sd2c
    data{c,6} = SD1C/SD2C; %% sd1c/sd2c
    data{c,7} = mean(X(:,2)); %% average
    set(handles.duzina, 'String', max(X(:,1))-min(X(:,1)));
    set(handles.minimalni, 'String', min(X(:,2)));
    set(handles.maksimalni, 'String', max(X(:,2)));
    set(handles.edit18, 'String', data{c,4});
    set(handles.text34, 'String', data{c,5});
    set(handles.text35, 'String', data{c,6});
    set(handles.text9,'String','[s]');
    set(handles.text7,'String','[bpm]');
    set(handles.text8,'String','[bpm]');
else
   c=get(handles.popup,'Value');
    time_cut{c,1}= round(get(handles.sl1,'Value'));
    set(handles.sl1,'Value',time_cut{c,1});
    set(handles.sl2,'Value',time_cut{c,2});
    set(handles.minsl,'String',num2str(time_cut{c,1}));
    plot(handles.fig1,signali{c,1},signali_rr{c}, 'Color','b','Marker','o','MarkerFaceColor','r','LineWidth',1.7); title(handles.fig1,'RR tachogram');
    xlabel(handles.fig1,'Time [s]');
    ylabel(handles.fig1,'RR intervals[s]');
    grid(handles.fig1,'on');
    legend(handles.fig1,strcat('Signal:',num2str(c)));
    xlim(handles.fig1,[time_cut{c,1} time_cut{c,2}]);

    rr = [];
    X1=signali{c,1};
    X2=signali_rr{c,1};
    k1 = find(X1 <= time_cut{c,1});k1=k1(end);
    k2 = find(X1 >= time_cut{c,2});k2=k2(1);
    X1=X1(k1:k2);
    X2=X2(k1:k2);
    rr=X2;
    x = rr; x(end) = [];
    y = rr; y(1) = [];
    L = length(x);

    SD1C = sqrt((1/L) * sum(((x - y) - mean(x - y)).^2)/2);
    SD2C = sqrt((1/L)*sum(((x+y) - mean(x+y)).^2)/2);
           
    xlim(handles.fig1,[X1(1) X1(end)]);
    ylim(handles.fig1,[min(X2)-1 max(X2)+1]);

    %data{c,1} = max(X(:,1)); %% trajanje signala
    %data{c,2} = max(X(:,2)); %% max amplituda
    %data{c,3} = min(X(:,2)); %% min amplituda
    data{c,4} = SD1C; %% sd1c
    data{c,5} = SD2C; %% sd2c
    data{c,6} = SD1C/SD2C; %% sd1c/sd2c
    data{c,7} = mean(X2); %% average
    set(handles.duzina, 'String', max(X1)-min(X1));
    set(handles.minimalni, 'String', min(X2));
    set(handles.maksimalni, 'String', max(X2));
    %set(handles.edit15, 'String', data{c,4});
    %set(handles.edit16, 'String', data{c,5});
    %set(handles.edit17, 'String', data{c,6});
    
    set(handles.text9,'String','[s]');
    set(handles.edit18, 'String', data{c,4});
    set(handles.text34, 'String', data{c,4});
    set(handles.text35, 'String', data{c,4});
    set(handles.text7,'String','[s]');
    set(handles.text8,'String','[s]');
end
% --- Executes during object creation, after setting all properties.
function sl1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sl1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sl2_Callback(hObject, eventdata, handles)
% hObject    handle to sl2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

global time_cut;
global signali;
global signali_rr;
global data;
global typerr;
if(typerr ~= 1)
    c=get(handles.popup,'Value');
    time_cut{c,2}= round(get(handles.sl2,'Value'));
    set(handles.maxsl,'String',num2str(time_cut{c,2}));
    set(handles.sl1,'Value',time_cut{c,1});
    set(handles.sl2,'Value',time_cut{c,2});
    plot(handles.fig1,signali{c,1},signali{c,2}, 'Color','b','Marker','o','MarkerFaceColor','r','LineWidth',1.7); title(handles.fig1,'HR tachogram');
    xlabel(handles.fig1,'Time [s]');
    ylabel(handles.fig1,'Heart rate [BPM]');
    grid(handles.fig1,'on');
    legend(handles.fig1,strcat('Signal:',num2str(c)));
    X1=signali{c,1};
    X2=signali{c,2};
    k1 = find(X1 <= time_cut{c,1});k1=k1(end);
    k2 = find(X1 >= time_cut{c,2});k2=k2(1);
    X1=X1(k1:k2);
    X2=X2(k1:k2);
    X=[X1 X2];
    xlim(handles.fig1,[X1(1) X1(end)]);
    ylim(handles.fig1,[min(X2) max(X2)]);
    rr = [];

    for i = 1:length(X(:,2))
        rr = [rr 60/X(i,2)];
    end
    x = rr; x(end) = [];
    y = rr; y(1) = [];
    L = length(x);

    SD1C = sqrt((1/L) * sum(((x - y) - mean(x - y)).^2)/2);
    SD2C = sqrt((1/L)*sum(((x+y) - mean(x+y)).^2)/2);

    %data{c,1} = max(X(:,1)); %% trajanje signala
    %data{c,2} = max(X(:,2)); %% max amplituda
    %data{c,3} = min(X(:,2)); %% min amplituda
    data{c,4} = SD1C; %% sd1c
    data{c,5} = SD2C; %% sd2c
    data{c,6} = SD1C/SD2C; %% sd1c/sd2c
    data{c,7} = mean(X(:,2)); %% average
    set(handles.duzina, 'String', max(X(:,1))-min(X(:,1)));
    set(handles.minimalni, 'String', min(X(:,2)));
    set(handles.maksimalni, 'String', max(X(:,2)));
    set(handles.edit18, 'String', data{c,4});
    set(handles.text34, 'String', data{c,5});
    set(handles.text35, 'String', data{c,6});
    set(handles.text9, 'String','[s]');
    set(handles.text7,'String','[bpm]');
    set(handles.text8,'String','[bpm]');
else
        c=get(handles.popup,'Value');
    time_cut{c,2}= round(get(handles.sl2,'Value'));
    set(handles.maxsl,'String',num2str(time_cut{c,2}));
    set(handles.sl1,'Value',time_cut{c,1});
    set(handles.sl2,'Value',time_cut{c,2});
    plot(handles.fig1,signali{c,1},signali_rr{c}, 'Color','b','Marker','o','MarkerFaceColor','r','LineWidth',1.7); title(handles.fig1,'RR tachogram');
    xlabel(handles.fig1,'Time [s]');
    ylabel(handles.fig1,'RR interval [s]');
    grid(handles.fig1,'on');
    legend(handles.fig1,strcat('Signal:',num2str(c)));
    X1=signali{c,1};
    X2=signali_rr{c,1};
    k1 = find(X1 <= time_cut{c,1});k1=k1(end);
    k2 = find(X1 >= time_cut{c,2});k2=k2(1);
    X1=X1(k1:k2);
    
    X2=X2(k1:k2);
 
    rr = [];
    for i = 1:length(X2)
        rr = [rr 60/X2(i)];
    end

    rr=X2;
 
    xlim(handles.fig1,[X1(1) X1(end)]);
    ylim(handles.fig1,[min(X2)-1 max(X2)+1]);

    x = rr; x(end) = [];
    y = rr; y(1) = [];
    L = length(x);

    SD1C = sqrt((1/L) * sum(((x - y) - mean(x - y)).^2)/2);
    SD2C = sqrt((1/L)*sum(((x+y) - mean(x+y)).^2)/2);

    %data{c,1} = max(X(:,1)); %% trajanje signala
    %data{c,2} = max(X(:,2)); %% max amplituda
    %data{c,3} = min(X(:,2)); %% min amplituda
    data{c,4} = SD1C; %% sd1c
    data{c,5} = SD2C; %% sd2c
    data{c,6} = SD1C/SD2C; %% sd1c/sd2c
    data{c,7} = mean(X2); %% average
    set(handles.duzina, 'String', max(X1)-min(X1));
    set(handles.minimalni, 'String', min(X2));
    set(handles.maksimalni, 'String', max(X2));
    set(handles.edit18, 'String', data{c,4});
    set(handles.text34, 'String', data{c,5});
    set(handles.text35, 'String', data{c,6});
    set(handles.text9,'String','[s]');
    set(handles.text7,'String','[s]');
    set(handles.text8,'String','[s]');
end
% --- Executes during object creation, after setting all properties.
function sl2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sl2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function minsl_Callback(hObject, eventdata, handles)
% hObject    handle to minsl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minsl as text
%        str2double(get(hObject,'String')) returns contents of minsl as a double


% --- Executes during object creation, after setting all properties.

global time_cut;
global signali;
global signali_rr;
global data;
global typerr;
if(typerr ~= 1)
    c=get(handles.popup,'Value');
    set(handles.sl1,'Value',round(str2double(get(handles.minsl,'String'))));
    set(handles.minsl,'String',round(str2double(get(handles.minsl,'String'))));
    time_cut{c,1}=round(str2double(get(handles.minsl,'String')));
    plot(handles.fig1,signali{c,1},signali{c,2}, 'Color','b','Marker','o','MarkerFaceColor','r','LineWidth',1.7);
    title(handles.fig1,'HR tachogram');
    xlabel(handles.fig1,'Time [s]');
    ylabel(handles.fig1,'Heart rate [BPM]');
    grid(handles.fig1,'on');
    legend(handles.fig1,strcat('Signal:',num2str(c)));
    X1=signali{c,1};
    X2=signali{c,2};
    k1 = find(X1 <= time_cut{c,1});
    k1=k1(end);
    k2 = find(X1 >= time_cut{c,2});k2=k2(1);
    X1=X1(k1:k2);
    X2=X2(k1:k2);
    X=[X1 X2];

    rr = [];
    for i = 1:length(X2)
        rr = [rr 60/X2];
    end
    xlim(handles.fig1,[X1(1) X1(end)]);
    ylim(handles.fig1,[min(X2) max(X2)]);

    x = rr; x(end) = [];
    y = rr; y(1) = [];
    L = length(x);

    SD1C = sqrt((1/L) * sum(((x - y) - mean(x - y)).^2)/2);
    SD2C = sqrt((1/L)*sum(((x+y) - mean(x+y)).^2)/2);

    %data{c,1} = max(X(:,1)); %% trajanje signala
    %data{c,2} = max(X(:,2)); %% max amplituda
    %data{c,3} = min(X(:,2)); %% min amplituda
    data{c,4} = SD1C; %% sd1c
    data{c,5} = SD2C; %% sd2c
    data{c,6} = SD1C/SD2C; %% sd1c/sd2c
    data{c,7} = mean(X(:,2)); %% average
    set(handles.duzina, 'String', max(X(:,1))-min(X(:,1)));
    set(handles.minimalni, 'String', min(X(:,2)));
    set(handles.maksimalni, 'String', max(X(:,2)));
    set(handles.edit18, 'String', data{c,4});
    set(handles.text34, 'String', data{c,5});
    set(handles.text35, 'String', data{c,6});
    set(handles.text9,'String','[s]');
    set(handles.text7,'String','[bpm]');
    set(handles.text8,'String','[bpm]');
    %set(handles.sl2,'Value',get(str2num(handles.maxsl,'String')));
else
    c=get(handles.popup,'Value');
    set(handles.sl1,'Value',round(str2double(get(handles.minsl,'String'))));
    set(handles.minsl,'String',round(str2double(get(handles.minsl,'String'))));
    time_cut{c,1}=round(str2double(get(handles.minsl,'String')));
    plot(handles.fig1,signali{c,1},signali_rr{c}, 'Color','b','Marker','o','MarkerFaceColor','r','LineWidth',1.7); 
    title(handles.fig1,'RR tachogram');
    xlabel(handles.fig1,'Time [s]');
    ylabel(handles.fig1,'RR interval [s]');
    grid(handles.fig1,'on');
    legend(handles.fig1,strcat('Signal:',num2str(c)));
    X1=signali{c,1};
    X2=signali_rr{c};
    k1 = find(X1 <= time_cut{c,1});
    k1=k1(end);
    k2 = find(X1 >= time_cut{c,2});k2=k2(1);
    X1=X1(k1:k2);
    X2=X2(k1:k2);



    rr = X2;
    xlim(handles.fig1,[X1(1) X1(end)]);
    ylim(handles.fig1,[min(X2)-1 max(X2)+1]);
    x = rr; x(end) = [];
    y = rr; y(1) = [];
    L = length(x);

    SD1C = sqrt((1/L) * sum(((x - y) - mean(x - y)).^2)/2);
    SD2C = sqrt((1/L)*sum(((x+y) - mean(x+y)).^2)/2);

    %data{c,1} = max(X(:,1)); %% trajanje signala
    %data{c,2} = max(X(:,2)); %% max amplituda
    %data{c,3} = min(X(:,2)); %% min amplituda
    data{c,4} = SD1C; %% sd1c
    data{c,5} = SD2C; %% sd2c
    data{c,6} = SD1C/SD2C; %% sd1c/sd2c
    data{c,7} = mean(X2); %% average
    set(handles.duzina, 'String', max(X1)-min(X1));
    set(handles.minimalni, 'String', min(X2));
    set(handles.maksimalni, 'String', max(X2));
    set(handles.edit18, 'String', data{c,4});
    set(handles.text34, 'String', data{c,5});
    set(handles.text35, 'String', data{c,6});
    set(handles.text9,'String','[s]');
    set(handles.text7,'String','[s]');
    set(handles.text8,'String','[s]');
    %set(handles.sl2,'Value',get(str2num(handles.maxsl,'String')));
end
function minsl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minsl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function maxsl_Callback(hObject, eventdata, handles)
% hObject    handle to maxsl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxsl as text
%        str2double(get(hObject,'String')) returns contents of maxsl as a double


global time_cut;
global signali;
global signali_rr;
global data;
global typerr;
if(typerr ~= 1)
    c=get(handles.popup,'Value');
    set(handles.sl2,'Value',round(str2double(get(handles.maxsl,'String'))));
    time_cut{c,2}=round(str2double(get(handles.maxsl,'String')));
    set(handles.maxsl,'String',round(str2double(get(handles.maxsl,'String'))));
    plot(handles.fig1,signali{c,1},signali{c,2}, 'Color','b','Marker','o','MarkerFaceColor','r','LineWidth',1.7);
    title(handles.fig1,'HR tachogram');
    xlabel(handles.fig1,'Time [s]');
    ylabel(handles.fig1,'Heart rate [BPM]');
    grid(handles.fig1,'on');
    legend(handles.fig1,strcat('Signal:',num2str(c)));

    X1=signali{c,1};
    X2=signali{c,2};
    k1 = find(X1 <= time_cut{c,1});k1=k1(end);
    k2 = find(X1 >= time_cut{c,2});k2=k2(1);
    X1=X1(k1:k2);
    X2=X2(k1:k2);
    xlim(handles.fig1,[X1(1) X1(end)]);
    ylim(handles.fig1,[data{c,3} data{c,2}]);
    X=[X1 X2];
    rr=[];
    for i = 1:length(X(:,2))
        rr = [rr 60/X(i,2)];
    end
    x = rr; x(end) = [];
    y = rr; y(1) = [];
    L = length(x);

    SD1C = sqrt((1/L) * sum(((x - y) - mean(x - y)).^2)/2);
    SD2C = sqrt((1/L)*sum(((x+y) - mean(x+y)).^2)/2);

    %data{c,1} = max(X(:,1)); %% trajanje signala
    %data{c,2} = max(X(:,2)); %% max amplituda
    %data{c,3} = min(X(:,2)); %% min amplituda
    data{c,4} = SD1C; %% sd1c
    data{c,5} = SD2C; %% sd2c
    data{c,6} = SD1C/SD2C; %% sd1c/sd2c
    data{c,7} = mean(X(:,2)); %% average
    set(handles.duzina, 'String', max(X(:,1))-min(X(:,1)));
    set(handles.minimalni, 'String', min(X(:,2)));
    set(handles.maksimalni, 'String', max(X(:,2)));
    set(handles.edit18, 'String', data{c,4});
    set(handles.text34, 'String', data{c,5});
    set(handles.text35, 'String', data{c,6});
    set(handles.text9,'String','[s]');
    set(handles.text7,'String','[bpm]');
    set(handles.text8,'String','[bpm]');
else
     c=get(handles.popup,'Value');
    set(handles.sl2,'Value',round(str2double(get(handles.maxsl,'String'))));
    time_cut{c,2}=round(str2double(get(handles.maxsl,'String')));
    set(handles.maxsl,'String',round(str2double(get(handles.maxsl,'String'))));
    plot(handles.fig1,signali{c,1},signali_rr{c}, 'Color','b','Marker','o','MarkerFaceColor','r','LineWidth',1.7); 
    title(handles.fig1,'RR tachogram');
    xlabel(handles.fig1,'Time [s]');
    ylabel(handles.fig1,'RR interval [s]');
    grid(handles.fig1,'on');
    legend(handles.fig1,strcat('Signal:',num2str(c)));

    X1=signali{c,1};
    X2=signali_rr{c,1};
    k1 = find(X1 <= time_cut{c,1});k1=k1(end);
    k2 = find(X1 >= time_cut{c,2});k2=k2(1);
    X1 = X1(k1:k2);
    rr = X2;
    xlim(handles.fig1,[X1(1) X1(end)]);
    ylim(handles.fig1,[min(X2)-1 max(X2)+1]);

    
    x = rr; x(end) = [];
    y = rr; y(1) = [];
    L = length(x);

    SD1C = sqrt((1/L) * sum(((x - y) - mean(x - y)).^2)/2);
    SD2C = sqrt((1/L)*sum(((x+y) - mean(x+y)).^2)/2);

    %data{c,1} = max(X(:,1)); %% trajanje signala
    %data{c,2} = max(X(:,2)); %% max amplituda
    %data{c,3} = min(X(:,2)); %% min amplituda
    data{c,4} = SD1C; %% sd1c
    data{c,5} = SD2C; %% sd2c
    data{c,6} = SD1C/SD2C; %% sd1c/sd2c
    data{c,7} = mean(X2); %% average
    set(handles.duzina, 'String', max(X1)-min(X1));
    set(handles.minimalni, 'String', min(X2));
    set(handles.maksimalni, 'String', max(X2));
    set(handles.edit18, 'String', data{c,4});
    set(handles.text34, 'String', data{c,5});
    set(handles.text35, 'String', data{c,6});
    set(handles.text9,'String','[s]');
    set(handles.text7,'String','[s]');
    set(handles.text8,'String','[s]');
end
% --- Executes during object creation, after setting all properties.
function maxsl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxsl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgdatroundColor','white');
end


% --- Executes on selection change in popup.
function popup_Callback(hObject, eventdata, handles)
c=get(handles.popup,'Value');
global data;
global time_cut;
global signali_rr;
global typerr;
global signali;
if(typerr ~= 1)
    plot(handles.fig1,signali{c,1},signali{c,2}, 'Color','b','Marker','o','MarkerFaceColor','r','LineWidth',1.7); title(handles.fig1,'HR tachogram');
    xlabel(handles.fig1,'Time [s]');
    ylabel(handles.fig1,'Heart rate[BPM]');
    grid(handles.fig1,'on');
    legend(handles.fig1,strcat('Signal:',num2str(c)));
    set(handles.sl1,'Value',time_cut{c,1});
    set(handles.sl2,'Value',time_cut{c,2});
    X1=signali{c,1};
    X2=signali{c,2};
    k1 = find(X1 <= time_cut{c,1});
    k1=k1(end);
    k2 = find(X1 >= time_cut{c,2});
    k2=k2(1);
    X1=X1(k1:k2);
    X2=X2(k1:k2);
    rr=[];
    for i = 1:length(X2)
        rr = [rr 60/X2];
    end
    signali_rr{c}=rr;
    xlim(handles.fig1,[X1(1) X1(end)]);
    ylim(handles.fig1,[data{c,3} data{c,2}]);


    ylim(handles.fig1,[data{c,3} data{c,2}]);
    set(handles.duzina, 'String', data{c,1}-data{c,8});
    set(handles.minimalni, 'String', data{c,2});
    set(handles.maksimalni, 'String', data{c,3});
    set(handles.edit18, 'String', data{c,4});
    set(handles.text34, 'String', data{c,5});
    set(handles.text35, 'String', data{c,6});
    set(handles.text9, 'String','[s]');
    set(handles.text7,'String','[bpm]');
    set(handles.text8,'String','[bpm]');
    set(handles.sl1,'Min',data{c,8});
    set(handles.sl2,'Min',data{c,8});
    set(handles.sl1,'Max',data{c,1});
    set(handles.sl2,'Max',data{c,1});
    set(handles.minsl,'String',num2str(time_cut{c,1}));
    set(handles.maxsl,'String',num2str(time_cut{c,2}));
else
    plot(handles.fig1,signali{c,1},signali_rr{c}, 'Color','b','Marker','o','MarkerFaceColor','r','LineWidth',1.7); title(handles.fig1,'RR tachogram');
    xlabel(handles.fig1,'Time [s]');
    ylabel(handles.fig1,'RR interval [s]');
    grid(handles.fig1,'on');
    legend(handles.fig1,strcat('Signal:',num2str(c)));
    set(handles.sl1,'Value',time_cut{c,1});
    set(handles.sl2,'Value',time_cut{c,2});
    X1=signali{c,1};
    X2=signali_rr{c,1};
    k1 = find(X1 <= time_cut{c,1});
    k1=k1(end);
    k2 = find(X1 >= time_cut{c,2});
    k2=k2(1);
    X1=X1(k1:k2);
    X2=X2(k1:k2);
     rr=[];
    for i = 1:length(X2)
        rr = [rr 60/X2(i)];
    end
    xlim(handles.fig1,[X1(1) X1(end)]);


    ylim(handles.fig1,[min(X2)-1 max(X2)+1]);
    set(handles.duzina, 'String', data{c,1}-data{c,8});
    set(handles.minimalni, 'String', data{c,2});
    set(handles.maksimalni, 'String', data{c,3});
    set(handles.edit18, 'String', data{c,4});
    set(handles.text34, 'String', data{c,5});
    set(handles.text35, 'String', data{c,6});
    set(handles.text9,'String','[s]');
    set(handles.text7,'String','[s]');
    set(handles.text8,'String','[s]');
    set(handles.sl1,'Min',data{c,8});
    set(handles.sl2,'Min',data{c,8});
    set(handles.sl1,'Max',data{c,1});
    set(handles.sl2,'Max',data{c,1});
    set(handles.minsl,'String',num2str(time_cut{c,1}));
    set(handles.maxsl,'String',num2str(time_cut{c,2}));
end

% hObject    handle to popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup


% --- Executes during object creation, after setting all properties.
function popup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton17.
function pushbutton17_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global time_cut;
global signali;
fi = uigetdir;
d=cell(0);
for c=1:size(signali,1)
X1=signali{c,1};
X2=signali{c,2};
k1 = find(X1 <= time_cut{c,1});k1=k1(end);
k2 = find(X1 >= time_cut{c,2});k2=k2(1);
X1=X1(k1:k2);
X2=X2(k1:k2);
X=[X1 X2];
d{c}=X;
end
ime=get(handles.studyname,'String');

save(fullfile(fi,strcat(ime,'.mat')),'d');


function studyname_Callback(hObject, eventdata, handles)
% hObject    handle to studyname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of studyname as text
%        str2double(get(hObject,'String')) returns contents of studyname as a double


% --- Executes during object creation, after setting all properties.
function studyname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to studyname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton18.
function pushbutton18_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pth = uigetdir;
info = dir(fullfile(pth, '*.mat'));
data=cell(0);

leg = cell(0);
N = max(size(info));
M = -1;
showhr = get(handles.showhr,'Value');
showrr = get(handles.showrr,'Value');
for qw=1:max(size(info))
    
    ia=load(fullfile(pth,info(qw).name));
    M = max(size(ia.d));
    leg= [leg strrep(info(qw).name,'.mat','')];
    data{qw}=ia.d;
end
if(showhr)
figure;
for i=1:N
    T=[];
    subplot(N,1,i);
    labels=cell(0);
    for j=1:M
        t=data{i};
        t=t(j);
        labels{j}=num2str(j);
        t=t{1};
        T=[T t(:,2)];
    end
       boxplot(T,'labels',labels);
       title(leg{i});
       ylabel('HR [BPM]');
       xlabel('Subject number');
       grid on;
end
minSt=zeros(N,M);
maxSt=zeros(N,M);
meanSt=zeros(N,M);
ranSt=zeros(N,M);
for i=1:N
    for j=1:M
        t=data{i};
        t=t(j);
        t=t{1};
        minSt(i,j) = min(t(:,2));
        maxSt(i,j) = max(t(:,2));
        meanSt(i,j) = mean(t(:,2));
        ranSt(i,j) = max(t(:,2)) - min(t(:,2));
    end
end
figure;
ind = 1:M;
color = [0 0 1; 1 0 0 ; 0 1 0; 1 0.5 1; 0 1 1; 1 1 0;  1 0.5 0; 0.5 1 0;0.5 0 1; 0.7 0 0.7; .3 .3 .3];

subplot(2,2,1);
for i=1:N
    plot(ind,minSt(i,:),'Color',color(i,:),'LineWidth',1.7,'Marker','o','MarkerFaceColor','k');
   
    set(gca,'XTick',ind);
    hold on;
end
xlabel('Subject number');
ylabel('Minimum HR [BPM]');
title('Minimum heart rate');
grid on;
legend(leg);
subplot(2,2,2);
for i=1:N
    plot(ind,maxSt(i,:),'Color',color(i,:),'LineWidth',1.7,'Marker','o','MarkerFaceColor','k');
    
    set(gca,'XTick',ind);
    hold on;
end
xlabel('Ispitanik broj');
ylabel('Maximum HR [BPM]');
title('Maximum heart rate');
legend(leg);
grid on;
subplot(2,2,3);
for i=1:N

    plot(ind,meanSt(i,:),'Color',color(i,:),'LineWidth',1.7,'Marker','o','MarkerFaceColor','k');
    
    set(gca,'XTick',ind);
    hold on;
end
xlabel('Subject number');
ylabel('Average HR [BPM]');
title('Averaged heart rate');
grid on;
legend(leg);
subplot(2,2,4);
for i=1:N
    plot(ind,ranSt(i,:),'Color',color(i,:),'LineWidth',1.7,'Marker','o','MarkerFaceColor','k');
   
    set(gca,'XTick',ind)
    hold on;
end
xlabel('Subject number');
ylabel('Range of HR [BPM]');
legend(leg);
title('Range of heart rate');
grid on;
end

if(showrr)
figure;
for i=1:N
    T=[];
    subplot(N,1,i);
    labels=cell(0);
    for j=1:M
        t=data{i};
        t=t(j);
        
        labels{j}=num2str(j);
        t=t{1};
        rr = [];
        for ww = 1:length(t(:,2))
            rr = [rr 60/t(ww,2)];
        end
        T=[T rr'];
    end
       boxplot(T,'labels',labels);
       title(leg{i});
       ylabel('RR [s]');
       xlabel('Ispitanik broj');
       grid on;
end
minSt=zeros(N,M);
maxSt=zeros(N,M);
meanSt=zeros(N,M);
ranSt=zeros(N,M);
for i=1:N
    for j=1:M
        t=data{i};
        t=t(j);
        t=t{1};
        rr = [];
        for ww = 1:length(t(:,2))
            rr = [rr 60/t(ww,2)];
        end
        minSt(i,j) = min(rr);
        maxSt(i,j) = max(rr);
        meanSt(i,j) = mean(rr);
        ranSt(i,j) = max(rr) - min(rr);
    end
end
figure;
subplot(2,2,1);
ind = 1:M;
color = [0 0 1; 1 0 0 ; 0 1 0; 1 0.5 1; 0 1 1; 1 1 0;  1 0.5 0; 0.5 1 0;0.5 0 1; 0.7 0 0.7; .3 .3 .3];

subplot(2,2,1);
for i=1:N
    plot(ind,minSt(i,:),'Color',color(i,:),'LineWidth',1.7,'Marker','o','MarkerFaceColor','k');
    hold on;
end
xlabel('Subject number');
ylabel('Minimum RR interval [s]');
title('Minimum RR interval');
grid on;
legend(leg);
subplot(2,2,2);
for i=1:N
    plot(ind,maxSt(i,:),'Color',color(i,:),'LineWidth',1.7,'Marker','o','MarkerFaceColor','k');
    hold on;
end
xlabel('Subject number');
ylabel('Maximum RR interval[s]');
title('Maximum RR interval');
legend(leg);
grid on;
subplot(2,2,3);
for i=1:N

    plot(ind,meanSt(i,:),'Color',color(i,:),'LineWidth',1.7,'Marker','o','MarkerFaceColor','k');
    hold on;
end
xlabel('Subject number');
ylabel('Averaged RR interval [s]');
title('Averaged RR interval');
grid on;
legend(leg);
subplot(2,2,4);
for i=1:N
    plot(ind,ranSt(i,:),'Color',color(i,:),'LineWidth',1.7,'Marker','o','MarkerFaceColor','k');
    hold on;
end
xlabel('Subject number');
ylabel('Range of RR intervals [s]');
legend(leg);
title('Range of RR intervals');
grid on;
end
header = cell(0);
for i=1:N
    header{i} = strcat('Studija :',num2str(i), ' minimalna');
end
for i=1:N
    header{i+N} = strcat('Studija :',num2str(i), ' maksimalna');
end
for i=1:N
    header{i+2*N} = strcat('Studija :',num2str(i), ' srednja');
end
for i=1:N
    header{i+3*N} = strcat('Studija :',num2str(i), ' opseg');
end
% csvwrite(fullfile(pth,'tabelarni_prikaz.csv'),header);
% dlmwrite(fullfile(pth,'tabelarni_prikaz.csv'),[minSt' maxSt' meanSt' ranSt'],'delimiter',',','-append');
% --- Executes on button press in openstudy.
function openstudy_Callback(hObject, eventdata, handles)
% hObject    handle to openstudy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname,pname] = uigetfile;


function locfile_Callback(hObject, eventdata, handles)
% hObject    handle to locfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of locfile as text
%        str2double(get(hObject,'String')) returns contents of locfile as a double


% --- Executes during object creation, after setting all properties.
function locfile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to locfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function saveloc_Callback(hObject, eventdata, handles)
% hObject    handle to saveloc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of saveloc as text
%        str2double(get(hObject,'String')) returns contents of saveloc as a double


% --- Executes during object creation, after setting all properties.
function saveloc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to saveloc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in hrsig.
function hrsig_Callback(hObject, eventdata, handles)
% hObject    handle to hrsig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of hrsig


% --- Executes on button press in sigrr.
function sigrr_Callback(hObject, eventdata, handles)
% hObject    handle to sigrr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of sigrr


% --- Executes on button press in showrr.
function showrr_Callback(hObject, eventdata, handles)
% hObject    handle to showrr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showrr


% --- Executes on button press in showhr.
function showhr_Callback(hObject, eventdata, handles)
% hObject    handle to showhr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showhr


% --- Executes on button press in deltarr.
function deltarr_Callback(hObject, eventdata, handles)
% hObject    handle to deltarr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global time_cut;
global signali;
global signali_rr
c=get(handles.popup,'Value');
X1=signali{c,1};
X2=signali_rr{c};
k1 = find(X1 <= time_cut{c,1});k1=k1(end);
k2 = find(X1 >= time_cut{c,2});k2=k2(1);
X1=X1(k1:k2);
X2=X2(k1:k2);


y = [0];
for i = 1:length(X2) - 1
    y = [y X2(i+1)-X2(i)];
end
[mn,~] = min(y);
[mx,~] = max(y);
figure;
plot(X1, y, 'Color','b', 'Marker','o','MarkerFaceColor','r','LineWidth',1.7);
title('Delta RR');
xlabel('Time [s]');
ylabel('Delta RR [s]');
hold on
grid on
axis([X1(1) X1(end) mn mx]);



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.

function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
