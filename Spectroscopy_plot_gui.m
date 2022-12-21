function varargout = Spectroscopy_plot_gui(varargin)
% SPECTROSCOPY_PLOT_GUI MATLAB code for Spectroscopy_plot_gui.fig
%      SPECTROSCOPY_PLOT_GUI, by itself, creates a new SPECTROSCOPY_PLOT_GUI or raises the existing
%      singleton*.
%
%      H = SPECTROSCOPY_PLOT_GUI returns the handle to a new SPECTROSCOPY_PLOT_GUI or the handle to
%      the existing singleton*.
%
%      SPECTROSCOPY_PLOT_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPECTROSCOPY_PLOT_GUI.M with the given input arguments.
%
%      SPECTROSCOPY_PLOT_GUI('Property','Value',...) creates a new SPECTROSCOPY_PLOT_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Spectroscopy_plot_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Spectroscopy_plot_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Spectroscopy_plot_gui

% Last Modified by GUIDE v2.5 14-Sep-2022 17:09:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Spectroscopy_plot_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @Spectroscopy_plot_gui_OutputFcn, ...
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


% --- Executes just before Spectroscopy_plot_gui is made visible.
function Spectroscopy_plot_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Spectroscopy_plot_gui (see VARARGIN)

% Choose default command line output for Spectroscopy_plot_gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Spectroscopy_plot_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Spectroscopy_plot_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% hObject    handle to loadfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in loadfile.
function loadfile_Callback(hObject, eventdata, handles)

global Wavelength
global Time
global Transient
global t1
global t2
global lambda_1
global lambda_2
global figure1
global C2
global C1
global fajl
global clrbar
global axesfont
global fontsize

fajl = uigetfile('*.dat');
timefajl = uigetfile('*.dat');
lambdafajl = uigetfile('*.dat');
Transient=dlmread(fajl);
Time=dlmread(timefajl);
Wavelength=dlmread(lambdafajl);
Time(end) = [];
Transient(:,end) = [];

current_folder = pwd;
set(handles.statictext1, 'String', current_folder);

current_folder = fajl;
set(handles.static2, 'String', fajl);

fontsize=10;
clrbar=10;
axesfont = 10;


set( handles.slider1, 'Max', Wavelength(end),'Min', Wavelength(1),'Value', Wavelength(1) );
set( handles.slider_spectrum, 'Max', Time(end),'Min', Time(1),'Value', Time(1) );

axes(handles.axes2);
figure1 = pcolor(handles.axes2,Time,Wavelength,Transient);
view(2)
colormap jet;
set(gca,'TickLength',[0.03 0.03]);
C1 =str2double(get(handles.C1,'String'));
C2 = str2double(get(handles.C2,'String'));
caxis([C1 C2]);
t1 = str2double(get(handles.Time_1,'String'));
t2 = str2double(get(handles.Time_2,'String'));
lambda_1 = str2double(get(handles.Wavelength1,'String'));
lambda_2 = str2double(get(handles.Wavelength2,'String'));
axis([t1 t2 lambda_1 lambda_2]);
set(gca,'XColor', 'k');
set(gca,'Layer','top');
set(gca,'FontSize',fontsize);
set(gca,'GridLineStyle','none');
colorbar('FontSize',clrbar);
shading interp;
xlabel('Delay [fs]','FontSize',axesfont,'Color','k');
ylabel('Wavelength (nm)','FontSize',axesfont,'Color','k');
addToolbarExplorationButtons(gcf)
hold on;
box on;

valna1 = Wavelength(round(end/2));
I = find(abs(Wavelength-valna1)<((Wavelength(end)-Wavelength(end-1))/2));a=round(size(I)/2);a = a(1);I = I(a);Wavelength_I = Wavelength(I)
axes(handles.axes3);
h1=plot(handles.axes3,Time,Transient(I,:));set([h1],'LineWidth',2)
grid on;
axis([t1 t2 -inf inf]);
addToolbarExplorationButtons(gcf)
legend({['\lambda_1 =' num2str(valna1) 'nm' ]},'FontSize',fontsize,'Location','best')
xlabel(['Delay[fs]'],'FontSize',axesfont);
ylabel('\DeltaA','FontSize',axesfont);
xt = get(gca, 'XTick');set(gca, 'FontSize', fontsize);box on;grid on;

K1 = round(size(Wavelength)/2);
K1 = K1(1);
y=20; %KOLIKO NM SE SMOOTHA
%SMOOTHANI SPEKTRI ODKOMENTIRATI
axes(handles.axes4);
g1=plot(handles.axes4,smooth(Transient(:,K1),y),Wavelength);
set([g1],'LineWidth',2)
set(gca,'linewidth',2)
line(xlim(), [0,0], 'LineWidth', 2, 'Color', 'k');
axis([-inf inf lambda_1 lambda_2]); %340 650
xlabel(['\DeltaA'],'FontSize',axesfont);
ylabel(['Wavelength (nm)'],'FontSize',axesfont);
xt = get(gca, 'XTick');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 40 40]);
set(gca, 'FontSize', fontsize)
box on;
grid on;

% hObject    handle to loadfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)

function Time_1_Callback(hObject, eventdata, handles)
plot2d(handles)

% --- Executes during object creation, after setting all properties.
function Time_1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Time_2_Callback(hObject, eventdata, handles)
plot2d(handles)


% --- Executes during object creation, after setting all properties.
function Time_2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Wavelength1_Callback(hObject, eventdata, handles)
global Time
plot_spectrum(Time(1),handles)
plot2d(handles)
W1 = str2double(get(handles.Wavelength1,'String'));
W2 = str2double(get(handles.Wavelength2,'String'));
set( handles.slider1, 'Max', W2,'Min', W1,'Value', W1 );


% --- Executes during object creation, after setting all properties.
function Wavelength1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Wavelength2_Callback(hObject, eventdata, handles)
global Time
plot_spectrum(Time(1),handles)
plot2d(handles)

W1 = str2double(get(handles.Wavelength1,'String'));
W2 = str2double(get(handles.Wavelength2,'String'));
set( handles.slider1, 'Max', W2,'Min', W1,'Value', W1 );

% --- Executes during object creation, after setting all properties.
function Wavelength2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function C1_Callback(hObject, eventdata, handles)
plot2d(handles)
% --- Executes during object creation, after setting all properties.
function C1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function C2_Callback(hObject, eventdata, handles)
plot2d(handles)
% --- Executes during object creation, after setting all properties.
function C2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%FUNKCIJA ZA PLOTANJE 2D
function plot2d(handles)
global Wavelength
global Time
global Transient
global t1
global t2
global lambda_1
global lambda_2
global figure1
global C2
global C1
global clrbar
global axesfont
global fontsize

axes(handles.axes2);
figure1 = pcolor(handles.axes2,Time,Wavelength,Transient);
view(2)
colormap jet;
set(gca,'TickLength',[0.03 0.03]);
C1 =str2double(get(handles.C1,'String'));
C2 = str2double(get(handles.C2,'String'));
caxis([C1 C2]);
t1 = str2double(get(handles.Time_1,'String'));
t2 = str2double(get(handles.Time_2,'String'));
lambda_1 = str2double(get(handles.Wavelength1,'String'));
lambda_2 = str2double(get(handles.Wavelength2,'String'));
axis([t1 t2 lambda_1 lambda_2]);
set(gca,'XColor', 'k');
set(gca,'Layer','top');
set(gca,'FontSize',fontsize);
set(gca,'GridLineStyle','none');
colorbar('FontSize',clrbar);
shading interp;
xlabel('Delay [fs]','FontSize',axesfont,'Color','k');
ylabel('Wavelength (nm)','FontSize',axesfont,'Color','k');
addToolbarExplorationButtons(gcf)
hold on;
box on;

%PLOTANJE SPEKTARA
function plot_spectrum(time1,handles)

global Wavelength
global Time
global Transient
global lambda_1
global lambda_2
global axesfont
global fontsize
global C1
global C2

B=transpose(Transient);

lambda_1 = str2double(get(handles.Wavelength1,'String'));
lambda_2 = str2double(get(handles.Wavelength2,'String'));


K1 = find(abs(Time-time1)<((abs(Time(end)-Time(end-1)))));
K1 = round(mean(K1));
Time(K1)
y=10; %KOLIKO NM SE SMOOTHA
%SMOOTHANI SPEKTRI ODKOMENTIRATI
axes(handles.axes4);
g1=plot(handles.axes4,smooth(Transient(:,K1),y),Wavelength);
set([g1],'LineWidth',2)
set(gca,'linewidth',2)
line(xlim(), [0,0], 'LineWidth', 2, 'Color', 'k');
axis([-inf inf lambda_1 lambda_2]); %340 650
xlabel(['\DeltaA'],'FontSize',axesfont);
ylabel(['Wavelength [nm]'],'FontSize',axesfont);
xt = get(gca, 'XTick');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 40 40]);
set(gca, 'FontSize', fontsize)
box on;
grid on;

%funkcija za plotanje cutova
function plot_cut(valna1,handles)

global Wavelength
global Time
global Transient
global t1
global t2
global fontsize
global axesfont
global h1



t1 =str2double(get(handles.Time_1,'String'));
t2 =str2double(get(handles.Time_2,'String'));
I = find(abs(Wavelength-valna1)<((Wavelength(end)-Wavelength(end-1))/2));a=round(size(I)/2);a = a(1);I = I(a);Wavelength_I = Wavelength(I);
axes(handles.axes3);
h1=plot(handles.axes3,Time,Transient(I,:),'DisplayName',[num2str(round(valna1)) ' nm']);set([h1],'LineWidth',2)
grid on;
axis([t1 t2 -inf inf]);
addToolbarExplorationButtons(gcf)
lgd = legend();

% legend({[ h1 ]},'FontSize',fontsize,'Location','best')
xlabel(['Delay [fs]'],'FontSize',axesfont);
ylabel('\DeltaA','FontSize',axesfont);
xt = get(gca, 'XTick');set(gca, 'FontSize', fontsize);box on;grid on;


% --- Executes on button press in Background.
function Background_Callback(hObject, eventdata, handles)
global Time
global Transient

average = str2double(get(handles.bckg,'String'));
a = size(Time);
b = size(Transient);
vektor = zeros( [b(1) 1 ]);
r = 0;
for i=1:b(1)
    for j=1:average
        r = r + Transient(i,j);
    end
    vektor(i) = r/average;
    r = 0;
end
for j=1:b(2)
    Transient(:,j) = Transient(:,j)-vektor;
end
% dlmwrite('C:\Users\Femtolab\Desktop\rhodamin basaric gusce.dat',Transient)
figure(2)
plot(vektor,'LineWidth',2);
xlabel('Wavelength','FontSize',30,'Color','k');
ylabel('\DeltaA','FontSize',30,'Color','k');
set(gca,'FontSize',30);
addToolbarExplorationButtons(gcf)
grid on;
box on;

PLot_2d_Callback(hObject, eventdata, handles)

% hObject    handle to Background (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function bckg_Callback(hObject, eventdata, handles)
% hObject    handle to bckg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bckg as text
%        str2double(get(hObject,'String')) returns contents of bckg as a double


% --- Executes during object creation, after setting all properties.
function bckg_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bckg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider_spectrum_Callback(hObject, eventdata, handles)

time1 = get(hObject,'Value');

plot_spectrum(time1,handles);
set(handles.Spec_time2,'String', [num2str(round(time1/1000)), 'ps']);
% hObject    handle to slider_spectrum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_spectrum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_spectrum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
valna1 = get(hObject,'Value');
plot_cut(valna1, handles);
% set(hObject,'Min') = str2double(get(handles.Wavelength1,'String'));
% set(hObject,'Max') = str2double(get(handles.Wavelength2,'String'));


% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in Click_cut.
function Click_cut_Callback(hObject, eventdata, handles)

[click_x, click_y] = getpts(handles.axes2); %biranje tocaka mišem

%ako je multiple curves oznaceno plotaj ih sve, inace samo prvi klik

    
for i = 1:size(click_y)
    plot_cut(click_y(i), handles);
    hold on;
end

hold off;
% hObject    handle to Click_cut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in plot_cut_button.
function plot_cut_button_Callback(hObject, eventdata, handles)
global h1

valna1 = str2double(get(handles.valna_cut,'String'));
a = get(handles.Multiple_cuts,'Value');


if a == 1
    plot_cut(valna1,handles);
    hold on;
else
    hold off;
    plot_cut(valna1,handles);
    hold off;
end
% hObject    handle to plot_cut_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Multiple_cuts.
function Multiple_cuts_Callback(hObject, eventdata, handles)
% hObject    handle to Multiple_cuts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Multiple_cuts



function valna_cut_Callback(hObject, eventdata, handles)
% hObject    handle to valna_cut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of valna_cut as text
%        str2double(get(hObject,'String')) returns contents of valna_cut as a double


% --- Executes during object creation, after setting all properties.
function valna_cut_CreateFcn(hObject, eventdata, handles)
% hObject    handle to valna_cut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Spectrum_button.
function Spectrum_button_Callback(hObject, eventdata, handles)

global Transient

B=transpose(Transient);

% time1 = 10000; 
time1 = str2double(get(handles.Vrijeme_za_spektar,'String'));
a = get(handles.Multiple_spectra,'Value');

if a == 1
    plot_spectrum(time1,handles);
    hold on;
else
    hold off;
    plot_spectrum(time1,handles);
    hold off;
end

set(handles.Spec_time2,'String', [num2str(round(time1/1000)), 'ps']);

% hObject    handle to Spectrum_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Spectral_click.
function Spectral_click_Callback(hObject, eventdata, handles)

[click_x, click_y] = getpts(handles.axes2); %biranje tocaka mišem

% time1 = click_x(end);

for i = 1:size(click_y)
    plot_spectrum(click_x(i), handles);
    hold on;
end
hold off;

set(handles.Spec_time2,'String', [num2str(round(click_x(end)/1000)), 'ps']);

% hObject    handle to Spectral_click (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Multiple_spectra.
function Multiple_spectra_Callback(hObject, eventdata, handles)
% hObject    handle to Multiple_spectra (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Multiple_spectra



function Vrijeme_za_spektar_Callback(hObject, eventdata, handles)
% hObject    handle to Vrijeme_za_spektar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Vrijeme_za_spektar as text
%        str2double(get(hObject,'String')) returns contents of Vrijeme_za_spektar as a double


% --- Executes during object creation, after setting all properties.
function Vrijeme_za_spektar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Vrijeme_za_spektar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_chirp_click.
function pushbutton_chirp_click_Callback(hObject, eventdata, handles)

global Wavelength
global Time
global Transient
global t1
global t2
global lambda_1
global lambda_2
global Disperzija


% plot2d(handles);

[Chirp_x, Chirp_y] = getpts(handles.axes2); %biranje tocaka mišem
scatter(Chirp_x,Chirp_y,250,'x','MarkerEdgeColor',[1 1 1],'MarkerFaceColor',[0 0 0],'LineWidth',2);
hold on
% plot(Chirp_x,Chirp_y,'LineWidth',2); %plot koji spaja križiæe koji se dobiju nakon klikanja
hold on;

yfit1 = fittype('a0+a1*x+a2*x^2+a3*x^3','ind','x');
opts=fitoptions(yfit1);
set(opts,'robust','on','Algorithm', 'Levenberg-Marquardt','tolfun',1e-12,'tolx',1e-6,'maxiter',800,'MaxFunEvals',2000, 'StartPoint',[Chirp_x(1) 1 0.1 0.01]);
[fr1, gof1]=fit(Chirp_x, Chirp_y, yfit1,opts);
hold on;


x = Time;
a0 = fr1.a0; a1 = fr1.a1; a2 = fr1.a2; a3 = fr1.a3;
y = a0+a1*x+a2*x.^2+a3*x.^3;
plot(x, y,'--', 'LineWidth',2,  'Color', [1 1 1]);
axis([t1 t2 lambda_1 lambda_2]); %MORA BITI ISTO KAO I POCETNI AXIS, NA OVIM OSIMA CRTA POLINOM
addToolbarExplorationButtons(gcf)
hold on;
box on;

% TREBA IMPLEMENTIRATI

% X0 = find(abs(x-(Chirp_x(1)-200))<((x(2)-x(1)))); a=round(size(X0)/2);a = a(1);X0 = X0(a);x(X0) %kako ne bi uzeo slucajno tocku polinoma koja je na npr 15000fs jer se on i do tamo proteze (kao sinus npr), gleda polinom od prve kliknute tocke - 500 pa do zadnja kliknuta tocka(gornje valne duljine) + 500
% Xn = find(abs(x-(Chirp_x(end)+500))<((x(2)-x(1)))); bb=round(size(Xn)/2);bb = bb(1);Xn = Xn(bb);x(Xn)
% x = x(X0:Xn);
% y = y(X0:Xn);

Polinom = [x y];
Polinom_interp = interp1(Polinom(:,2), Polinom(:,1), Wavelength);
Disperzija = Polinom_interp;
Disperzija = floor(Disperzija); %zaokruzivanje disperzije, disperzija su zapravo sada koraci od 1fs, tj stupci matrice, ako u disperziji stoji 13 to znaci da tu valnu duljinu treba pomaknuti ulijevo za 13 stupaca
u=Disperzija-Disperzija(1); %varijabla za provjeru disperzije, treba biti otprilike 10-20 pa na vecim pikselima do 2000 npr


b=size(Transient);
matrica = zeros([b(1) b(2)]);

for i=1:b(1)
    
    Time2 = Time + (Disperzija(i)-Disperzija(1));
    a = interp1(Time, Transient(i,:), Time2);
    matrica(i,:) = a;
%     Ravnanje = i;
    
end


% CRTANJE IZRAVNATE
figure(3);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 36 36]);
axes('color','black');
addToolbarExplorationButtons(gcf)
set(gca,'XColor', 'k','YColor', 'k');
surf(Time,Wavelength,matrica);
axis([t1 t2 lambda_1 lambda_2]);
% axis([-inf inf -inf inf]);
set(gca,'FontSize',40);
colormap jet;
set(gca,'TickLength',[0.03 0.03]);
caxis([-0.0008 0.0038]);
% caxis([-inf inf]);
set(gca,'XColor', 'k');
set(gca,'Layer','top');
set(gca,'GridLineStyle','none');
colorbar('FontSize',40);
shading interp;
% xticks([0 3000 13000 130000 1300000])
% xticklabels({'})
xlabel('Delay [fs]','FontSize',40,'Color','k');
ylabel('Wavelength (nm)','FontSize',40,'Color','k');
box on

% hObject    handle to pushbutton_chirp_click (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Save_button.
function Save_button_Callback(hObject, eventdata, handles)

global fajl
global matrica
global Time
global Wavelength

name_matrix = strcat('CORRECTED',fajl);
name_times = strcat('CORRECTED_Times_',fajl);
name_wv = strcat('CORRECTED_Wavelength_',fajl);

dlmwrite(name_matrix, matrica,'delimiter',' ');
dlmwrite(name_times, Time,'delimiter',' ');
dlmwrite(name_wv, Wavelength,'delimiter',' ');

% hObject    handle to Save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function W1_Callback(hObject, eventdata, handles)
% hObject    handle to W1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of W1 as text
%        str2double(get(hObject,'String')) returns contents of W1 as a double


% --- Executes during object creation, after setting all properties.
function W1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to W1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function W2_Callback(hObject, eventdata, handles)
% hObject    handle to W2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of W2 as text
%        str2double(get(hObject,'String')) returns contents of W2 as a double


% --- Executes during object creation, after setting all properties.
function W2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to W2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cut_matrix.
function cut_matrix_Callback(hObject, eventdata, handles)

global Transient
global Wavelength
global W1
global W2
global I_cut
global J_cut
% global Time
% global t1
% global t2
% global figure1
% global C2
% global C1
% global fontsize
% global clrbar
% global axesfont

W1 = str2double(get(handles.W1,'String'));
W2 = str2double(get(handles.W2,'String'));
I_cut = find(abs(Wavelength-W1)<((Wavelength(end)-Wavelength(end-1))/2));a=round(size(I_cut)/2);a = a(1);I_cut = I_cut(a);Wavelength(I_cut);J_cut = find(abs(Wavelength-W2)<((Wavelength(end)-Wavelength(end-1))/2));a=round(size(J_cut)/2);a = a(1);J_cut = J_cut(a);Wavelength(J_cut)
assignin('base','Transient',Transient(I_cut:J_cut,:));
assignin('base','Wavelength',Wavelength(I_cut:J_cut));

Transient = Transient(I_cut:J_cut,:);
Wavelength = Wavelength(I_cut:J_cut);

% plot2d(handles)
set(handles.Wavelength1, 'String', W1);
set(handles.Wavelength2, 'String', W2);
Wavelength2_Callback(hObject, eventdata, handles);




% hObject    handle to cut_matrix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object creation, after setting all properties.
function Spectrum_text_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Spectrum_text (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function Spec_time2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Spec_time2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
