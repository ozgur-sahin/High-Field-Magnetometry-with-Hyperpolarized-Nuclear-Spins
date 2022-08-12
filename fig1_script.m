%% 1) Resonant AC field and pulses
clear;clc;close all;
L=8; %number of cycles or pulses
%pvec=[-0.5+(1:L) L]; %equally spaced pulses
pvec=(1:L); %equally spaced pulses

if contains(pwd,'AjoyLabBlueOak')
    sys=1; %1 for lab PC
    %dirname='C:\Users\AjoyLabBlueOak\OneDrive\Documents\MATLAB\Berkeley_Lab\phase paper figs\fig3\';
else %for laptop
    sys=0;
    %dirname='D:\OneDrive\Documents\MATLAB\lab data\';
end


% WO sequence
T=73; %period
phase=0;  %set between 0 and pi

hcon=start_fig(7,[4 1]);
x=linspace(0,(L)*T,500);
y=30*sin(2*pi/(4*T)*x)+100;
plot_preliminaries(x,y,1,'nomarker');
set(gca,'ylim',[-10 160]);
hold on
yline(0);
for j=0:L-1
    p1=plot_yline(j*T+39,5);set(p1,'linestyle','--');
    p1=plot_yline(j*T+71,5);set(p1,'linestyle','--');
    rectangle('Position',[j*T 0 30 60]);
    hold on
end

%set(gca,'xlim',[0.5 L+1]);
plot_labels('Time','Field [au]');

if sys==0 %for laptop
    dirname='D:\OneDrive\Documents\MATLAB\Berkeley_Lab\mag paper figs\fig1\';
else %for blue oak
    dirname='C:\Users\AjoyLabBlueOak\OneDrive\Documents\MATLAB\Berkeley_Lab\mag paper figs\fig1\';   
end
savefig(hcon,[dirname 'concept_figure.fig']);

%% 8 raw data windows

explimits{1}='2021-06-08-000001_Proteus_chunk2_1';
explimits{2}='2021-06-08-000001_Proteus_chunk2_8';

if contains(pwd,'AjoyLabBlueOak')
    sys=1; %1 for lab PC
else
    sys=0; %laptop
end

selected_files1= transpose(agilent_read_proteus(explimits,sys));
if sys==0 %for laptop
    dirname2='D:\OneDrive\Documents\MATLAB\lab data\';
else %for blue oak
    dirname2='C:\Users\AjoyLabBlueOak\OneDrive\Documents|MATLAB\lab data\';
end

for j=1:length(selected_files1)
    A=load([dirname2 selected_files1{j}]);
    y3{j}=A.pulsechunk;
end
x3=(1:length(y3{1}))*1e-9;

%%
for ind=1:length(selected_files1)
    h3(ind)=start_fig(ind+2,[2 2]);
    p3(j)=plot_preliminaries(x3,y3{ind},1,'nomarker');
    plot_labels('Time [s]','Signal [au]');
    ylim([2050 2350]);
    savefig(h3(ind),[dirname 'raw8osc_window_' num2str(ind) '.fig']);
    
    Y3{ind}=abs(fftshift(fft(y3{ind}-mean(y3{ind}))));
    F3{ind}=fftfreq(x3,true);
    Ymax(ind)=max(Y3{ind});
end

h4=start_fig(14,[2 2]);
p4=plot_preliminaries(1:length(Ymax),Ymax,1);
% SP=splinefit(1:length(Ymax),Ymax,5,2);
% pp=ppval(SP,1:length(Ymax));
% pspl=plot_preliminaries(1:length(Ymax),pp,2,'nomarker');
plot_labels('Window count','Signal [au]');
savefig(h4,[dirname 'windowFT.fig']);

%% raw data and its FT
%clear;clc;close all;

if sys==1 %for blue oak
    A=load('C:\Users\AjoyLabBlueOak\Documents\MATLAB\lab data\2021-05-30-172033_Proteus_chunk8');
else %for laptop
    A=load('D:\OneDrive\Documents\MATLAB\lab data\2021-05-30-172033_Proteus_chunk8');
end
y=A.pulsechunk(:);
clear A;

yt=y(1:1000); %truncated
x=(1:length(yt));

h=start_fig(1,[2 2]);
p=plot_preliminaries(x,yt-mean(yt),2,'nomarker');
plot_labels('Time [ns]','Signal [au]');

if sys==1 %blue oak
    dirname='C:\Users\AjoyLabBlueOak\OneDrive\Documents\MATLAB\Berkeley_Lab\mag paper figs\fig1\';
else %laptop
    dirname='D:\OneDrive\Documents\MATLAB\Berkeley_Lab\mag paper figs\fig1\';
end
savefig(h,[dirname 'raw_oscillation.fig']);


ytt=y(1:100000);
%ytt=y;
h2=start_fig(2,[2 2]);
Ts=1e-9;
Y=abs(fft(ytt-mean(ytt)))*Ts*1e3;

bw=1/Ts;
n=length(ytt);
F=linspace(0,bw,n);
pf=plot_preliminaries(F/1e6,Y,2,'nomarker'); %plotting in the MHz scale
xlim([0 40]);
ylim([-max(Y)*3e-2 1.1*max(Y)]);
plot_labels('Frequency [MHz]','Signal [au]');
savefig(h2,[dirname 'FT_raw_oscillation.fig']);

