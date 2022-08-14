clear;clc;close all;

optim=optimoptions('lsqcurvefit','display','off');
explimits1{1}='2021-05-19-235457_Proteus';
explimits1{2}='2021-05-20-054229_Proteus';
inputArray1=[10:10:60]; %pulse width us

explimits2{1}='2021-05-21-234002_Proteus';
explimits2{2}='2021-05-22-042842_Proteus';
inputArray2=[15:10:55]; %pulse width us

explimits3{1}='2021-05-29-022651_Proteus';
explimits3{2}='2021-05-29-121110_Proteus';
inputArray3=[12.5:5:57.5]; %pulse width us

inputArray=[inputArray1 inputArray2 inputArray3];
[inputArray, IX]=sort(inputArray);

ex=[min(inputArray) max(inputArray)]; %min max points for input
rep=10;

IX2=[];
for ind=1:length(IX)
    IX2=[IX2 (IX(ind)-1)*rep+1:IX(ind)*rep]; 
end

if contains(pwd,'AjoyLabBlueOak')
    sys=1; %1 for lab PC
else
    sys=0; %laptop
end

selected_files1= transpose(agilent_read_proteus(explimits1,sys));
selected_files2= transpose(agilent_read_proteus(explimits2,sys));
selected_files3= transpose(agilent_read_proteus(explimits3,sys));

selected_files=[selected_files1' selected_files2' selected_files3'];
selected_files=selected_files(IX2);
clear selected_files1 selected_files2 selected_files3

for j=1:length(selected_files)
    [S(j),xfull{j},yfull{j},x{j},y{j},baseline{j}] = process_data_proteus(selected_files{j},sys);
end
%%
IX3=find(inputArray==30);
IX4=[(IX3-1)*rep+2:IX3*rep];
%IX4=(IX3-1)*rep+21;
Yavg=0;
xavg=x{IX4(1)};
Ts=xavg(2)-xavg(1);
for ind=IX4
    Y{ind}=abs(fftshift(fft(y{ind}-baseline{ind})))*Ts;
    F=fftfreq(x{ind},true);
    %selected_files{ind}
    Yavg=Yavg+Y{ind};
end
Yavg=Yavg/length(IX4);


%% fig5a
if sys==1
    dirname='C:\Users\AjoyLabBlueOak\OneDrive\Documents\MATLAB\Berkeley_Lab\paper figs\fig5\';
else
    dirname='D:\OneDrive\Documents\MATLAB\Berkeley_Lab\mag paper figs\fig5\';
end

h1=start_fig(1,[6 2]); 
p1=plot_preliminaries(F,smooth(Yavg,100),1,'nomarker');
xline(1000);
xline(4000);
xlim([0 max(F)]);
plot_labels('Frequency [Hz]','Signal [au]');
savefig(h1,[dirname 'fig5a_main.fig']);

%% inset on how the experiment is conducted
L=8; %number of cycles or pulses
%pvec=[-0.5+(1:L) L]; %equally spaced pulses
pvec=(1:L); %equally spaced pulses


% WO sequence
T=73; %period
phase=0;  %set between 0 and pi

hpuls=start_fig(2,[4 1]);
hold on
yline(0);
for j=0:L-1
    p1=plot_yline(j*T+39,5);set(p1,'linestyle','--');
    p1=plot_yline(j*T+71,5);set(p1,'linestyle','--');
    rectangle('Position',[j*T 0 30 60]);
    hold on
end

%set(gca,'xlim',[0.5 L+1]);
ylim([-5 70]);
plot_labels('Time','Field [au]');

savefig(hpuls,[dirname 'fig5a_insetpulse.fig']);

hchirp=start_fig(3,[3 1]);
xchirp=linspace(0,20,500);
ychirp=100*sin((10+1.5*xchirp).*xchirp/12);
xchirp=[xchirp xchirp+20 40+xchirp];
ychirp=[ychirp ychirp ychirp];
pchirp=plot_preliminaries(xchirp,ychirp,2,'nomarker');

savefig(hchirp,[dirname 'fig5a_insetchirp.fig']);

%% fig5b

hfirst=start_fig(4,[3 1]); %first harmonic zoom
relind=find((F>2300)&(F<3200));
Ftrunc=F(relind);
ytrunc=smooth(Yavg(relind),200);
pfirst=plot_preliminaries(Ftrunc,ytrunc,1,'noline');
set(pfirst,'MarkerSize',2);

SF=splinefit(Ftrunc,ytrunc,20,3);
yy=ppval(SF,Ftrunc);
pspl=plot_preliminaries(Ftrunc,yy,2,'nomarker');

xlim([2300 3200]);
plot_labels('Frequency [Hz]','Signal [au]');
savefig(hfirst,[dirname 'fig5b_first_harmonic.fig']);

%%
hsec=start_fig(5,[3 1]);  %second harmonic zoom
relind2=find((F>4800)&(F<6200));
psec=plot_preliminaries(F(relind2),smooth(Yavg(relind2),200),1,'noline');
set(psec,'MarkerSize',3);

xlim([4800 6200]);

gaussfun= @(par,x) par(1)*exp(-(x-par(4)).^2/par(2))+par(3);
lorefun= @(par,x) par(1)./(1+2/par(3)*(x-par(2)).^2)+par(4); %par1 Amplitude; par2 center; par3 FWHM
par0=[0.6 5520 200 0]; % for lorentzian
%par0=[0.6 200 0 5520]; %for gaussian

Ysmoothed=smooth(Yavg(relind2),1000)';
%fitgau=lsqcurvefit(gaussfun,par0,F(relind),Ysmoothed);
fitlor=lsqcurvefit(lorefun,par0,F(relind2),Ysmoothed);
xfit=linspace(4800,6200,5000);
%yfit=gaussfun(fitgau,xfit);
yfit=lorefun(fitlor,xfit);
plot_preliminaries(xfit,yfit,2,'nomarker');

maxp=max(yfit);
HM=maxp/2;
FWHM1=fitlor(2)+sqrt(fitlor(3)/2*(fitlor(1)/(HM-fitlor(4))-1));
FWHM2=fitlor(2)-sqrt(fitlor(3)/2*(fitlor(1)/(HM-fitlor(4))-1));
FWHM=[FWHM1 FWHM2];
plot_preliminaries(FWHM,lorefun(fitlor,FWHM),3);

plot_labels('Frequency [Hz]','Signal [au]');
savefig(hsec,[dirname 'fig5b_second_harmonic.fig']);

%% fig5c time domain signal

h6=start_fig(6,[6 2]);
p6=plot_preliminaries(x{85},y{85},1,'nomarker');
xlim([0 20]);
plot_labels('Time [s]','Signal [au]');
savefig(h6,[dirname 'fig5c.fig']);

%% fig5d oscillations

clear relind
relind{1}=find(x{85}>6.45 & x{85}<6.6);
relind{2}=find(x{85}>10.45 & x{85}<10.6);
relind{3}=find(x{85}>16.35 & x{85}<16.5);

[h7, ax7]=start_fig(7,[6 2]); %dimensions were [6 2] for the paper and [2 2] for the bottom figures
osc=y{85}-baseline{85};
p7=plot_preliminaries(x{85},osc,1,'nomarker');
set(p7,'Color','#FF4500');
yline(0);
rectangle('Position',[x{85}(relind{1}(1)) -5 x{85}(relind{1}(end))-x{85}(relind{1}(1)) 10],...
    'LineStyle','--','EdgeColor','b');
rectangle('Position',[x{85}(relind{2}(1)) -15 x{85}(relind{2}(end))-x{85}(relind{2}(1)) 25],...
    'LineStyle','--','EdgeColor','b');
rectangle('Position',[x{85}(relind{3}(1)) -5 x{85}(relind{3}(end))-x{85}(relind{3}(1)) 10],...
    'LineStyle','--','EdgeColor','b');
plot_labels('Time [s]','Signal [au]');
mid=(max(osc)+min(osc))/2;
rang=30;
ylim([mid-rang mid+rang]);
xlim([0 20]);
savefig(h7,[dirname 'fig5d.fig']);
saveas(h7,[dirname 'fig5d.pdf'],'pdf');

%% partial FTs

gaussfun= @(par,x) par(1)*exp(-(x-par(2)).^2/par(3))+par(4); %par1 amplitude; par2 center; par3 FWHM
lorefun= @(par,x) par(1)./(1+2/par(3)*(x-par(2)).^2)+par(4); %par1 amplitude; par2 center; par3 FWHM
fitif=false;

h8=start_fig(8,[2 2]);
relind=find(x{85}>6.45 & x{85}<6.6);
%p8=plot_preliminaries(x{85}(relind),osc(relind),1,'nomarker');
Y8=abs(fftshift(fft(osc(relind))))*Ts;
F8=fftfreq(x{85}(relind),true);
relind2=find(F8<4000 & F8>0);
Y8=Y8(relind2);F8=F8(relind2);
if ~fitif
    p8f=plot_preliminaries(F8,Y8,1);
else
    p8f=plot_preliminaries(F8,Y8,1);
    par0=[2.6e-3 2082 400 5e-4];
    lb=[1e-2 2000 10 3.9e-4];
    ub=[2e-1,2200,inf,1.2e-3];
    %lorfit=lsqcurvefit(lorefun,par0,F8,Y8,lb,ub,optim);
    gaussfit=lsqcurvefit(gaussfun,par0,F8,Y8,lb,ub,optim);
    xfit=linspace(0,4000,1000);
    %yfit=lorefun(lorfit,xfit);
    yfit=gaussfun(gaussfit,xfit);
    p8fit=plot_preliminaries(xfit,yfit,2,'nomarker');
end
set(p8f,'Color','b');
set(p8f,'MarkerFaceColor','b');
xlim([1900 2300]);
%xlim([2000 3600]);
plot_labels('Frequency [Hz]','Signal [au]');
savefig(h8,[dirname 'fig5e_1.fig']);

%%
h9=start_fig(9,[2 2]);
relind=find(x{85}>10.45 & x{85}<10.6);
%p9=plot_preliminaries(x{85}(relind),osc(relind),1,'nomarker');
Y9=abs(fftshift(fft(osc(relind))))*Ts;
F9=fftfreq(x{85}(relind),true);
relind2=find(F9<4000 & F9>0);
Y9=Y9(relind2);F9=F9(relind2);
if ~fitif
    p9f=plot_preliminaries(F9,Y9,1);
end
if fitif
    p9f=plot_preliminaries(F9,Y9,1);
    par0=[1e-2 2654 200 5e-4];
    lb=[0 2500 0 3.9e-4];
    ub=[1e-1,2700,1000^2,0.5e-3];
    lorfit=lsqcurvefit(lorefun,par0,F9,Y9,lb,ub,optim);
    gaussfit=lsqcurvefit(gaussfun,par0,F9,Y9,lb,ub,optim);
    xfit=linspace(0,4000,1000);
    %yfit=lorefun(lorfit,xfit);
    yfit=gaussfun(gaussfit,xfit);
    p9fit=plot_preliminaries(xfit,yfit,2,'nomarker');
end
set(p9f,'Color','b');
set(p9f,'MarkerFaceColor','b');
xlim([2400 2900]);
%xlim([2000 3600]);
ylim([-0.005 0.1]);
plot_labels('Frequency [Hz]','Signal [au]');
savefig(h9,[dirname 'fig5e_2.fig']);

%%
h10=start_fig(10,[2 2]);
relind=find(x{85}>16.35 & x{85}<16.5);
%p10=plot_preliminaries(x{85}(relind),osc(relind),1,'nomarker');
Y10=abs(fftshift(fft(osc(relind))))*Ts;
F10=fftfreq(x{85}(relind),true);
if ~fitif
    p10f=plot_preliminaries(F10,Y10,1);
end
if fitif
    p10f=plot_preliminaries(F10,Y10,1);
    relind2=find(F10<4000 & F10>3300);
    Y10=Y10(relind2);F10=F10(relind2);
    par0=[1e-2 3500 200 1.5e-3];
    lb=[3e-3 3400 10 6e-4];
    ub=[3e-2,4000,1000^2,5e-3];
    lorfit=lsqcurvefit(lorefun,par0,F10,Y10,lb,ub,optim);
    gaussfit=lsqcurvefit(gaussfun,par0,F10,Y10,lb,ub,optim);
    xfit=linspace(0,4000,1000);
    %yfit=lorefun(lorfit,xfit);
    yfit=gaussfun(gaussfit,xfit);
    p10fit=plot_preliminaries(xfit,yfit,2,'nomarker');
end
set(p10f,'Color','b');
set(p10f,'MarkerFaceColor','b');
xlim([3200 3800]);
%xlim([2000 3600]);
plot_labels('Frequency [Hz]','Signal [au]');
savefig(h10,[dirname 'fig5e_3.fig']);

%% short-time FT for SI

%shortFT=abs(stft(osc,3));
% d = seconds(1e-3);
% figure
% stft(osc,d,'OverlapLength',98);
% figure
% plot(shortFT(1,:))
windowno=500; %number of windows of FT
pointsno=floor(length(osc)/windowno); %no of points in each FT
oscmat=reshape(osc(1:pointsno*windowno),pointsno,windowno);
STFTmat=zeros(pointsno,windowno);
freqax=fftfreq(x{85}(1:pointsno),false);
%  figure(11);
% hold on

for ind=1:windowno
    STFTmat(:,ind)=abs(fft(oscmat(:,ind)));
    %plot(freqax,STFTmat(:,ind));
end
[h12,ax12]=start_fig(12,[4 2]);
image([0,max(x{85})],[min(freqax/1000),max(freqax/1000)],STFTmat);
xlabel(ax12,'Time [s]');
ylabel(ax12,'Frequency [kHz]');
set(ax12,'Ydir','normal')
ylim(ax12,[0 max(freqax/1000)/2]);
map=flip([uint8(floor(linspace(0,255,256)))' uint8(floor(linspace(0,255,256)))' 255*ones(256,1)]);
colormap('hot');
cb=colorbar();
title(cb,'Signal Intensity [au]');
hold on;

saveas(h12,[dirname 'freq_chirp_response.svg'],'svg'); 
