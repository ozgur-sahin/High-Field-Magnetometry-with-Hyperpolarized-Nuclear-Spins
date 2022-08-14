% %% fig6 graphs
%
% clear;clc;close all;
%
% explimits1{1}='2021-05-19-235457_Proteus';
% explimits1{2}='2021-05-20-054229_Proteus';
% inputArray1=[10:10:60]; %pulse width us
%
% explimits2{1}='2021-05-21-234002_Proteus';
% explimits2{2}='2021-05-22-042842_Proteus';
% inputArray2=[15:10:55]; %pulse width us
%
% explimits3{1}='2021-05-29-022651_Proteus';
% explimits3{2}='2021-05-29-121110_Proteus';
% inputArray3=[12.5:5:57.5]; %pulse width us
%
% %inputArray=[inputArray1 inputArray2 inputArray3];
% inputArray=[inputArray1 inputArray2];
% [inputArray, IX]=sort(inputArray);
%
% ex=[min(inputArray) max(inputArray)]; %min max points for input
% rep=10;
%
% IX2=[];
% for ind=1:length(IX)
%     IX2=[IX2 (IX(ind)-1)*rep+1:IX(ind)*rep];
% end
%
% if contains(pwd,'AjoyLabBlueOak')
%     sys=1; %1 for lab PC
% else
%     sys=0; %laptop
% end
%
% selected_files1= transpose(agilent_read_proteus(explimits1,sys));
% selected_files2= transpose(agilent_read_proteus(explimits2,sys));
% selected_files3= transpose(agilent_read_proteus(explimits3,sys));
%
% %selected_files=[selected_files1' selected_files2' selected_files3'];
% selected_files=[selected_files1' selected_files2'];
% selected_files=selected_files(IX2);
% clear selected_files1 selected_files2 selected_files3
%
% for j=1:length(selected_files)
%     [S(j),xfull{j},yfull{j},x{j},y{j},baseline{j}] = process_data_proteus(selected_files{j},sys);
% end
%
% clear S xfull yfull;
%
% dirname='C:\Users\AjoyLabBlueOak\OneDrive\Documents\MATLAB\Berkeley_Lab\paper figs\fig6\';
% dirname2='C:\Users\AjoyLabBlueOak\OneDrive\Documents\MATLAB\Berkeley_Lab\paper figs\fig6\fig6a\';
%
% %% fig6a
%
% ex=[min(inputArray) max(inputArray)]; %min max points for input
% rep=10;
%
% IX2=[];
% for ind=1:length(IX)
%     IX2=[IX2 (IX(ind)-1)*rep+1:IX(ind)*rep];
% end
%
%
% h1=start_fig(1,[4 2]);
%
% for ind=1:length(inputArray)
%     for j=1:rep
%         IX4=(ind-1)*rep+j;
%         %Smooth data
%         baseline = smooth(y{IX4},1000)';
%
%         % mean(y)
%         y_sub=y{IX4}-baseline;
%
%         x_sub=x{IX4};
%         Ts=x_sub(2)-x_sub(1);
%         Y=fft(y_sub);
%         Yabs(j,:)=abs(fftshift(Y))*Ts;
%         Ybase=abs(fftshift(fft(baseline)))*Ts;
%
%         n=size(Y,2);
%         fs=1/max(x_sub);
%         %         F=(-n/2+1:n/2)*(fs);
%         F{ind}=linspace(-1/2/Ts,1/2/Ts-fs,n);
%     end
%     y_avg{ind}=smooth((sum(Yabs(1:rep,:),1)/rep),1000)';
%     [lol peakind(ind)]=max(y_avg{ind});
%     peaks(ind)=F{ind}(peakind(ind));
%     p=plot_preliminaries(F{ind},y_avg{ind},ind);
%     %p=plot_preliminaries(F,sum(Yabs,1)/length(selected_files1),1);
%     set(p,'MarkerSize',0.5);
%     xlim([-8e3,8e3]);
%     ylim([-0.1 1])
%     %ylim([-2.5 -0.1])
%     %plot_labels('Frequency [Hz]','log(Signal)');
%     plot_labels('Frequency [Hz]','Signal');
%     title(['Pulse Width ' num2str(inputArray(ind)) ' \mus']);
%     pline=plot_yline(1000,2);
%     pline2=plot_yline(4000,2);
%     set(pline,'LineStyle','--');
%     set(pline2,'LineStyle','--');
%     %gifdraw(h,ind,filename,2);
%     filename=['pw_sweep_average_pw_' num2str(inputArray(ind)) '_us.fig'];
%     savefig(h1,[dirname2 filename]);
%     hold off;
%     clear Yabs;
% end
%
% peaks=abs(peaks);
% %% fig6b
%
% peakssec=peaks;
% %if pulse width is larger than 35 we see folding back
% foldind=find(inputArray>35);
% for ind=foldind
%     peakssec(ind)=2*max(F{ind})-peakssec(ind);
% end
%
% peaksfir=peakssec/2;
%
% optim=optimoptions('lsqcurvefit','display','off');
% lb=[-inf -inf];
% ub=[inf inf];
% harfun= @(par,x) par(1)*x./(x+par(2));
% h3=start_fig(3,[2 2]);
% ppeak=plot_preliminaries(inputArray,peakssec/2,1,'noline');
% set(ppeak,'DisplayName','Experimental Harmonic Freq');
%
% par0=[2000,50];
% harfit=lsqcurvefit(harfun,par0,inputArray,peakssec/2,lb,ub,optim);
% xfit=linspace(min(inputArray),max(inputArray),500);
% yfit=harfun(harfit,xfit);
% pfit=plot_preliminaries(xfit,yfit,2,'nomarker');
% set(pfit,'DisplayName','Theory Fit');
%
% tdead=43; %us
% xtheo=linspace(min(inputArray),max(inputArray),500);
% ytheo=xtheo./(xtheo+tdead)/120*1e6; %Hz
% ptheo=plot_preliminaries(xtheo,ytheo,3,'nomarker');
% set(ptheo,'DisplayName','Theoretical Harmonic Freq');
%
% plot_labels('Pulse Width [us]','First Harmonic Freq');
% legend('Location','best');
% savefig(h3,[dirname 'fig6b.fig']);
%
% %% fig6c
%
% for ind=1:length(F)
%     linewidth=400;
%     [lol indfirhar(ind)]=min(abs(F{ind}-peaks(ind)/2));
%     [lol indsechar(ind)]=min(abs(F{ind}-peaks(ind)));
%      [lol indfirharl(ind)]=min(abs(F{ind}-peaks(ind)/2+linewidth));
%     [lol indfirharr(ind)]=min(abs(F{ind}-peaks(ind)/2-linewidth));
%     [lol indsecharl(ind)]=min(abs(F{ind}-peaks(ind)+linewidth));
%     [lol indsecharr(ind)]=min(abs(F{ind}-peaks(ind)-linewidth));
%     intensity1(ind)=sum(y_avg{ind}(indfirharl(ind):indfirharr(ind)));
%     intensity2(ind)=sum(y_avg{ind}(indsecharl(ind):indsecharr(ind)));
% end
%
% h4=start_fig(4,[2 2]);
%
% % p4(1)=plot_preliminaries(inputArray,intensity1,2);
% % p4(2)=plot_preliminaries(inputArray,intensity2,2);
%
% IXlol=find(inputArray==35);
% %p4(3)=plot_preliminaries(inputArray,intensity2./intensity1,1);
% p4(3)=plot_preliminaries(inputArray(1:IXlol),intensity2(1:IXlol)./intensity1(1:IXlol),1);
% plot_labels('Pulse Width [us]','Signal Intensity [au]');
% %legend('First Harmonic','Second Harmonic','Ratio','Location','best');
% savefig(h4,[dirname 'fig6c.fig']);
%
%% fig6d processing
clear;clc;close all

explimits11{1}='2021-05-28-124811_Proteus';
explimits11{2}='2021-05-29-014441_Proteus';
inputArray11=[5:100]*1.695/4; %voltage mVpp multiplied by uT/mVpp
ACfreq=2760;

if contains(pwd,'AjoyLabBlueOak')
    sys=1; %1 for lab PC
else
    sys=0; %laptop
end

if sys==1
    dirname='C:\Users\AjoyLabBlueOak\OneDrive\Documents\MATLAB\Berkeley_Lab\mag paper figs\fig6\';
else
    dirname='D:\OneDrive\Documents\MATLAB\Berkeley_Lab\mag paper figs\fig6\';
end
ex=[min(inputArray11) max(inputArray11)]; %min max points for input

fourier_plot2=false;
decay_plot=false;
intensity_plot=true;



selected_files11= transpose(agilent_read_proteus(explimits11,sys));
%selected_files11=explimits11';
y1 = [1:1:size(selected_files11,1)];

inputArray11 = inputArray11(1:size(selected_files11,1));
for j=1:length(selected_files11)
    [S(j),x2,y2{j},x_truncated,y_truncated{j}] = process_data_proteus(selected_files11{j},sys);
end

clear S x2 y2;



%% fig6d
%truncate to 5s
[v IX3]=min(abs(x_truncated-5));
x_truncated=x_truncated(1:IX3);
for ind=1:length(inputArray11)
    y_truncatedfull{ind}=(y_truncated{ind}(1:IX3));
end

if fourier_plot2
    h=start_fig(3,[2 2]);
    plot_labels('Frequency [Hz]','Signal');
    %pline=plot_yline(inp,3);set(pline,'linestyle','--');
end
for ind=1:length(selected_files11)
    %Smooth data
    baseline = smooth(y_truncatedfull{ind},1000)';
    
    % mean(y)
    y_sub=y_truncatedfull{ind}-baseline;
    
    x_sub2=x_truncated;
    Ts=x_sub2(2)-x_sub2(1);
    Y2=fft(y_sub);
    Y2abs(ind,:)=abs(fftshift(Y2))*Ts;
    Y2base=abs(fftshift(fft(baseline)))*Ts;
    
    n=size(Y2,2);
    fs=1/max(x_sub2);
    %F=(-n/2+1:n/2)*(fs);
    F3=linspace(-1/2/Ts,1/2/Ts,n);
    %[v,b]=min(abs(F-0));
    %S=Y2abs(b);
    
    [v IX3]=min(abs(F3-ACfreq));
    voltintensity(ind)=max(Y2abs(ind,IX3-10:IX3+10));
    
    [v IX3]=min(abs(F3-2*ACfreq));
    voltintensity2(ind)=max(Y2abs(ind,IX3-10:IX3+10));
    
    if fourier_plot2
        %p=plot_preliminaries(F,sum(Yabs,1)/ind,1);
        p=plot_preliminaries(F3,Y2abs(ind,:),1);
        ylim([-1 15]);
        xlim([0 max(F3)]);
        hold off
        set(p,'MarkerSize',0.5);
        
        
        title([num2str(inputArray11(ind)) ' mVpp']);
        ylim([-1 6])
        %xlim([5515 5525]);
        xlim([2755 2765]);
        gifdraw(h,ind,'fourier_plot_for_paperfirsthar.gif',0.1)
    end
end


%%
if intensity_plot
    h5=start_fig(11,[2 2]);
    
    %     for ind=1:length(selected_files11)
    %         [v IX3]=min(abs(F-ACfreq));
    %
    %         voltintensity(ind)=Yabs(ind,IX3);
    %         [v IX3]=min(abs(F-2*ACfreq));
    %
    %         voltintensity2(ind)=Yabs(ind,IX3);
    %     end
    
    %     voltintensity=smooth(voltintensity,100);
    %     voltintensity2=smooth(voltintensity2,100);
    
    p1=plot_preliminaries(inputArray11,voltintensity,1,'noline');
    set(p1,'DisplayName','First Harmonic');
    p2=plot_preliminaries(inputArray11,voltintensity2,2,'noline');
    set(p2,'DisplayName','Second Harmonic');
    %plot_labels('Voltage (mVpp)','Signal Intensity');
    plot_labels('Field Strength |B_{AC}| (\muT)','Signal Intensity');
    
    
    linfit1=polyfit(inputArray11,voltintensity,1);
    linfit2=polyfit(inputArray11,voltintensity2,2);
    xfit=linspace(min(inputArray11),max(inputArray11),500);
    yfit1=polyval(linfit1,xfit);
    yfit2=polyval(linfit2,xfit);
    plot_preliminaries(xfit,yfit1,1,'nomarker');
    plot_preliminaries(xfit,yfit2,2,'nomarker');
    
    legend([p1 p2],'Location','best');
    
    savefig(h5,[dirname 'fig6d.fig']);
end

%% sensitivity calculation
wingdata=Y2abs(1:11,end-100:end)';
noisedev=std(wingdata);
noisemean=mean(wingdata);

Bfields=inputArray11(1:11)*1695/4; %nT
voltinttrun=voltintensity(1:11)-noisemean;
[linfittrun,Sfit]=polyfit(Bfields,voltinttrun,1);
xfittrun=linspace(0,max(Bfields),500);
[yfittrun]=polyval(linfittrun,xfittrun,Sfit);

start_fig(55,[2 2]);
plot(xfittrun,yfittrun,'DisplayName','Line fit');
hold on
% plot(xfittrun,yfittrun,'r--','DisplayName','68% Confidence Interval');
% plot(xfittrun,yfittrun,'r--','DisplayName','68% Confidence Interval');
errorbar(Bfields,voltinttrun,noisedev,'DisplayName','Data');
% xlim([0 17]);
% ylim([0 1.25]);
xlabel('Magnetic Field [nT]');
legend();
hold off

%noiselev=noisemean+noisedev;

% minfield=(noiselev-linfittrun(2))/linfittrun(1);
% minfielddev=mean(delta);
sens=noisedev/linfittrun(1); %nT


%% fig6e processing

explimits12{1}='2021-04-30-191102_Proteus';
explimits12{2}='2021-05-01-020329_Proteus';
% explimits12{1}='2021-05-03-134435_Proteus';
% explimits12{2}='2021-05-03-160309_Proteus';
inputArray12=[5:90]*1.695/4; %voltage mVpp
ACfreq=1750;

ex=[min(inputArray12) max(inputArray12)]; %min max points for input

fourier_plot2=false;
decay_plot=false;
intensity_plot=true;

if contains(pwd,'AjoyLabBlueOak')
    sys=1; %1 for lab PC
else
    sys=0; %laptop
end

selected_files12= transpose(agilent_read_proteus(explimits12,sys));
%selected_files12=explimits12';
y1 = [1:1:size(selected_files12,1)];

inputArray12 = inputArray12(1:size(selected_files12,1));
for j=1:length(selected_files12)
    [S(j),x2,y2{j},x_truncated,y_truncated{j}] = process_data_proteus(selected_files12{j},sys);
end

clear S x2 y2;

%% fig6e

[v IX3]=min(abs(x_truncated-5));
x_truncated=x_truncated(1:IX3);
for ind=1:length(inputArray12)
    y_truncatedfull{ind}=(y_truncated{ind}(1:IX3));
end

if fourier_plot2
    h=start_fig(3,[2 2]);
    plot_labels('Frequency [Hz]','Signal');
    %pline=plot_yline(inp,3);set(pline,'linestyle','--');
end
for ind=1:length(selected_files12)
    %Smooth data
    baseline = smooth(y_truncatedfull{ind},1000)';
    
    % mean(y)
    y_sub11=y_truncatedfull{ind}-baseline;
    
    x_sub21=x_truncated;
    Ts=x_sub21(2)-x_sub21(1);
    Y21=fft(y_sub11);
    Y21abs(ind,:)=abs(fftshift(Y21))*Ts;
    Y21base=abs(fftshift(fft(baseline)))*Ts;
    
    n=size(Y21,2);
    fs=1/max(x_sub21);
    %F=(-n/2+1:n/2)*(fs);
    F31=linspace(-1/2/Ts,1/2/Ts,n);
    %[v,b]=min(abs(F-0));
    %S=Y2abs(b);
    
    [v IX31]=min(abs(F31-ACfreq));
    voltintensity11(ind)=max(Y21abs(ind,IX31-10:IX31+10));
    
    [v IX31]=min(abs(F31-2*ACfreq));
    voltintensity21(ind)=max(Y21abs(ind,IX31-10:IX31+10));
    
    if fourier_plot2
        %p=plot_preliminaries(F,sum(Yabs,1)/ind,1);
        p=plot_preliminaries(F31,Y21abs(ind,:),1);
        ylim([-1 15]);
        xlim([0 max(F31)]);
        hold off
        set(p,'MarkerSize',0.5);
        
        
        title([num2str(inputArray12(ind)) ' mVpp']);
        ylim([-1 6])
        %xlim([5515 5525]);
        xlim([2755 2765]);
        %gifdraw(h,ind,'fourier_plot_for_paperfirsthar.gif',0.1)
    end
end


%%
if intensity_plot
    h6=start_fig(12,[2 2]);
    
    %     for ind=1:length(selected_files12)
    %         [v IX31]=min(abs(F-ACfreq));
    %
    %         voltintensity(ind)=Yabs(ind,IX31);
    %         [v IX31]=min(abs(F-2*ACfreq));
    %
    %         voltintensity2(ind)=Yabs(ind,IX31);
    %     end
    
    %     voltintensity=smooth(voltintensity,100);
    %     voltintensity2=smooth(voltintensity2,100);
    
    p11=plot_preliminaries(inputArray12,voltintensity11,1,'noline');
    set(p11,'DisplayName','First Harmonic');
    p21=plot_preliminaries(inputArray12,voltintensity21,2,'noline');
    set(p21,'DisplayName','Second Harmonic');
    %plot_labels('Voltage (mVpp)','Signal Intensity');
    plot_labels('Field Strength |B_{AC}| (\muT)','Signal Intensity');
    
    
    linfit11=polyfit(inputArray12,voltintensity11,1);
    linfit21=polyfit(inputArray12,voltintensity21,2);
    xfit=linspace(min(inputArray12),max(inputArray12),500);
    yfit11=polyval(linfit11,xfit);
    yfit21=polyval(linfit21,xfit);
    plot_preliminaries(xfit,yfit11,1,'nomarker');
    plot_preliminaries(xfit,yfit21,2,'nomarker');
    
    legend([p11 p21],'Location','best');
    %dirname='C:\Users\AjoyLabBlueOak\OneDrive\Documents\MATLAB\Berkeley_Lab\paper figs\fig6\';
    savefig(h6,[dirname 'fig6e.fig']);
    saveas(h6,[dirname 'fig6e.svg'],'svg');
end
