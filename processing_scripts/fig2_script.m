%% paper figure: frequency sweep intensity
clear;clc;close all;

explimits1{1}='2021-05-23-233035_Proteus';
explimits1{2}='2021-05-24-052241_Proteus';
inputArray1=[2600:5:2900]; %frequency

explimits2{1}='2021-05-24-234222_Proteus';
explimits2{2}='2021-05-25-055753_Proteus';
inputArray2=[2880:5:3200];

explimits3{1}='2021-05-25-160643_Proteus';
explimits3{2}='2021-05-25-210046_Proteus';
inputArray3=[2000:40:4000];

explimits4{1}='2021-05-26-141927_Proteus';
explimits4{2}='2021-05-26-161108_Proteus';
inputArray4=[100:100:2000];

explimits5{1}='2021-05-26-173924_Proteus';
explimits5{2}='2021-05-26-201813_Proteus';
inputArray5=[4100:100:6800];

explimits6{1}='2021-06-01-150601_Proteus';
explimits6{2}='2021-06-01-150601_Proteus';
inputArray6=[0];



%ex=[min(inputArray1) max(inputArray1)]; %min max points for input
fourier_plot=true;
decay_plot=true;
T2_plot=true;

if contains(pwd,'AjoyLabBlueOak')
    sys=1; %1 for lab PC
else
    sys=0; %laptop
end

selected_files1= transpose(agilent_read_proteus(explimits1,sys));
%selected_files1=explimits1';
y1 = [1:1:size(selected_files1,1)];

inputArray1 = inputArray1(1:size(selected_files1,1));

selected_files2= transpose(agilent_read_proteus(explimits2,sys));
%selected_files1=explimits1';
yx2 = [1:1:size(selected_files2,1)];
inputArray2 = inputArray2(1:size(selected_files2,1));
IX=find(~(inputArray2<=2900));
inputArray2=inputArray2(IX);
selected_files2=selected_files2(IX);

selected_files3= transpose(agilent_read_proteus(explimits3,sys));
%selected_files1=explimits3';
yx3 = [1:1:size(selected_files3,1)];
inputArray3 = inputArray3(1:size(selected_files3,1));
IX=find(~((inputArray3>=2600)&(inputArray3<=3200)|inputArray3==2000));
inputArray3=inputArray3(IX);
selected_files3=selected_files3(IX);

selected_files4= transpose(agilent_read_proteus(explimits4,sys));
%selected_files4=explimits4';
yx4 = [1:1:size(selected_files4,1)];
inputArray4 = inputArray4(1:size(selected_files4,1));

selected_files5= transpose(agilent_read_proteus(explimits5,sys));
%selected_files5=explimits5';
yx5 = [1:1:size(selected_files5,1)];
inputArray5 = inputArray5(1:size(selected_files5,1));

selected_files6= transpose(agilent_read_proteus(explimits6,sys));
%selected_files6=explimits6';
yx5 = [1:1:size(selected_files6,1)];
inputArray6 = inputArray6(1:size(selected_files6,1));


for j=1:size(selected_files1,1)
    [S(j),x{j},y{j},x_truncated,y_truncated{j}] = process_data_proteus(selected_files1{j},sys);
end

for j=1:size(selected_files2,1)
    [S2(j),x2{j},y2{j},x_truncated2,y_truncated2{j}] = process_data_proteus(selected_files2{j},sys);
end

for j=1:size(selected_files3,1)
    [S3(j),x3{j},y3{j},x_truncated3,y_truncated3{j}] = process_data_proteus(selected_files3{j},sys);
end

for j=1:size(selected_files4,1)
    [S4(j),x4{j},y4{j},x_truncated4,y_truncated4{j}] = process_data_proteus(selected_files4{j},sys);
end

for j=1:size(selected_files5,1)
    [S5(j),x5{j},y5{j},x_truncated5,y_truncated5{j}] = process_data_proteus(selected_files5{j},sys);
end

for j=1:size(selected_files6,1)
    [S6(j),x6{j},y6{j},x_truncated6,y_truncated6{j}] = process_data_proteus(selected_files6{j},sys);
end

inputArray=[inputArray1 inputArray2 inputArray3 inputArray4 inputArray5 inputArray6];
[inputArray IX]=sort(inputArray);
y_truncatedfull=[y_truncated y_truncated2 y_truncated3 y_truncated4 y_truncated5 y_truncated6];
y_truncatedfull=y_truncatedfull(IX);
% Sfull=[S S2 S3];
% Sfull=Sfull(IX);
selected_files=[selected_files1' selected_files2' selected_files3' selected_files4' selected_files5' selected_files6'];
selected_files=selected_files(IX);

clear S S2 S3 S4 S5 S6;


%truncate to 5s
[v IX]=min(abs(x_truncated-5));
x_truncated=x_truncated(1:IX);
for ind=1:length(inputArray)
    y_truncatedfull{ind}=(y_truncatedfull{ind}(1:IX));
end
%%
% if fourier_plot
%     h=start_fig(3,[4 2]);
% end
%pline=plot_yline(inp,3);set(pline,'linestyle','--');
for ind=1:length(selected_files)
    %Smooth data
    baseline = smooth(y_truncatedfull{ind},1000)';
    
    % mean(y)
    y_sub=y_truncatedfull{ind}-baseline;
    
    x_sub=x_truncated;
    Ts=x_sub(2)-x_sub(1);
    Y=fft(y_sub);
    Yabs(ind,:)=abs(fftshift(Y))*Ts;
    Ybase=abs(fftshift(fft(baseline)))*Ts;
    
    n=length(Y);
    fs=1/max(x_sub);
    %F=(-n/2+1:n/2)*(fs);
    F=linspace(-1/2/Ts,1/2/Ts,n);
    %[v,b]=min(abs(F-0));
    %S=Yabs(b);
    
%     if fourier_plot
%         p=plot_preliminaries(F,Yabs(ind,:),1);
%         %p=plot_preliminaries(F,Ybase,1);
%         
%         hold off
%         set(p,'MarkerSize',0.5);
%         
%         plot_labels('Frequency [Hz]','Signal');
%         title(['Frequency ' num2str(inputArray(ind)) ' Hz']);
%         ylim([-0.1 8])
%         xlim([0 7000]);
%         %gifdraw(h,ind,'fourier_plot.gif',0.1)
%     end
end

if decay_plot
    h=start_fig(2,[2 2]);
    plot_labels('Time [s]', 'Signal');
end
for j=1:length(selected_files)
    y_truncatedfull{j}=y_truncatedfull{j}/max(y_truncatedfull{j}); %normalization
    if decay_plot
        p(j)=plot_preliminaries(x_truncated,(y_truncatedfull{j}),j,'nomarker'); 
        hold off
        %ylim([0 1]);
        title([num2str(inputArray(j)) ' Hz']);
        
        % set(p(j),'MarkerSize',0.5);
        %p2(j)=plot_preliminaries(x_truncated_5{j},baseline,2,'nomarker');
        %gifdraw(h,j,'decay_plot.gif',0.1);
        ypro=log10(y_truncatedfull{j});
        xpro=sqrt(x_truncated);
        linfit(j,:)=polyfit(xpro,ypro,1);
    end
end

%plot important points fig2b
points=[0 2885 2800 5000];
%points=inputArray(34:116)
if sys==1
    dirname='C:\Users\AjoyLabBlueOak\OneDrive\Documents\MATLAB\Berkeley_Lab\paper figs\fig2\';
else
    dirname='D:\OneDrive\Documents\MATLAB\Berkeley_Lab\mag paper figs\fig2\';
end
h11=start_fig(11,[2 2]);
for ind=1:length(points)
    
    IX=find(inputArray==points(ind));
    ppoints(ind)=plot_preliminaries(x_truncated,smooth(y_truncatedfull{IX},100)',ind,'nomarker');
%     drawnow;
%     hold off;
    if inputArray(IX)==0
        set(ppoints(ind),'DisplayName','DC');
    else
        set(ppoints(ind),'DisplayName',[num2str(inputArray(IX)) ' Hz']);
    end
    
    
end
ylim([0 1]);
legend('Location','best'); 
plot_labels('Time [s]','Normalized Signal');
savefig(h11,[dirname 'fig2b_decay.fig']);

%plot important points in sqrt-log axis fig2c
h21=start_fig(21,[2 2]);
set(gca, 'YScale', 'log')
for ind=1:length(points)
    
    IX=find(inputArray==points(ind));
    ppoints(ind)=plot_preliminaries(sqrt(x_truncated),(smooth(y_truncatedfull{IX},100)),ind,'nomarker');
    xlim([0 sqrt(5)]);
    if inputArray(IX)==0
        set(ppoints(ind),'DisplayName','DC');
    else
        set(ppoints(ind),'DisplayName',[num2str(inputArray(IX)) ' Hz']);
    end
    
    
end

set(gca, 'YScale', 'log')
legend('Location','best'); 
plot_labels('sqrt(Time) [s]','Normalized Signal');
savefig(h21,[dirname 'fig2c_decaysqrtlog.fig']);


if T2_plot
    for ind=1:length(selected_files)
        Sfull(ind)=sum(y_truncatedfull{ind})*Ts; %integrated signal
    end
   
    h5=start_fig(5,[4 1]);
    %gaussian fitting
    %Ssmooth=smooth(Sfull,50)';
    %     gaussfun= @(par,x) par(1)*exp(-(x-par(4)).^2/par(2))+par(3);
    %     lorefun= @(par,x) par(1)./(1+2/par(3)*(x-par(2)).^2)+par(4) %par1 Amplitude; par2 center; par3 FWHM
    %     %par0=[1.8 2900 100 4]; % for lorentzian
    %     par0=[-4 100 4 2900]; %for gaussian
    %     fitgau=lsqcurvefit(gaussfun,par0,inputArray,Ssmooth);
    %     %fitlor=lsqcurvefit(lorefun,par0,inputArray,Ssmooth);
    %     xfit=linspace(min(inputArray),max(inputArray),5000);
    %     %fitgau(1)=-1.2;
    %     yfit=gaussfun(fitgau,xfit);
    %     %yfit=lorefun(fitlor,xfit);
    %     linewidth=sqrt(fitgau(2)/2)*2*sqrt(2*log(2));
    
   
    pinteg=plot_preliminaries(inputArray,Sfull,1,'noline');
    set(pinteg,'MarkerSize',3);
    %     pfit=plot_preliminaries(xfit,yfit,2,'nomarker'); %gaussian fit
    %     FWHM1=fitgau(4)+linewidth/2;
    %     FWHM2=fitgau(4)-linewidth/2;
    %     yFWHM=[gaussfun(fitgau,FWHM1) gaussfun(fitgau,FWHM2)];
    %     pFWHM=plot_preliminaries([FWHM1 FWHM2],yFWHM,3,'noline');
    
    %spline fit
    xspl=linspace(min(inputArray),max(inputArray),5000);
    SF=splinefit(inputArray,Sfull,8,3);
    yy = ppval(SF,xspl);
    pspl=plot_preliminaries(xspl,yy,2,'nomarker');
    
    %dirname='C:\Users\AjoyLabBlueOak\OneDrive\Documents\MATLAB\Berkeley_Lab\paper figs\fig2\';
    plot_labels('Frequency [Hz]','Integrated Signal');
    savefig([dirname 'fig2a.fig']);
    
    % legend('end','mid','quarter');
end

%colormap code
% xlim=find(freq_8>16.25e3 & freq_8<33.75e3);
% dnp_combined_1=flatten_spectrum(bs_spec_dnp8,abs(enhancementFactor8),xlim);
% start_fig(2,[1 2]);
% y=inputArray8;x=freq_8(xlim);
% pcolor(x,y,(dnp_combined_1));
% grid on;
% set(gca,'layer','top');
% set(gca,'Ydir','reverse')
% shading interp;
% colormap(NegativeEnhancingColormap(128, [min(dnp_combined_1(:)) max(dnp_combined_1(:))], ...
%         [0 0 1], [1 0 0], 1));
% colorbar;
% title('1%');

%% fig2d

truncation=linspace(0.1,5,100);
clear S S2 S3 S4 S5 S6 linewidth;
if fourier_plot
    h=start_fig(1,[2 2]);
end
for ind=1:length(truncation)
    [v IX]=min(abs(x_truncated-truncation(ind)));
    x_truncatedlol=x_truncated(1:IX);
    
    for j=1:length(inputArray)
        y_truncatedfulllol{j}=(y_truncatedfull{j}(1:IX));
        %plot_preliminaries(x_truncatedlol,y_truncatedfulllol{j},j,'nomarker');
        hold off
        Ts=x_truncatedlol(2)-x_truncatedlol(1);
        S(j)=sum(y_truncatedfulllol{j})*Ts; %integrated signal up to truncation point
    end
    Ssmooth=smooth(S,100)';
    if fourier_plot
    plot_preliminaries(inputArray,S,1);
    end
    
    
    gaussfun= @(par,x) par(1)*exp(-(x-par(4)).^2/par(2))+par(3);
    lorefun= @(par,x) par(1)./(1+2/par(3)*(x-par(2)).^2)+par(4); %par1 Amplitude; par2 center; par3 FWHM
    optim=optimoptions('lsqcurvefit','Display','off');
    %par0=[1.8 2900 100 4]; % for lorentzian
    par0=[0 100 4 2900]; %for gaussian
    lb=[-inf,50,-inf,-inf];
    ub=[inf,inf,inf,inf];
    fitgau=lsqcurvefit(gaussfun,par0,inputArray,S,lb,ub,optim);
    %fitlor=lsqcurvefit(lorefun,par0,inputArray,Ssmooth,lb,ub,optim);
    xfit=linspace(min(inputArray),max(inputArray),5000);
    yfit=gaussfun(fitgau,xfit);
    %yfit=lorefun(fitlor,xfit);
    if fourier_plot
        plot_preliminaries(xfit,yfit,2,'nomarker');
        xlim([2400 3400]);
        gifdraw(h,ind,'linewidth_graph.gif',0.2);
    end
    linewidth(ind)=sqrt(fitgau(2)/2)*2*sqrt(2*log(2));
end
%%
h2=start_fig(2,[2 2]);
p2=plot_preliminaries(truncation/73e-6/1e4,linewidth,1);

xlim([0 max(truncation/73e-6/1e4)]);
plot_labels('No. of Pulses (x10^4)','Linewidth [Hz]');
set(p2,'Color','m','MarkerFaceColor','m','MarkerEdgeColor','m');
%savefig(h2,[dirname 'fig2d.fig']);
%saveas(h2,[dirname 'fig2d.pdf'],'pdf');
saveas(h2,[dirname 'fig2d.svg'],'svg');

