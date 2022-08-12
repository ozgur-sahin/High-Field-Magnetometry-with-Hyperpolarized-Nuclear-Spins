clear;clc;close all;

explimits1{1}='2021-05-30-013506_Proteus';
explimits1{2}='2021-05-30-111737_Proteus';
inputArray=[1:100];
%explimits1{2}=explimits1{1};
ex=[min(inputArray) max(inputArray)]; %min max points for input
fourier_plot=false;
decay_plot=false;
intensity_plot=false;
normalize=false;

if contains(pwd,'AjoyLabBlueOak')
    sys=1; %1 for lab PC
else
    sys=0; %laptop
end

selected_files1= transpose(agilent_read_proteus(explimits1,sys));
%selected_files1=explimits1';
y1 = [1:1:size(selected_files1,1)];

inputArray = inputArray(1:size(selected_files1,1));
for j=1:size(selected_files1,1)
    [S_5(j),x_5{j},y_5{j},x_truncated{j},y_truncated{j}] = process_data_proteus(selected_files1{j},sys);   
end
%%
if fourier_plot
    h=start_fig(3,[4 2]);
    %pline=plot_yline(inp,3);set(pline,'linestyle','--');
end
for ind=1:length(selected_files1)
    %Smooth data
    baseline = smooth(y_truncated{ind},1000)';
    
    % mean(y)
    y_sub=y_truncated{ind}-baseline;
    
    x_sub=x_truncated{ind};
    Ts=x_sub(2)-x_sub(1);
    Y=fft(y_sub);
    Yabs(ind,:)=abs(fftshift(Y))*Ts;
    Ybase=abs(fftshift(fft(baseline)))*Ts;
    
    n=size(Y,2);
    fs=1/max(x_sub);
    %F=(-n/2+1:n/2)*(fs);
    F=linspace(-1/2/Ts,1/2/Ts,n);
    %[v,b]=min(abs(F-0));
    %S=Yabs(b);
    
    if fourier_plot
        p=plot_preliminaries(F,sum(Yabs,1)/ind,1);
        %p=plot_preliminaries(F,Ybase,1);
        
        hold off
        set(p,'MarkerSize',0.5);
        
        plot_labels('Frequency [Hz]','Signal');
        title(['Iteration No: ' num2str(inputArray(ind))]);
        ylim([-0.1 1])
        %gifdraw(h,ind,'fourier_plot.gif',0.1)
    end
end
if fourier_plot
    h=start_fig(3,[4 2]);
    p=plot_preliminaries(F,(sum(Yabs,1)/length(selected_files1)),1);
    %p=plot_preliminaries(F,sum(Yabs,1)/length(selected_files1),1);
    set(p,'MarkerSize',0.5);
    %ylim([-0.1 0.5])
    %ylim([-2.5 -0.1])
    plot_labels('Frequency [Hz]','Signal');
    %plot_labels('Frequency [Hz]','Signal');
    pline=plot_yline(1000,2);
    pline2=plot_yline(4000,2);
    set(pline,'LineStyle','--');
    set(pline2,'LineStyle','--');
end


if decay_plot
    h=start_fig(2,[2 2]);
    plot_labels('Time [s]', 'Signal');
    y_avg=0;
end
for j=1:length(selected_files1)
    if normalize
        y_truncated{j}=y_truncated{j}/max(y_truncated{j});
    end
    if decay_plot
        p(j)=plot_preliminaries(x_truncated{j},(y_truncated{j}),1,'nomarker'); %not normalized
        hold off
        %ylim([0 1]);
        
        % set(p(j),'MarkerSize',0.5);
        %p2(j)=plot_preliminaries(x_truncated{j},baseline,2,'nomarker');
        y_avg=y_avg+y_truncated{j};
        gifdraw(h,j,'decay_plot.gif',0.3);
        ypro=log10(y_truncated{j});
        xpro=sqrt(x_truncated{j});
        linfit(j,:)=polyfit(xpro,ypro,1);
    end
end
%     y_avg=y_avg/length(selected_files1);
%     h5=start_fig(4,[2 2]);
%     plot_preliminaries(x_truncated{1},y_avg,1,'nomarker');

%last decay curve (the one going on the paper)
if sys==1
    A=load('C:\Users\AjoyLabBlueOak\Documents\MATLAB\lab data\2021-05-30-152757_Proteus');%no AC field run
    A=load('C:\Users\AjoyLabBlueOak\Documents\MATLAB\lab data\2021-06-01-150601_Proteus');%DC field run
else
    A=load('D:\OneDrive\Documents\MATLAB\lab data\2021-05-30-152757_Proteus');%no AC field run
    A=load('D:\OneDrive\Documents\MATLAB\lab data\2021-06-01-150601_Proteus');%DC field run
end

if normalize
    noAC=A.pulseAmp/max(A.pulseAmp);
else
    noAC=A.pulseAmp;
end
noACx=A.time_axis;

    %%
    dirName='D:\OneDrive\Documents\MATLAB\Berkeley_Lab\mag paper figs\fig1b\';
    h4=start_fig(4,[6 2]);
    p41=plot_preliminaries(x_truncated{1},y_truncated{end},1,'nomarker');
    set(p41,'DisplayName','Resonant Frequency');
    p42=plot_preliminaries(noACx,noAC,5,'nomarker');
    set(p42,'DisplayName','No AC Field','LineStyle','--');
    plot_labels('Time [s]','Signal');
    title('Decay');
    xlim([0 max(x_truncated{1})]);
    legend('Location','best');
    savefig(h4,[dirName 'decay_signal.fig']);
    
    %%
    h5=start_fig(5,[6 2]);
    baseline=smooth(y_truncated{end},100)';
    p51=plot_preliminaries(x_truncated{1},baseline,2,'nomarker');
    set(p51,'DisplayName','Resonant Frequency')
    p52=plot_preliminaries(noACx,smooth(noAC,100)',5,'nomarker');
    set(p52,'DisplayName','No AC Field','LineStyle','--');
    x88=linspace(0,20,500);
    T2=31;
%     y88=exp(-sqrt(x88/T2));
%     p88=plot_preliminaries(x88,y88,4);
%     set(p88,'DisplayName','lol');
    legend('Location','best');
    title('Baseline');
    plot_labels('Time [s]','Signal');
    xlim([0 max(x_truncated{1})]);
    if normalize
        ylim([0 1]);
    end
    savefig(h5,[dirName 'baseline.fig']);
    
    ypro=log10(smooth(noAC,100)');
    xpro=sqrt(noACx);
    ypro2=log10(y88);
    xpro2=sqrt(x88);
    linfit=polyfit(xpro,ypro,1);
%     h8=start_fig(8,[2 2]);
%     p8=plot_preliminaries(xpro,ypro,1);
%     p9=plot_preliminaries(xpro2,ypro2,2);
%     T22=1/(linfit(1)*log(10))^2;
%     
    %%
    h6=start_fig(6,[6 2]);
    plot_preliminaries(x_truncated{1},y_truncated{end}-baseline,3,'nomarker');
    yline(0,'k--');
    plot_labels('Time [s]','Signal');
    title('Baseline Subtracted');
    xlim([0 max(x_truncated{1})]);
    savefig(h6,[dirName 'oscillations.fig']);
    
    %%
    h7=start_fig(7,[2 2]);
    plot_preliminaries(x_truncated{1}(1:100),y_truncated{end}(1:100)-baseline(1:100),3,'nomarker');
    yline(0,'k--');
    plot_labels('Time [s]','Signal');
    %title('Baseline Subtracted');
    xlim([x_truncated{1}(1) x_truncated{1}(100)]);
    savefig(h7,[dirName 'oscillations_inset.fig']);


if intensity_plot
    [lol IX]=min(abs(F-inp));
    intensity=Yabs(IX);
    
    if 2*inp<max(F)
        [lol IX]=min(abs(F-2*inp));
        intensity2=Yabs(IX);
    else
        [lol IX]=min(abs(F-(2*max(F)-2*inp)));
        intensity2=Yabs(IX);
        h=start_fig(1,[2 2]);
        p=plot_preliminaries(inputArray,S,1,'noline');
        %     linfit=polyfit(inputArray,intensity,1);
        %     xfit=linspace(min(inputArray),max(inputArray),500);
        %     yfit=polyval(linfit,xfit);
        set(p,'DisplayName','First Harmonic');
        %     pfit=plot_preliminaries(xfit,yfit,1,'nomarker');
        plot_labels('Pulse Width [\mus]','Signal');
    end
end
