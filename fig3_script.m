clear;clc;close all;

explimits{1}='2021-05-30-013506_Proteus';
explimits{2}='2021-05-30-111737_Proteus';
inputArray1=[1:100]; %iterations



if contains(pwd,'AjoyLabBlueOak')
    sys=1; %1 for lab PC
    dirname='C:\Users\AjoyLabBlueOak\OneDrive\Documents\MATLAB\Berkeley_Lab\mag paper figs\fig3\';
else
    sys=0; %laptop
    dirname='D:\OneDrive\Documents\MATLAB\Berkeley_Lab\mag paper figs\fig3\';
end

selected_files= transpose(agilent_read_proteus(explimits,sys));

for j=1:size(selected_files,1)
    [S(j),xfull{j},yfull{j},x,y{j},baseline{j}] = process_data_proteus(selected_files{j},sys);
end


%%
osc=y{1}-baseline{1};
h2=start_fig(2,[2 2]); %inset graph full oscillations
p2=plot_preliminaries(x,osc,1,'nomarker'); 
%set(p2,'MarkerSize',3);
plot_labels('Time [s]','Signal [au]');
yline(0);
xlim([0 20]);
ylim([-10 10]);
savefig(h2,[dirname 'fig3a.fig']);
saveas(h2,[dirname 'fig3a.svg'],'svg');


L=300; %truncation index
h1=start_fig(1,[4 1]); %truncated version
x_trunc=x(1:L);
osc_trunc=osc(1:L);
mid=(max(osc_trunc)+min(osc_trunc))/2;
rang=(max(osc_trunc)-min(osc_trunc))/2;

p1=plot_preliminaries(x_trunc*1e3,osc_trunc,1);
set(p1,'MarkerSize',3);

% SF=splinefit(x_trunc,osc_trunc,295,3);
% yspl=ppval(SF,x_trunc);
% pspl=plot_preliminaries(x_trunc*1e3,yspl,2,'nomarker');

plot_labels('Time [ms]','Signal [au]');
yline(0);
%yline(mean(osc_trunc));
ylim([mid-1.1*rang mid+1.1*rang]);
xlim([min(x_trunc*1e3) max(x_trunc*1e3)]);
savefig(h1,[dirname 'fig3a_inset.fig']);
saveas(h1,[dirname 'fig3a_inset.svg'],'svg');

%% FT
osctot=0;
for ind=1:length(inputArray1)
    osctot=osctot+y{ind}-baseline{ind};
end
osctot=osctot/length(inputArray1);
Ts=x(2)-x(1);
Yosc=abs(fftshift(fft(osctot)))*Ts;
F=fftfreq(x,true);

[h3,ax3]=start_fig(3,[6 2]);
p3=plot_preliminaries(F,Yosc,2,'nomarker');
%set(p3,'MarkerSize',2);
xline(2760);
xline(2760*2);
plot_labels('Frequency [Hz]','Signal [au]');
xlim([0 max(F)]);
ylim([-0.2 4]);
ax3.XTick=[ax3.XTick max(F)];
xtickformat(ax3,'%.0f');
savefig(h3,[dirname 'fig3b_mainFT.fig']);
saveas(h3,[dirname 'fig3b_mainFT.svg'],'svg');

[h4,ax4]=start_fig(4,[2 2]); %noise inset
p4=plot_preliminaries(F,Yosc,2);
set(p4,'MarkerSize',3);
plot_labels('Frequency [Hz]','Signal [au]');
%xlim([2755 2761]);
xlim([2055 2057]);
ylim([0 0.003]);
ax4.XTick=[2055:2057];

saveas(h4,[dirname 'fig3b_inset_noise.svg'],'svg');

h5=start_fig(5,[2 2]); %second harmonic inset
p5=plot_preliminaries(F,Yosc,2);
set(p5,'MarkerSize',5);
yline(3.43423/2);
plot_labels('Frequency [Hz]','Signal [au]');
rang=1;
xlim([5520-rang 5520+rang]);
ylim([0 4]);
savefig(h5,[dirname 'fig3b_inset_secondharmonic.fig']);
saveas(h5,[dirname 'fig3b_inset_secondharmonic.svg'],'svg');

h6=start_fig(6,[2 2]); %first harmonic inset
p6=plot_preliminaries(F,Yosc,2);
set(p6,'MarkerSize',3);
yline(0.535/2);
plot_labels('Frequency [Hz]','Signal [au]');
rang=1;
xlim([2760-rang 2760+rang]);
ylim([0 0.7]);
savefig(h6,[dirname 'fig3b_inset_firstharmonic.fig']);
saveas(h6,[dirname 'fig3b_inset_firstharmonic.svg'],'svg');

%% fig3c loglog graph
h13=start_fig(13,[6 2]);
p13=plot_preliminaries(F,Yosc,3,'nomarker');
%set(p13,'MarkerSize',3);
plot_labels('Frequency [Hz]','Signal [au]');
set(gca,'Yscale','log','Xscale','log');
savefig(h13,[dirname 'fig3c_mainFT.fig']);
saveas(h13,[dirname 'fig3c_mainFT.svg'],'svg');

h14=start_fig(14,[2 2]); %noise inset
p14=plot_preliminaries(F,Yosc,3);
set(p14,'MarkerSize',3);
plot_labels('Frequency [Hz]','Signal [au]');
%xlim([2755 2761]);
xlim([2055 2057]);
ylim([1e-4 0.005]);
set(gca,'Yscale','log','Xscale','log');
savefig(h14,[dirname 'fig3c_inset_noise.fig']);
saveas(h14,[dirname 'fig3c_inset_noise.svg'],'svg');

h15=start_fig(15,[2 2]); %second harmonic inset
p15=plot_preliminaries(F,Yosc,3);
set(p15,'MarkerSize',3);
yline(3.43423/2);
plot_labels('Frequency [Hz]','Signal [au]');
rang=2;
xlim([5520-rang 5520+rang]);
ylim([1e-2 4]);
set(gca,'Yscale','log','Xscale','log');
savefig(h15,[dirname 'fig3c_inset_secondharmonic.fig']);
saveas(h15,[dirname 'fig3c_inset_secondharmonic.svg'],'svg');

h16=start_fig(16,[2 2]); %first harmonic inset
p16=plot_preliminaries(F,Yosc,3);
set(p16,'MarkerSize',3);
yline(0.535/2);
plot_labels('Frequency [Hz]','Signal [au]');
rang=1;
xlim([2760-rang 2760+rang]);
ylim([5e-3 0.6]);
set(gca,'Yscale','log','Xscale','log');
savefig(h16,[dirname 'fig3c_inset_firstharmonic.fig']);
saveas(h16,[dirname 'fig3c_inset_firstharmonic.svg'],'svg');

