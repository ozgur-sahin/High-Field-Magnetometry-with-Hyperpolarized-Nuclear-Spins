%% fig7 processing
clear;clc;close all;

if contains(pwd,'AjoyLabBlueOak')
    sys=1; %1 for lab PC
else %for laptop
    sys=0;
end

if sys==1
    dirname='C:\Users\AjoyLabBlueOak\OneDrive\Documents\MATLAB\Berkeley_Lab\mag paper figs\fig7\';
else
    dirname='D:\OneDrive\Documents\MATLAB\Berkeley_Lab\mag paper figs\fig7';
end

explimits11{1}='2021-06-02-204801_Proteus';
explimits11{2}='2021-06-03-053149_Proteus';
inputArray11=[1:43]; %no of iterations

if contains(pwd,'AjoyLabBlueOak')
    sys=1; %1 for lab PC
else
    sys=0; %laptop
end

selected_files111= transpose(agilent_read_proteus(explimits11,sys));
selected_files11=selected_files111(1:length(inputArray11));
selected_files12=selected_files111(length(inputArray11)+1:2*length(inputArray11));

for j=1:length(selected_files11)
    [S21(j),x21full{j},y21full{j},x21{j},y21{j}] = process_data_proteus(selected_files11{j},sys);
end

for j=1:length(selected_files12)
    [S21(j),x21full{j},y21full{j},x22{j},y22{j}] = process_data_proteus(selected_files12{j},sys);
end

clear S21 x21full y21full selected_files111;

%% fig7 plotting
h31=start_fig(31,[2 5]);
    plot_labels('Time [s]', 'Signal');
    y_avg=0;
    for j=1:length(inputArray11)
        if j==1
            continue
        end
%         if j==43
%             'lol'
%         end
        y21{j}=y21{j}/max(y21{j});
        p21(j)=plot_preliminaries(x21{j},smooth(y21{j},1000)+0.05*j,j,'nomarker'); %normalized
        set(p21(j),'DisplayName','2800 Hz'); 
        title('2800 Hz');
        xlim([0 20]);
     
    end
    savefig(h31,[dirname 'fig7a.fig']);
    
h32=start_fig(32,[2 5]);
    plot_labels('Time [s]', 'Signal');
    y_avg=0;
    for j=1:length(inputArray11)
        if j==1
            continue
        end
%         if j==43
%             'lol'
%         end
        y22{j}=y22{j}/max(y22{j});
        p22(j)=plot_preliminaries(x22{j},smooth(y22{j},1000)+0.05*j,j,'nomarker'); %normalized
        set(p22(j),'DisplayName','2800 Hz'); 
        title('2885 Hz');
        xlim([0 20]);
     
    end
    savefig(h32,[dirname 'fig7b.fig']);