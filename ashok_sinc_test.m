clear;clc;close all;
%% 0) Setup the parameters
dt=1e-9; %sampling time 1ns
tacq=32e-6; %acqusition time
t=0:dt:tacq-dt;
wL=2*pi*20e6; %Larmor freq
wac=2*pi*2000; %AC frequency
%wac=2*pi/4/73e-6; %resonance for pi/2
g=1e-2; %strength of AC field
phi=0;

%% 1) Preliminary checks, not important for later
S=cos(wL*t);
S2=cos(wL*t + g*cos(wac*t+phi));

% Plot to check
start_fig(1,[3 2]);
plot_preliminaries(t,S,1,'nomarker');
plot_preliminaries(t,S2,2,'nomarker');
set(gca,'ylim',[-1.1 1.1]);
plot_labels('Time [s]','Signal');

% Take fft
Y=fft(S)/length(S);Y=abs(fftshift(Y));
Y2=fft(S2)/length(S2);Y2=abs(fftshift(Y2));

n=size(Y,2);
fs=1/max(t);
F=(-n/2:n/2-1)*(fs);
start_fig(2,[3 2]);
plot_preliminaries(F/1e6,log10(Y),1);
plot_preliminaries(F/1e6,log10(Y2),2);
p1=plot_yline(20,5);set(p1,'linestyle','--');
set(gca,'xlim',[19.8 20.2]);

plot_labels('Frequency [MHz]','Log Signal');


max(Y)-max(Y2)

%% 2) Make phi change in a loop

phi_vec=0:pi/32:4*pi;
for j=1:size(phi_vec,2)
    phi=phi_vec(j);
    S=cos(wL*t);
    S2=cos(wL*t + g*(sin(wac*t+phi)-sin(phi)));
    % Take fft
    Y=fft(S)/length(S);Y=abs(fftshift(Y));
    Y2=fft(S2)/length(S2);Y2=abs(fftshift(Y2));
    
    signal(j)=max(Y2);
    %signal(j)=max(Y)-max(Y2);
end

start_fig(3,[3 2]);
plot_preliminaries(phi_vec,signal,2);
plot_labels('Phase','Signal');


%% 3) Main simulation of AC magnetometry, keeping track of global phase
tau=73e-6; %center-to-center pulse distance
T=20e-3; %total time
phi_0=0;
g=1e-2; %strength of AC field in frequency units over fac
cycles=floor(T/tau);
for j=1:cycles
    phi_L=wL*(j-1)*tau; %keep track of the phase of Larmor frequency
    phi_AC=wac*(j-1)*tau+phi_0; %keep track of the phase of the AC field
    %S=cos(wL*t+phi_L);
    S2=cos(wL*t+phi_L + g*(sin(wac*t+phi_AC)-sin(phi_AC)));
    % Take fft
    %Y=fft(S)/length(S);Y=abs(fftshift(Y));
    Y2=fft(S2)/length(S2);Y2=abs(fftshift(Y2));
    
    %signal(j)=max(Y)-max(Y2);
    signal(j)=max(Y2);
    time_vector(j)=j*tau;
end

start_fig(4,[3 2]);
plot_preliminaries(time_vector*1e3,signal,1);
plot_labels('Time [ms]','Signal');

% Take fft
P=fft(signal-mean(signal))/length(signal);P=abs(fftshift(P));
np=size(P,2);
fs_p=1/(cycles*tau);
freq=(-np/2:np/2-1)*(fs_p);
start_fig(5,[3 2]);
p11=plot_preliminaries(freq,P,1);
p1=plot_yline(wac/(2*pi),5);set(p1,'linestyle','--');
p1=plot_yline(2*wac/(2*pi),5);set(p1,'linestyle','--');
set(p11,'DisplayName','Chunked');
plot_labels('Frequency [Hz]','Signal');

%% continuous part of 3
tcon=dt:dt:cycles*tau;
S3=cos(wL*tcon+phi_0 + g*(sin(wac*tcon+phi_AC)-sin(phi_AC)));
P3=abs(fftshift(fft(S3)/length(S3)));
F3=fftfreq(tcon,true);
[~, IX]=min(abs(F3-wL/2/pi));
P3(IX)=0;
p2=plot_preliminaries(F3-wL/2/pi,P3,2);
set(p2,'DisplayName','Continuous');
xlim([-5000 5000]);
%ylim([-5 30]);

legend([p11,p2])

%% 4) Check scaling with AC amplitude g

tau=73e-6; %center-to-center pulse distance
T=20e-3; %total time
phi_0=pi/2; %AC field phase
g=1e-3; %strength of AC field in frequency units
cycles=floor(T/tau);
tcon=dt:dt:cycles*tau; %continuous time array
gvec=0:1e-2:20e-2;

for k=1:size(gvec,2)
    g=gvec(k);
    for j=1:cycles
        phi_L=wL*(j-1)*tau; %keep track of the phase of Larmor frequency
        phi_AC=wac*(j-1)*tau+phi_0; %keep track of the phase of the AC field
        S=cos(wL*t+phi_L);
        S2=cos(wL*t+phi_L + g*(sin(wac*t+phi_AC)-sin(phi_AC)));
        % Take fft
        Y=fft(S)/length(S);Y=abs(fftshift(Y));
        Y2=fft(S2)/length(S2);Y2=abs(fftshift(Y2));
         
        signal(j)=max(Y2);
        %signal(j)=max(Y)-max(Y2);
         time_vector(j)=j*tau;
    end
    % Take fft
    P=fft(signal-mean(signal));P=abs(fftshift(P))/length(P);
    np=size(P,2);
    fs_p=1/T;
    freq=(-np/2:np/2-1)*(fs_p);
    
    
    start_fig(5,[3 2]);
    p11=plot_preliminaries(freq,P,1);
    p1=plot_yline(wac/(2*pi),5);set(p1,'linestyle','--');
    p1=plot_yline(2*wac/(2*pi),5);set(p1,'linestyle','--');
    set(p11,'DisplayName','Chunked');
    plot_labels('Frequency [Hz]','Signal');
    
    
    [v ,b]=min(abs(freq-wac/(2*pi)));
    [v2 ,b2]=min(abs(freq-2*wac/(2*pi)));
     
    primary(k)=P(b);secondary(k)=P(b2);
    
    % continuous part
    
    S3=cos(wL*tcon+phi_0 + g*(sin(wac*tcon+phi_0)-sin(phi_0)));
    P3=abs(fftshift(fft(S3)/length(S3)));
    F3=fftfreq(tcon,true);
    [~, IX]=min(abs(F3-wL/2/pi));
    P3(IX)=0;
    p22=plot_preliminaries(F3-wL/2/pi,P3,2);
    set(p22,'DisplayName','Continuous');
    xlim([-5000 5000]);
    %ylim([-5 30]);
    
    legend([p11,p22]);
    
    [~ ,bcon]=min(abs(F3-wL/2/pi-wac/(2*pi)));
    [~ ,bcon2]=min(abs(F3-wL/2/pi-2*wac/(2*pi)));
    primarycon(k)=P3(bcon); secondarycon(k)=P3(bcon2);
end
start_fig(6,[3 2]);
plot_preliminaries(gvec,primary,1);
plot_preliminaries(gvec,secondary,2);
plot_labels('AC field strength g [Hz]','Signal');
title('Chunked Case');

start_fig(7,[3 2]);
plot_preliminaries(gvec,primary./gvec,1);
plot_preliminaries(gvec,secondary./gvec,2);
plot_labels('AC field strength g [Hz]','Signal/g');
title('Chunked Case');

start_fig(8,[3 2]);
plot_preliminaries(gvec,primarycon,1);
plot_preliminaries(gvec,secondarycon,2);
plot_labels('AC field strength g','Signal');
title('Continuous Case');
legend('First Harmonic','Second Harmonic','Location','best');

start_fig(9,[3 2]);
plot_preliminaries(gvec,primarycon./gvec,1);
plot_preliminaries(gvec,secondarycon./gvec,2);
plot_labels('AC field strength g','Signal/g');
title('Continuous Case');
legend('First Harmonic','Second Harmonic','Location','best');

%% 5) Same as (3) but simulating also phase of the signal
tau=73e-6; %center-to-center pulse distance
T=20e-3; %total time
g=1e-2; %strength of AC field in frequency units
phi_0=pi/2;
cycles=floor(T/tau);
for j=1:cycles
    phi_L=wL*(j-1)*tau; %keep track of the phase of Larmor frequency
    phi_AC=wac*(j-1)*tau+phi_0; %keep track of the phase of the AC field
    S=cos(wL*t+phi_L);
    S2=cos(wL*t+phi_L + g*(sin(wac*t+phi_AC)-sin(phi_AC)));
    % Take fft
    Y=fft(S);Yamp=abs(fftshift(Y));Yangle=angle(fftshift(Y));
    Y2=fft(S2);Y2amp=abs(fftshift(Y2));Y2angle=angle(fftshift(Y2));
    n=size(Y,2);
    fs=1/max(t);
    F=(-n/2:n/2-1)*(fs);
    
    [v b]=min(abs(F-20e6));
    signal(j)=Yamp(b)-Y2amp(b);
    phase(j)=Yangle(b)-Y2angle(b);
    time_vector(j)=j*tau;
end

start_fig(9,[3 2]);
%plot_preliminaries(time_vector*1e3,signal/g,2);
plot_preliminaries(time_vector*1e3,phase/g,3);
plot_labels('Time [ms]','Signal/g');

% Take fft of signal and phase
P=fft(signal);P=abs(fftshift(P));
P2=fft(phase);P2=abs(fftshift(P2));
np=size(P,2);
fs_p=1/(cycles*tau);
freq=(-np/2:np/2-1)*(fs_p);

start_fig(10,[3 2]);
%plot_preliminaries(freq,P,1);
plot_preliminaries(freq,P2,3);
p1=plot_yline(wac/(2*pi),5);set(p1,'linestyle','--');
p1=plot_yline(2*wac/(2*pi),5);set(p1,'linestyle','--');
plot_labels('Frequency [Hz]','Signal');


%% 6) Main simulation of AC magnetometry with resonant pi case, keeping track of global phase
tau=73e-6; %center-to-center pulse distance
wac=1/2/tau*2*pi;
T=20e-3; %total time
g=1e-2; %strength of AC field in frequency units over fac
phi_0=pi/2;
cycles=floor(T/tau);
for j=1:cycles
     phi_L=wL*(j-1)*tau+(-1)^j*2*(j-1)*g*sin(phi_0); %keep track of the phase of Larmor frequency
    phi_AC=wac*(j-1)*tau; %keep track of the phase of the AC field
    %S=cos(wL*t+phi_L);
    S2=cos(wL*t+phi_L + g*(sin(wac*t+phi_AC)-sin(phi_AC)));
    % Take fft
    %Y=fft(S)/length(S);Y=abs(fftshift(Y));
    Y2=fft(S2)/length(S2);Y2=abs(fftshift(Y2));
    
    signal(j)=max(Y2);
    time_vector(j)=j*tau;
end

start_fig(4,[3 2]);
plot_preliminaries(time_vector*1e3,signal/g,2);
plot_labels('Time [ms]','Signal/g');

% Take fft
P=fft(signal-mean(signal))/length(signal);P=abs(fftshift(P));
np=size(P,2);
fs_p=1/(cycles*tau);
freq=(-np/2:np/2-1)*(fs_p);
start_fig(5,[3 2]);
plot_preliminaries(freq,P,1);
p1=plot_yline(wac/(2*pi),5);set(p1,'linestyle','--');
p1=plot_yline(2*wac/(2*pi),5);set(p1,'linestyle','--');
plot_labels('Frequency [Hz]','Signal');
