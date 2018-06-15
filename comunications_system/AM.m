fs=1e6; %采样率
N=8192;
t=0:1/fs:(N-1)/fs; %时域采样点
f=(0:(N-1))*fs/N-fs/2; %频域采样点
fm=2000; %基带频率
fss=20000; %载波频率
m=cos(fm*2*pi*t); %基带信号

subplot(4,1,1);
plot(t,m);
title(' 基带信号','FontWeight','bold');
xlabel('t/s','FontSize',12);
axis([0.001,0.004,-1,1]);

m_w=fft(m); %基带频谱
subplot(4,1,2);
plot(f,fftshift(abs(m_w)));
title(' 基带频谱','FontWeight','bold');
xlabel('f/Hz','FontSize',12);
axis([-10000,10000,0,2000]);

s=cos(fss*2*pi*t); %载波信号
subplot(4,1,3);
plot(t,s);
title(' 载波信号','FontWeight','bold');
xlabel('t/s','FontSize',12);
axis([0.001,0.004,-1,1]);

s_w=fft(s); %载波频谱
subplot(4,1,4);
plot(f,fftshift(abs(s_w)));
title(' 载波频谱','FontWeight','bold');
xlabel('f/Hz','FontSize',12);
axis([-30000,30000,0,2000]);

A=2;
sm=s.*(m+A); %调制信号
figure;
subplot(2,1,1);
plot(t,sm);
title(' 调制信号','FontWeight','bold');
xlabel('t/s','FontSize',12);
axis([0.001,0.004,-5,5]);

sm_w=fft(sm); %调制信号频谱
subplot(2,1,2);
plot(f,fftshift(abs(sm_w)));
title(' 调制信号频谱','FontWeight','bold');
xlabel('f/Hz','FontSize',12);
axis([-40000,40000,0,2000]);

figure;
sm1=awgn(sm,5); %加噪声后的调制信号
subplot(4,1,1);
plot(t,sm1);
title(' 加噪声后的调制信号','FontWeight','bold');
xlabel('t/s','FontSize',12);
axis([0.001,0.004,-5,5]);



%滤波器设计？？？？
fsamp = 1e6; %采样频率为1MHz
fcuts = [16000 17500 22500 24000];
mags = [0 1 0];
devs = [0.05 0.01 0.05];
[n,Wn,beta,ftype] =kaiserord(fcuts,mags,devs,fsamp);
hh =fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale');
sm2=fftfilt(hh,sm1);

subplot(4,1,2); % 滤掉部分噪声后的调制信号
plot(t,sm2);
title(' 滤掉部分噪声后的调制信号','FontWeight','bold');
xlabel('t/s','FontSize',12);
axis([0.001,0.004,-5,5]);

sp=2*sm2.*s; % 与本地载波相乘后的信号
subplot(4,1,3);
plot(t,sp);
title(' 与本地载波相乘后的信号','FontWeight','bold');
xlabel('t/s','FontSize',12);
axis([0.001,0.004,-5,5]);

fsamp = 1e6; %采样频率为1MHz
fcuts = [3000 20000];
mags = [1 0];
devs = [0.01 0.05]; % 通带波动1%，阻带波动5%
[n,Wn,beta,ftype] =kaiserord(fcuts,mags,devs,fsamp);
hh1=fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale');
sd=fftfilt(hh1,sp); % 低通滤波后的信号
subplot(4,1,4);
plot(t,sd);
title(' 低通滤波后的信号','FontWeight','bold');
xlabel('t/s','FontSize',12);
axis([0.001,0.004,0,4]);


figure;
sm1_w=fft(sm1); % 加噪声后的调制信号频谱
subplot(4,1,1);
plot(f,fftshift(abs(sm1_w)));
title(' 加噪声后的调制信号频谱','FontWeight','bold');
xlabel('f/Hz','FontSize',12);
axis([-40000,40000,0,2000]);

sm2_w=fft(sm2); % 滤掉部分噪声后的调制信号频谱
subplot(4,1,2);
plot(f,fftshift(abs(sm2_w)));
title(' 滤掉部分噪声后的调制信号频谱','FontWeight','bold');
xlabel('f/Hz','FontSize',12);
axis([-40000,40000,0,2000]);

sp_w=fft(sp); % 与本地载波相乘后的信号频谱
subplot(4,1,3);
plot(f,fftshift(abs(sp_w)));
title(' 与本地载波相乘后的信号频谱','FontWeight','bold');
xlabel('f/Hz','FontSize',12);
axis([-45000,45000,0,1200]);

sd_f=fft(sd); % 低通滤波后的信号频谱
subplot(4,1,4);
plot(f,fftshift(abs(sd_f)));
title(' 低通滤波后的信号频谱','FontWeight','bold');
xlabel('f/Hz','FontSize',12);
axis([-40000,40000,0,3500]);
