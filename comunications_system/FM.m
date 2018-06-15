clear;
fs=1e6;                   %采样率
N=8192;
fm=2000;                  %基带频率
fss=20000;                %载波频率
t=0:1/fs:(N-1)/fs;        %时域采样点
f=(0:(N-1))*fs/N-fs/2;    %频域采样点
fmk=5;

m=cos(fm*2*pi*t);         	        %基带信号
subplot(4,1,1);
plot(t,m);
title('基带信号','FontWeight','bold');
xlabel('t/s');
axis([0.001,0.004,-1.5,1.5]);
m_pp=fft(m);              	        %基带频谱
subplot(4,1,2);
plot(f,fftshift(abs(m_pp)));
title('基带频谱','FontWeight','bold');
xlabel('f/Hz');
axis([-10000,10000,0,3500]);

w1 = 0 ;w2 = 0; 
for n = 1:length(t) 
w1 = m(n) + w2;    
w2 = m(n) + w1; 
f0(n) = w1/(2*fss);    
end 
f0 = f0*2*pi/max(abs(f0)*2*pi);        %积分后信号
subplot(4,1,3);
plot(t,f0);
title('积分后信号','FontWeight','bold');
xlabel('t/s');
axis([0.001,0.004,-1.5,1.5]);
f0_pp=fft(f0);                          %积分后信号频谱
subplot(4,1,4);
plot(f,fftshift(abs(f0_pp)));
title('积分后信号频谱','FontWeight','bold');
xlabel('f/Hz');
axis([-10000,10000,0,3500]);
fi=f0*fmk;
figure;

s=cos(fss*2*pi*t);        	       %载波信号
subplot(4,1,1);
plot(t,s);
title('载波信号','FontWeight','bold');
xlabel('t/s');
axis([0.001,0.004,-2,2]);
s_pp=fft(s);             	       %载波频谱
subplot(4,1,2);
plot(f,fftshift(abs(s_pp)));
title('载波频谱','FontWeight','bold');
xlabel('f/Hz');
axis([-30000,30000,0,3500]);

sfm=s.*cos(fi)-sin(fss*2*pi*t).*sin(fi);   %调制后信号
subplot(4,1,3);
plot(t,sfm);
title('调制后信号','FontWeight','bold');
xlabel('f/Hz');
axis([0.001,0.005,-2,2]);
sfm_pp=fft(sfm);          	     %调制后信号频谱
subplot(4,1,4);
plot(f,fftshift(abs(sfm_pp)));
title('调制后信号频谱','FontWeight','bold');
xlabel('f/Hz');
axis([-40000,40000,0,1500]);
figure;

sfm1=awgn(sfm,30);                   %加噪后调制信号
subplot(2,1,1);
plot(t,sfm1);
title('加噪后调制信号','FontWeight','bold');
xlabel('t/s');
axis([0.001,0.005,-2,2]);
sfm1_pp=fft(sfm1);                     %加噪后调制信号频谱
subplot(2,1,2);
plot(f,fftshift(abs(sfm1_pp)));
title('加噪后调制信号频谱','FontWeight','bold');
xlabel('f/Hz');
axis([-40000,40000,0,1500]);
figure;

fsamp = 1e6; 			     %采样频率为1MHz 
fcuts = [1000 10000 37000 38000];  
mags = [0 1 0]; 
devs = [0.05 0.01 0.05]; 
[n,Wn,beta,ftype] = kaiserord(fcuts,mags,devs,fsamp); 
hh = fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale'); 
sfm2=fftfilt(hh,sfm1);              %带通滤波后调制信号
subplot(2,1,1);                     
plot(t,sfm2);
title('带通滤波后调制信号','FontWeight','bold');
xlabel('t/s');
axis([0.001,0.004,-2,2]);
sfm2_pp=fft(sfm2);                   %带通滤波后调制信号频谱
subplot(2,1,2);
plot(f,fftshift(abs(sfm2_pp)));
title('带通滤波后调制信号频谱','FontWeight','bold');
xlabel('f/Hz');
axis([-40000,40000,0,1500]);
figure;

for i=1:length(t)-1                 %接收信号通过微分器处理
sd0(i)=(sfm2(i+1)-sfm2(i))/(1/fs); 
end 
sd= abs(hilbert(sd0));
subplot(2,1,1);
plot(sd);
title('恢复出的信号','FontWeight','bold');
xlabel('t/s');
axis([3733,8000,0,200000]);
fsamp = 1e6;                       %采样频率为1MHz 
fcuts = [0 1000 7000 8500];  
mags = [0 1 0]; 
devs = [0.05 0.01 0.05]; 
[n,Wn,beta,ftype] = kaiserord(fcuts,mags,devs,fsamp); 
hh = fir1(n,Wn,ftype,kaiser(n+1,beta),'noscale'); 
sd1=fftfilt(hh,sd);                %滤掉直流分量的调制信号
subplot(2,1,2);
plot(sd1);
title('滤掉直流分量的信号','FontWeight','bold');
xlabel('t/s');
axis([3733,8000,-8e4,8e4]);
