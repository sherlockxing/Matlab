clear;
Fs=1e4;                        %采样频率
len=20;                        %码元长度
in=randint(1,len,4);          %产生初始码元序列
sig=[];
out=[];
for    t=1:2000              %产生基带信号
n=fix(t/100);
if n==0
in_a(t)=0;
else
in_a(t)=in(n);
end
end
subplot(2,1,1);              %基带信号
plot(in_a,'LineWidth',3);
title('基带信号','FontWeight','bold','FontSize',20);
xlabel('t/s','FontSize',18);
axis([100,2100,-0.5,3.5]);
set(gca,'XTick',0:100:2000);

grid on;
cxn=xcorr(in_a,'unbiased'); %%计算序列的自相关函数
nfft=1024;
CXk=fft(cxn,nfft);
Pxx=abs(CXk);
index=0:round(nfft/2-1);
k=index*Fs/nfft;
subplot(2,1,2);
plot_Pxx=10*log10(Pxx(index+1));
plot(k,plot_Pxx,'LineWidth',2);
title('基带信号功率谱','FontWeight','bold','FontSize',20);
axis([0,5000,-10,40]);
xlabel('Hz','FontSize',18,'FontSize',18);
for    i=1:len        %产生单极性归零码信号
if  in(i)==0
ins=[0,0];
elseif in(i)==1
ins=[1,0];
elseif in(i)==2
ins=[2,0];
else
ins=[3,0];
end
sig=[sig,ins];
end
for    t=1:4000
n=fix(t/100);
if n==0
s(t)=0;
else
s(t)=sig(n);
end
end
figure;
subplot(2,1,1);            %单极性归零码
plot(s,'LineWidth',3);
title('单极性归零码','FontWeight','bold','FontSize',20);
xlabel('t/s','FontSize',18);
axis([100,4100,-0.5,3.5]);
set(gca,'XTick',0:200:4100);
grid on;
cxn=xcorr(s,'unbiased');
nfft=1024;
CXk=fft(cxn,nfft);
Pxx=abs(CXk);
index=0:round(nfft/2-1);
k=index*Fs/nfft;
subplot(2,1,2);
plot_Pxx=10*log10(Pxx(index+1));
plot(k,plot_Pxx,'LineWidth',2);
title('单极性归零码功率谱','FontWeight','bold','FontSize',20);
axis([0,5000,-10,40]);
xlabel('Hz','FontSize',18);
s1=awgn(s,20);              %添加噪声
figure;
subplot(2,1,1);
plot(s1);
title('添加噪声后的信号','FontWeight','bold','FontSize',20);
xlabel('t/s','FontSize',18);
axis([100,4100,-0.5,3.5]);
set(gca,'XTick',0:500:4100);
cxn=xcorr(s1,'unbiased');
nfft=1024;
CXk=fft(cxn,nfft);
Pxx=abs(CXk);
index=0:round(nfft/2-1);
k=index*Fs/nfft;
subplot(2,1,2);
plot_Pxx=10*log10(Pxx(index+1));
plot(k,plot_Pxx,'LineWidth',2);
title('添加噪声后的信号功率谱','FontWeight','bold','FontSize',20);
axis([0,5000,-10,40]);
xlabel('Hz','FontSize',18);    %滤波器设计
fp=500;                            %通带截止
fs= 550;                        %阻带截止
ws=fs*2/Fs;
wp=fp*2/Fs;
[N, Wp] = ellipord(wp,ws,1,40);
[b,a]=ellip(N,1,40,Wp);
sf0=filter(b,a,s1) ;    %滤掉部分噪声后的信号
figure;
subplot(2,1,1);
plot(sf0);
title('滤掉部分噪声后的信号','FontWeight','bold','FontSize',20);
xlabel('t/s','FontSize',18);
set(gca,'XTick',0:500:4100);
cxn=xcorr(sf0,'unbiased'); nfft=1024;
CXk=fft(cxn,nfft);
Pxx=abs(CXk);
index=0:round(nfft/2-1);
k=index*Fs/nfft;
subplot(2,1,2);
plot_Pxx=10*log10(Pxx(index+1));
plot(k,plot_Pxx,'LineWidth',2);
title('滤掉部分噪声后的信号功率谱','FontWeight','bold','FontSize',20);
axis([0,5000,-10,40]);
xlabel('Hz','FontSize',18,'FontSize',18);
for m=1:20                  %抽样判决
p=round(sf0(200*m-20));
out=[out,p];
end
figure;
subplot(2,1,1);        %原始码元序列
stairs(in,'LineWidth',3);
title('原始码元序列','FontWeight','bold','FontSize',20);
axis([1,21,-0.5,3.5]);
set(gca,'XTick',0:1:20);
grid on;
subplot(2,1,2);  %抽样判决后恢复的信号序列
stairs(out,'LineWidth',3);
title('抽样判决后恢复的信号序列','FontWeight','bold','FontSize',20);
axis([1,21,-0.5,3.5]);
set(gca,'XTick',0:1:20);
grid on;
