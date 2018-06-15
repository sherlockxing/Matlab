clear;clc;
%信息源-m序列码
N=1e4; msg=randint(1,N);
%--------------------------------------------------------------------------

%信道编码-卷积编码
trellis=poly2trellis(3,[7 5]);
code=convenc(msg,trellis);
%poly2strellis将卷积码多项式转换成MATLAB的trellis网格表达式的函数
%内部参数前者是卷积码的约束长度N，后者是根据输入输出连线情况的一个m*n矩阵，m为输入信号的个数，n为输出信号的个数
%trellis结构是将多项式用一种matlab能够理解的一种数据结构表示，节典型的用空间换取时间的方式。换了一种数据结构。
%convenc，使用卷积编码器，对二进制信息msg进行编码，trellis是编码器的trellis结构（网格结构）
%--------------------------------------------------------------------------

%QPSK调制
pcode=reshape(code,[2,length(code)/2]);
%B = reshape(A,size)是指返回一个和A元素相同的n维数组。
bitA=2*pcode(1,:)-1;
bitB=2*pcode(2,:)-1;
ABcomplex=bitA+1i*bitB;
%双极性序列a，b。
scatterplot(ABcomplex,1,0,'r*');%画散点图
axis([-2,2,-2,2]);
ax=gca;
ax.XTick=(-2:0.5:2);
ax.YTick=(-2:0.5:2);
title('a点信号星座图','Color','r');
%QPSK散点图
%axis([-2,2,-2,2]);
xlabel('实部-I','Color','r');
ylabel('虚部-Q','Color','r');
grid on;
%--------------------------------------------------------------------------

%成型滤波
nSamples=10;
%IIR数字滤波器截断深度8,滚降系数0.5,根升
ht=rcosdesign(0.5,8,nSamples,'sqrt');
%ht=rcosfir(0.5,4,10,0.4,'sqrt');
%归一化系统的冲激响应
htmax=max(ht(:));
ht=(ht)/max(ht(:));
len=length(ht);
vxright=(len-1)/2/nSamples;
vxleft=-vxright;
px=(vxleft:1/nSamples:vxright);

%插值填充0向量,提高频域分辨率
vec=[1;zeros(nSamples-1,1)];
%摘取bitA的前几个字符做示范
ABfill=vec*ABcomplex;
ABfill=ABfill(:);
ABshaped=conv(ht,ABfill);
len=length(ht);
row=size(ABcomplex,2);
px=((len-1)/2:1:(row-1)*nSamples+(len-1)/2)+1;
vx=(px-px(1))/nSamples;
%--------------------------------------------------------------------------

%b点眼图
eyediagram( ABshaped(px),2*nSamples,2,0,'b-');
subplot(2,1,1);
title('I路眼图');
xlabel('时间(s)');
ylabel('幅度');
grid on;
subplot(2,1,2);
title('Q路眼图');
xlabel('时间(s)');
ylabel('幅度');
grid on;  
%--------------------------------------------------------------------------

%添加噪声
ABnoise=awgn(ABshaped,10);
%--------------------------------------------------------------------------

%模拟调制信号在信道中的传输情况接收信号并通过接收滤波器解调并抽样
ABdemodTest=conv(ABshaped,ht)*(htmax)*(htmax);
ABdemod=conv(ABnoise,ht)*(htmax)*(htmax);
px=len:1:(row-1)*nSamples+len;
vx=(px-px(1))/nSamples;

%加噪后星座图
ABtemp=ABdemod(px);
ABdemodVal=ABtemp((1:nSamples:(row-1)*nSamples+1));
scatterplot(ABdemodVal,1,0,'b*');
title('加噪信号的QPSK信号的散点图','Color','b');
axis([-2,2,-2,2]);
xlabel('实部I','Color','b');
ylabel('虚部Q','Color','b');
grid on;
%--------------------------------------------------------------------------

%d点眼图
eyediagram(ABdemod(px),2*nSamples,2,0,'b-');
subplot(2,1,1);
title('I路眼图');
xlabel('时间(s)');
ylabel('幅度');
grid on;
subplot(2,1,2);
title('Q路眼图');
xlabel('时间(s)');
ylabel('幅度');
grid on;
%--------------------------------------------------------------------------

%抽样判决及解码
ABtemp=[real(ABdemodVal),imag(ABdemodVal)]';
PbitGet=(ABtemp>0);
bitget=PbitGet(:);
bitEnd=vitdec(bitget,trellis,5*2,'term','hard');
Z=double(bitEnd');
%--------------------------------------------------------------------------

%计算误码率
npcode=reshape(msg,[2,length(msg)/2]);
ndata=npcode'; 
nbitA=2*npcode(1,:)-1;
nbitB=2*npcode(2,:)-1;
nABcomplex=nbitA+1i*nbitB;
nN=N/2;
SNR_DB=[0:1:12]; 
a=1; 
Tb=1; 
Eb=a*a*Tb; 
P_signal=Eb/Tb; 
NO=Eb./(10.^(SNR_DB/10)); 
P_noise=P_signal*NO; 
sigma=sqrt(P_noise); 
for Eb_NO_id=1:length(sigma) 
    noise1=sigma(Eb_NO_id)*randn(1,nN); 
    noise2=sigma(Eb_NO_id)*randn(1,nN); 
    receive=nABcomplex+noise1+noise2*j; 
    resum=0; 
    total=0; 
    m1=find(angle(receive)<=pi/2&angle(receive)>0); 
    remessage(1,m1)=1+j; 
    redata(m1,1)=1; 
    redata(m1,2)=1; 
    m2= find( angle(receive)>pi/2&angle(receive)<=pi); 
    remessage(1,m2)=-1+j; 
    redata(m2,1)=0; 
    redata(m2,2)=1; 
    m3=find( angle(receive)>-pi&angle(receive)<=-pi/2); 
    remessage(1,m3)=-1-j; 
    redata(m3,1)=0; 
    redata(m3,2)=0; 
    m4=find( angle(receive)>-pi/2&angle(receive)<=0); 
    remessage(1,m4)=1-j; 
    redata(m4,1)=1; 
    redata(m4,2)=0; 
    [resum,ratio1]=symerr(ndata,redata); 
    pbit(Eb_NO_id)=resum/(nN*2); 
    [total,ratio2]=symerr(nABcomplex,remessage); 
    pe(Eb_NO_id)=total/nN; 
end 
Pe=1-(1-1/2*erfc(sqrt(10.^(SNR_DB/10)/2))).^2; 
Pbit=1/2*erfc(sqrt(10.^(SNR_DB/10)/2)); 
figure(11) 
data=pcode'; 
SNR_DB=[0:1:12]; 
a=1; 
Tb=1; 
Eb=a*a*Tb; 
P_signal=Eb/Tb; 
NO=Eb./(10.^(SNR_DB/10)); 
P_noise=P_signal*NO; 
sigma=sqrt(P_noise); 
for Eb_NO_id=1:length(sigma) 
    noise1=sigma(Eb_NO_id)*randn(1,N); 
    noise2=sigma(Eb_NO_id)*randn(1,N); 
    receive=ABcomplex+noise1+noise2*j; 
    resum=0; 
    total=0; 
    m1=find(angle(receive)<=pi/2&angle(receive)>0); 
    remessage(1,m1)=1+j; 
    m2= find( angle(receive)>pi/2&angle(receive)<=pi); 
    remessage(1,m2)=-1+j; 
    m3=find( angle(receive)>-pi&angle(receive)<=-pi/2); 
    remessage(1,m3)=-1-j; 
    m4=find( angle(receive)>-pi/2&angle(receive)<=0); 
    remessage(1,m4)=1-j;
    h=1;
    for i=1:length(remessage)
        total(h)=(real(remessage(i))+1)/2;
        total(h+1)=(imag(remessage(i))+1)/2;
        h=h+2;
    end
    hama = vitdec(total,trellis,1,'cont','hard'); 
    pe2(Eb_NO_id)=symerr(hama(2:N),msg(1:N-1))/(N-1); 
end
semilogy(SNR_DB,pe,':s',SNR_DB,Pe,'-*',SNR_DB,pe2,'-.') 
legend('QPSK仿真误码率','QPSK理论误码率','卷积后误码率',1) 
xlabel('信噪比/dB') 
ylabel('概率P') 
