%modified by jiangxin,2017-07-29
%E-mail: jiangxin_cdut@163.com
clear;clc;
tic;

%************************�׿��Ӳ�****************************
dt=0.001;t=0:dt:0.04;f0=40;
R=(1-2*(pi*f0*t).^2).*exp(-(pi*f0*t).^2);
%plot(t,R);
%************************************************************

%***********************ģ�Ͳ�������*************************
dx=4;dy=4;%�������
x=0:dx:400;y=0:dy:400;%��������
% x=-200:dx:200;y=-200:dy:200;%��������
x0=200;y0=200;%�ڵ㼤���ص�
dx0=10;%��ձ߽�
z0=(12+10*sin(5*pi*x0/100/4)*sin(5*pi*y0/100/4)+3);
y0=floor(z0/2)*2;
%y0=200;
t0=0.5;%����ʱ��
l=length(t);m=length(x);n=length(y);
% Setting Velocity & Density
xt = 0:dx/2:400+dx/2;yt = 0:dy/2:400+dy/2;
% xt = -200:dx/2:200+dx/2;yt = -200:dy/2:200+dy/2;
mt = length(xt);nt = length(yt);
vp=zeros(mt,nt);vs=zeros(mt,nt);pp=zeros(mt,nt);
y1=200;m1=find(yt==y1);%����λ��
y2=320;m2=find(yt==y2);%����λ��
% ************ģ��2**********
for i=2:length(x)
    mm0(i-1)=(12+round(10*sin(5*pi*x(i)/4/100)*sin(5*pi*50/4/100)+0.5));
    mm0(i-1)=floor(mm0(i-1)/2)*2;
    m3(i-1)=find(yt==mm0(i-1));%����λ��
end
%plot(x(2:101),m0);
y1=200;m1=find(yt==y1);%����λ��
y2=320;m2=find(yt==y2);%����λ��
for i=1:mt
    for j=1:nt
        if j<=m1
            vp(i,j)=1200;vs(i,j)=vp(i,j)/sqrt(3);pp(i,j)=1200;
        elseif j>m1&&j<=m2
            vp(i,j)=1700;vs(i,j)=vp(i,j)/sqrt(3);pp(i,j)=1700;
        else
            vp(i,j)=2.3e3;vs(i,j)=vp(i,j)/sqrt(3);pp(i,j)=2.3e3;
        end
    end
end
%for i=1:mt
%    for j=1:nt
%        if floor(i/4)==100
%            if j<=m1&&j>=m3(100)
%                vp(i,j)=1200;vs(i,j)=vp(i,j)/sqrt(3);pp(i,j)=1200;
%            elseif j>m1&&j<=m2
%                vp(i,j)=1700;vs(i,j)=vp(i,j)/sqrt(3);pp(i,j)=1700;
%            elseif j>m2
%                vp(i,j)=2.3e3;vs(i,j)=vp(i,j)/sqrt(3);pp(i,j)=2.3e3;
%            end
%        else
%            if j<=m1&&j>=m3(floor(i/4)+1)
%                vp(i,j)=1200;vs(i,j)=vp(i,j)/sqrt(3);pp(i,j)=1200;
%            elseif j>m1&&j<=m2
%                vp(i,j)=1700;vs(i,j)=vp(i,j)/sqrt(3);pp(i,j)=1700;
%            elseif j>m2
%                vp(i,j)=2.3e3;vs(i,j)=vp(i,j)/sqrt(3);pp(i,j)=2.3e3;
%            end
%        end
%    end
%end
N = pp.*vs.^2;%���Բ���
L = pp.*(vp.^2-2*vs.^2);
r=2;C = DCoef(r,'s');% ��ȡ�ռ�߽׵�ʦϵ��
%************************************************************

%************************�ȶ�������**************************
maxvp=max(max(vp));maxvs=max(max(vs));
if dt*sqrt(maxvp^2/dx^2+maxvs^2/dy^2)>1
%if dt*maxvp*sqrt(1/dx^2+1/dy^2)>(1/(abs(C(1))+abs(C(2))))
    error('�������������⣬��ָ��ȶ�');
end
%************************************************************

%**********************��ֵ������FCT��Ƶɢ*******************
P=zeros(m,n);Q=P;S=P;U=P;W=P;
P1=P;Q1=Q;S1=S;U1=U;W1=W;
m0=find(x==x0);n0=find(y==y0);l0=round(t0/dt)+1;
Param.pp=pp;Param.N=N;Param.L=L;% ���ڱ߽纯����������
Param.dt=dt;Param.dx=dx;Param.dy=dy;Param.r=r;% ���ڱ߽纯����������
dataX=zeros(m,n,l0);dataY=zeros(m,n,l0);VSP=zeros(m,l0);
for k=1:l0
    imagesc(W);colormap('gray');drawnow;
    if k<=length(t)
        P1(m0,n0)=R(k);Q1(m0,n0)=R(k);%S1(m0,n0)=R(k);
    end
    for i=r+1:m-r
        for j=r+1:n-r
%            if j<=dx0
%                Q1(i,j-1)=-Q1(i,j+2);Q1(i,j)=-Q1(i,j+1);
%                S1(i,j-1)=-S1(i,j+1);
%                S1(i,j)=0;
%            end
            U(i,j)=U1(i,j)+C(1)*dt/pp(2*i-1,2*j-1)*(1/dx*(P1(i,j)-P1(i-1,j))+1/dy*(S1(i,j)-S1(i,j-1)))+...
                C(2)*dt/pp(2*i-1,2*j-1)*(1/dx*(P1(i+1,j)-P1(i-2,j))+1/dy*(S1(i,j+1)-S1(i,j-2)));
            W(i,j)=W1(i,j)+C(1)*dt/pp(2*i,2*j)*(1/dx*(S1(i+1,j)-S1(i,j))+1/dy*(Q1(i,j+1)-Q1(i,j)))+...
                C(2)*dt/pp(2*i,2*j)*(1/dx*(S1(i+2,j)-S1(i-1,j))+1/dy*(Q1(i,j+2)-Q1(i,j-1)));
        end
    end
%     U = VxBoundary(U,U1,P1,S1,Param);%Vx���ձ߽�����
%     W = VyBoundary(W,W1,Q1,S1,Param);%Vy���ձ߽�����
%     U=FCTforEW(U1,U);W=FCTforEW(W1,W);
    for i=r+1:m-r
        for j=r+1:n-r
            P(i,j)=P1(i,j)+C(1)*dt*((L(2*i,2*j-1)+2*N(2*i,2*j-1))/dx*(U(i+1,j)-U(i,j))+L(2*i,2*j-1)/dy*(W(i,j)-W(i,j-1)))+...
                dt*C(2)*((L(2*i,2*j-1)+2*N(2*i,2*j-1))/dx*(U(i+2,j)-U(i-1,j))+L(2*i,2*j-1)/dy*(W(i,j+1)-W(i,j-2)));
            Q(i,j)=Q1(i,j)+C(1)*dt*(L(2*i,2*j-1)/dx*(U(i+1,j)-U(i,j))+(L(2*i,2*j-1)+2*N(2*i,2*j-1))/dy*(W(i,j)-W(i,j-1)))+...
                dt*C(2)*(L(2*i,2*j-1)/dx*(U(i+2,j)-U(i-1,j))+(L(2*i,2*j-1)+2*N(2*i,2*j-1))/dy*(W(i,j+1)-W(i,j-2)));
            S(i,j)=S1(i,j)+C(1)*dt*N(2*i-1,2*j)*(1/dx*(W(i,j)-W(i-1,j))+1/dy*(U(i,j+1)-U(i,j)))+...
                dt*C(2)*N(2*i-1,2*j)*(1/dx*(W(i+1,j)-W(i-2,j))+1/dy*(U(i,j+2)-U(i,j-1)));
%            if j<=dx0
%                S(i,j)=0;
%            end
         end
    end
%      P = TxxBoundary(U,W,P,P1,Q1,Param);%Txx���ձ߽�����
%      Q = TyyBoundary(U,W,Q,P1,Q1,Param);%Tyy���ձ߽�����
%      S = TxyBoundary(U,W,S,S1,Param);%Txy���ձ߽�����

     U1=U;W1=W;P1=P;Q1=Q;S1=S;
     VSP(:,k)=W(:,4);
     dataX(:,:,k)=U;dataY(:,:,k)=W;
end
  %**************************���ݳ�ͼ**************************
 figure(1);
 subplot(1,2,1);
[Y,X]=meshgrid(y,x);
surf(X',Y',U');shading interp;view(0,90);colormap('gray');
% hold on;
% plot([-400,0],[50,50],'r',[0 0],[50 400],'r');
% hold off;
axis square;axis tight;axis ij;
title('X����');xlabel('X');ylabel('Y');
subplot(1,2,2);
[Y,X]=meshgrid(y+dy/2,x+dx/2);
surf(X',Y',W');shading interp;view(0,90);colormap('gray');
% hold on;
% plot([-400,0],[50,50],'r',[0 0],[50 400],'r');
% hold off;
axis square;axis tight;axis ij;
title('Y����');xlabel('X');ylabel('Y');
set(gcf,'color','w');
%************************************************************
%**************************��������**************************
 [Y,X]=meshgrid(y,x);
 figure(1);
 %M1=moviein(10); 
 for k=1:l0
     h=surf(X',Y',dataX(:,:,k)');
     shading interp;view(0,90);colormap gray;
     axis square;axis tight;axis ij;material dull;
     hold on;
     plot([0,400],[200,200],'r');
     plot([0,400],[320,320],'r');
     hold off;
     title('X����');xlabel('X');ylabel('Y');
     M1(:,k)=getframe;
 end
 %movie(M1,2)
 [Y,X]=meshgrid(y+dy/2,x+dx/2);
 for k=1:l0
     h=surf(X',Y',dataY(:,:,k)');
     shading interp;view(0,90);colormap gray;
     axis square;axis tight;axis ij;material dull;
     hold on;
     plot([0,400],[200,200],'r');
     plot([0,400],[320,320],'r');
     hold off;
     title('Y����');xlabel('X');ylabel('Y');
     M2(:,k)=getframe;
 end
 movie2avi(M1,'Model2XVector','fps',10);
 movie2avi(M2,'Model2YVector','fps',10);
%************************************************************
%************************************************************
%************************************************************
pickup=find(y==y0);%�첨��λ��
 for k=1:l0
     for i=2:m
        dataplotX(k,i-1)=dataX(i,pickup,k);
        dataplotY(k,i-1)=dataY(i,pickup,k);
     end
 end
 %for i=1:m-1
 %   datamaxX(i)=max(abs(dataplotX(:,i)));
 %end
%��һ������
%datamaxX=max(max(abs(dataplotX)));
%datamaxY=max(max(abs(dataplotY)));
%dataplotX=dataplotX./datamaxX.*2+2;
%dataplotY=dataplotY./datamaxY.*2+2;
dataplotTT=(0:(l0-1))*dt;
figure(1);
subplot(1,2,1);
for i=1:m-1
    datamaxX=max(abs(dataplotX(:,i)));
    if datamaxX==0
        dataplotX(:,i)=dataplotX(:,i)+0.5;
    else
        dataplotX(:,i)=dataplotX(:,i)./datamaxX./2+0.5;
    end
    plot(dataplotX(:,i)'+i-1,dataplotTT,'b');
    hold on;
    for k=1:l0
        if dataplotX(k,i)>0.5
            plot([0.5+i-1,dataplotX(k,i)+i-1],[dataplotTT(k),dataplotTT(k)]);
            hold on;
        end
    end
    hold on;
end
axis([0,101,0,t0]);axis ij;
title('x����');ylabel('ʱ��');xlabel('����');
subplot(1,2,2);
for i=1:m-1
    datamaxY=max(abs(dataplotY(:,i)));
    if datamaxY==0
        dataplotY(:,i)=dataplotY(:,i)+0.5;
    else
        dataplotY(:,i)=dataplotY(:,i)./datamaxY./2+0.5;
    end
    plot(dataplotY(:,i)'+i-1,dataplotTT,'b');
    for k=1:l0
        if dataplotX(k,i)>0.5
            plot([0.5+i-1,dataplotY(k,i)+i-1],[dataplotTT(k),dataplotTT(k)]);
            hold on;
        end
    end
    hold on;
end
axis([0,101,0,t0]);axis ij;
title('Y����');ylabel('ʱ��');xlabel('����');
