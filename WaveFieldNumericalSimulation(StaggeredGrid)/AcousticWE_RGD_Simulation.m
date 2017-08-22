% �˳�������ģ������������δ�ӱ߽�������δ����ֵƵɢ��
clear;clc;figure(2);
% �׿��Ӳ�
tic;
dt=0.001;t=(0:40)*dt;f=30;
R=(1-2*(pi*f*t).^2).*exp(-(pi*f*t).^2);

% ģ�Ͳ�������
dx=2;dy=2;v=1000;%��������������ٶ�
x=-100:dx:100;y=-100:dy:100;%��������
x0=0;y0=0;t0=0.08;%�ڵ㼤���ص��ʱ��
m=length(x);n=length(y);

%��ֵ������FCT��Ƶɢ
data0=zeros(m,n);data1=data0;data2=data0;
m0=find(x==x0);n0=find(y==y0);l0=round(t0/dt)+1;
for k=1:l0
    if k<=length(t);
       data1(m0,n0)=R(k);
    end
    for i=2:m-1
        for j=2:n-1
            data2(i,j)=2*data1(i,j)-data0(i,j)+(v*dt)^2*...
               ((data1(i+1,j)-2*data1(i,j)+data1(i-1,j))/dx^2+...
               (data1(i,j+1)-2*data1(i,j)+data1(i,j-1))/dy^2);
        end
    end
%     data2 = FCTforAW(data0,data1,data2);
    data0=data1;data1=data2;
end

%���ݳ�ͼ
[X,Y]=meshgrid(x,y);
surf(X,Y,data2);shading interp;
h=title('Numerical Stimulation without FCT for Acoustic Wave');
set(h,'FontSize',14,'Color','b');
xlabel('X(m)','FontSize',14,'Color','b');
ylabel('Y(m)','FontSize',14,'Color','b');
view(0,90);%colormap('gray');
axis square;

toc;