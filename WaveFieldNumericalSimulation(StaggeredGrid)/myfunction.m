x=1:100;
y=1:100;
for i=1:length(x)
    for j=1:length(y)
        z(i,j)=-(12+round(10*sin(5*pi*x(i)/100)*sin(5*pi*y(j)/100)+0.5));
    end
end
x1=1:10:100;
y1=1:10:100;
for i=1:length(x1)
    for j=1:length(y1)
        z1(i,j)=-(12+round(10*sin(5*pi*x1(i)/100)*sin(5*pi*y1(j)/100)+0.5));
    end
end
[Y,X]=meshgrid(y,x);
surf(X',Y',z');
shading interp;colormap gray;
axis square;axis tight;grid on
xlabel('x(网格点数)');ylabel('y(网格点数)');zlabel('高程');
for i=1:length(x1)
    for j=1:length(y1)
        hold on;
        plot3(x1(i), y1(j), z1(i,j),'r.-');
    end
end
legend('地表','检波点');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot([0:100],zeros(1,101));
hold on;
plot([0:100],zeros(1,101)+50);
hold on;
plot([0:100],zeros(1,101)+80);
hold on;
plot([0:100],zeros(1,101)+100);
hold on;
plot(zeros(1,101),[0:100]);
hold on;
plot(zeros(1,101)+100,[0:100]);
hold on;
fill([[0:100],fliplr([0:100])],[zeros(1,101),fliplr(zeros(1,101)+50)],[1 1 1]);
hold on;
fill([[0:100],fliplr([0:100])],[zeros(1,101)+50,fliplr(zeros(1,101)+80)],[.5 1 1]);
hold on;
fill([[0:100],fliplr([0:100])],[zeros(1,101)+80,fliplr(zeros(1,101)+100)],[.5 .5 .5]);
hold on;
text(45,25,'1200m/s');
hold on;
text(45,65,'2000m/s');
hold on;
text(45,90,'3500m/s');
axis ij;
xlabel('x(网格点数)');ylabel('z(网格点数)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x=1:100;
y=1:100;
for i=1:length(x)
    for j=1:length(y)
        z(i,j)=-(12+round(10*sin(5*pi*x(i)/100)*sin(5*pi*y(j)/100)+0.5));
    end
end
x1=1:10:100;
y1=1:10:100;
for i=1:length(x1)
    for j=1:length(y1)
        z1(i,j)=-(12+round(10*sin(5*pi*x1(i)/100)*sin(5*pi*y1(j)/100)+0.5));
    end
end
x2=4:40:400;
y2=4:40:400;
x0=50;y0=50;
z0=-(12+10*sin(5*pi*x0/100)*sin(5*pi*y0/100)+3);
plot3(x0*4, y0*4, z0,'k*-');
hold on;
for i=1:length(x2)
    for j=1:length(y2)
        plot3(x2(i), y2(j), z1(i,j),'r.-');
        hold on;
    end
end
legend('检波点','炮点');
[Y,X]=meshgrid(4:4:400,4:4:400);
surf(X',Y',z');
hold on;
line([400 400],[400 400],[z(100,100),-200]);
hold on;
line([400 400],[0 0],[z(100,1),-200]);
hold on;
line([0 0],[400 400],[z(1,100),-200]);
hold on;
fill3([0,400,400,0],[0,0,400,400],[-200,-200,-200,-200],[.5 1 1]);
hold on;
fill3([0,400,400,0],[0,0,400,400],[-400,-400,-400,-400],[.5 .5 .5]);
hold on;
fill3([400,400,400,400],[0,0,400,400],[-200,-320,-320,-200],[.5 1 1]);
hold on;
fill3([0,0,0,0],[0,0,400,400],[-200,-320,-320,-200],[.5 1 1]);
hold on;
fill3([400,400,400,400],[0,0,400,400],[-320,-400,-400,-320],[.5 .5 .5]);
hold on;
fill3([0,0,0,0],[0,0,400,400],[-320,-400,-400,-320],[.5 .5 .5]);
hold on;
fill3([0,400,400,0],[0,0,0,0],[-200,-200,-320,-320],[.5 1 1]);
hold on;
fill3([0,400,400,0],[400,400,400,400],[-200,-200,-320,-320],[.5 1 1]);
hold on;
fill3([0,400,400,0],[0,0,0,0],[-320,-320,-400,-400],[.5 .5 .5]);
hold on;
fill3([0,400,400,0],[400,400,400,400],[-320,-320,-400,-400],[.5 .5 .5]);
hold on;
text(200,400,-150,'1200m/s');
hold on;
text(200,400,-260,'2000m/s');
hold on;
text(200,400,-360,'3500m/s');
%shading interp;colormap gray;
axis square;axis tight;grid on
xlabel('x(m)');ylabel('y(m)');zlabel('高程');


