clear;clc;

%************************1:the surface definition****************************
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
xlabel('x(grids)');ylabel('y(grids)');zlabel('digital elevation');
for i=1:length(x1)
    for j=1:length(y1)
        hold on;
        plot3(x1(i), y1(j), z1(i,j),'r.-');
    end
end
legend('surface','detection position');

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
xlabel('x(grids)');ylabel('z(grids)');


