clear;
clc;
close all;
a = 1;
ff = 1;
% x1 = [0.1; 0.2]; x2 = [0.1; -0.1]; x3 = [0.3; 0.2]; x4 = [-0.2; 0.2]; x5 = [3; 0.2]; x6 = [-0.2; 1.2];
x1 = [0.1; 0.5]; x2 = [0.1; -0.1]; x3 = [0.3; 0.1]; x4 = [-0.2; 0.2]; x5 = [3; -0.2]; x6 = [-0.2; 1.2];
x1 = [2.5; -2.5]; x2 = [-1.5; 2]; x3 = [-2; -3]; x4 = [-2; -2]; x5 = [1.5; 1.5]; x6 = [2; -3]; x7 = [2; -3];

% x1 = [1; 2]; x2 = [-0.5; -0.1]; x3 = [0.3; 2]; x4 = [0.8; 0.2]; x5 = [-3; -2];
N = 141;
A = [0.995 -0.194;0.194 0.995];
% A=[0.995,0.09983;-0.09983,0.995];
% A=[0,1;0.1,1.05];
% A=[0,1;-1,1];
% E = [0.0621; -0.6784;0.1655;0.4247];
% B = [-1;1;1;3];
% A=[0.9801,-0.1987;0.1987,0.9801];
B = [1;1];

% A = [1 1; 0 0.9];
% B = [1;0.0];

K = [0.0937; 0.073];
% K= [0.085580893593885;0.067038484978488];

% K =[0.132416429211607;0.110656332931837];


% t = 1:1:N;
% z = 1:1:N-1;
% figure(3)
% plot(t,x1)
% figure(4)
% plot(z,u1)

% x_bar(:,1) = x1(:,1);
for kk = 1:N-1
%     if kk == 1
%         x_bar(:,kk)= A*[0.1; 1; 1] + B*K(:,kk)'*[0.1; 1; 1] + E*L(:,kk)'*[0.1; 1; 1];
%     else
    epslon1(:,kk)=1*(x1(:,kk)-x2(:,kk))+1*(x1(:,kk)-x3(:,kk)); 
    epslon2(:,kk)=1*(x2(:,kk)-x1(:,kk))+1*(x2(:,kk)-x3(:,kk))+1*(x2(:,kk)-x4(:,kk)); 
    epslon3(:,kk)=1*(x3(:,kk)-x1(:,kk))+1*(x3(:,kk)-x2(:,kk))+1*(x3(:,kk)-x4(:,kk)); 
    epslon4(:,kk)=1*(x4(:,kk)-x2(:,kk))+1*(x4(:,kk)-x3(:,kk)); 
    epslon5(:,kk)=1*(x5(:,kk)-x4(:,kk)); 

    u1(:,kk) = K'*epslon1(:,kk); u2(:,kk) = K'*epslon2(:,kk); u3(:,kk) = K'*epslon3(:,kk); 
    u4(:,kk) = K'*epslon4(:,kk); u5(:,kk) = K'*epslon5(:,kk); 

    x1(:,kk+1)= A*x1(:,kk) - B*u1(:,kk);
    x2(:,kk+1)= A*x2(:,kk) - B*u2(:,kk);
    x3(:,kk+1)= A*x3(:,kk) - B*u3(:,kk);
    x4(:,kk+1)= A*x4(:,kk) - B*u4(:,kk);
    x5(:,kk+1)= A*x5(:,kk) - B*u5(:,kk);
%     end
end
% figure(2)
% plot(t,x_bar)
epslon1(:,N)=1*(x1(:,N)-x2(:,N))+1*(x1(:,N)-x3(:,N)); 
epslon2(:,N)=1*(x2(:,N)-x1(:,N))+1*(x2(:,N)-x3(:,N))+1*(x2(:,N)-x4(:,N)); 
epslon3(:,N)=1*(x3(:,N)-x1(:,N))+1*(x3(:,N)-x2(:,N))+1*(x3(:,N)-x4(:,N)); 
epslon4(:,N)=1*(x4(:,N)-x2(:,N))+1*(x4(:,N)-x3(:,N)); 
epslon5(:,N)=1*(x5(:,N)-x4(:,N));

sumU1 = sum(u1.^2);
sumU2 = sum(u2.^2);
sumU3 = sum(u3.^2);
sumU4 = sum(u4.^2);
sumU5 = sum(u5.^2);

sum_of_squares_epslon1 = sum(sum(epslon1(N-20:N-1).^2, 1));
sum_of_squares_epslon2 = sum(sum(epslon2(N-20:N-1).^2, 1));
sum_of_squares_epslon3 = sum(sum(epslon3(N-20:N-1).^2, 1));
sum_of_squares_epslon4 = sum(sum(epslon4(N-20:N-1).^2, 1));
sum_of_squares_epslon5 = sum(sum(epslon5(N-20:N-1).^2, 1));

sum_of = (sum_of_squares_epslon1 + sum_of_squares_epslon2+ sum_of_squares_epslon3+sum_of_squares_epslon4+sum_of_squares_epslon5)/20


sum_IAE = sum(sum(abs(epslon1(N-20:N-1))))+sum(sum(abs(epslon2(N-20:N-1))))+sum(sum(abs(epslon3(N-20:N-1))))+...
    sum(sum(abs(epslon4(N-20:N-1))))+sum(sum(abs(epslon5(N-20:N-1))));

sum_IAE_ave = sum_IAE/5;

z = 0:1:N-2;
figure(1)
plot(z,u1)
hold on;
plot(z,u2)
hold on;
plot(z,u3)
hold on;
plot(z,u4)
hold on;
plot(z,u5)

xlabel('Time Step'),ylabel('$u_{i}$','interpreter','latex');
h=legend('$u_{1}$','$u_{2}$','$u_{3}$','$u_{4}$','$u_{5}$');
set(h,'interpreter','latex','FontName','Times New Roman','FontSize',14);
set(gca,'FontSize',14,'Fontname', 'Times New Roman');


figure(2)
% set(gcf,'unit','centimeters','position',[10,10,8,5])    % 图形窗口在电脑屏幕上的位置和尺寸[左 下 宽 高]
% linewidth_line = 1.2;      % 图形线条宽度
% markersize = 2.5;          % 图形标记点大小
% linewidth_gca = 0.7;      % 横纵坐标轴宽度
% fontsize_gca = 7;           % 横纵坐标轴刻度字体大小
% fontsize_label = 9;         % 横纵坐标轴字体大小
% fontsize_legend = 7;      % 图例字体大小
% plot(0:N-1, x0(1,:))
% hold on;

plot(0:N-1, x1(1,:))
hold on;
plot(0:N-1, x2(1,:))
hold on;
plot(0:N-1, x3(1,:))
hold on;
plot(0:N-1, x4(1,:))
hold on;
plot(0:N-1, x5(1,:))

xlabel('Time Step'),ylabel('$x_{1,i}$','interpreter','latex');
h=legend('$x_{1,1}$','$x_{1,2}$','$x_{1,3}$','$x_{1,4}$','$x_{1,5}$');
set(h,'interpreter','latex','FontName','Times New Roman','FontSize',14);
set(gca,'FontSize',14,'Fontname', 'Times New Roman');

figure(3)
% set(gcf,'unit','centimeters','position',[10,10,8,5])    % 图形窗口在电脑屏幕上的位置和尺寸[左 下 宽 高]
% linewidth_line = 1.2;      % 图形线条宽度
% markersize = 2.5;          % 图形标记点大小
% linewidth_gca = 0.7;      % 横纵坐标轴宽度
% fontsize_gca = 7;           % 横纵坐标轴刻度字体大小
% fontsize_label = 9;         % 横纵坐标轴字体大小
% fontsize_legend = 7;      % 图例字体大小
% plot(0:N-1, x0(1,:))
% hold on;

plot(0:N-1, x1(2,:))
hold on;
plot(0:N-1, x2(2,:))
hold on;
plot(0:N-1, x3(2,:))
hold on;
plot(0:N-1, x4(2,:))
hold on;
plot(0:N-1, x5(2,:))

xlabel('Time Step'),ylabel('$x_{2,i}$','interpreter','latex');
h=legend('$x_{2,1}$','$x_{2,2}$','$x_{2,3}$','$x_{2,4}$','$x_{2,5}$');
set(h,'interpreter','latex','FontName','Times New Roman','FontSize',14);
set(gca,'FontSize',14,'Fontname', 'Times New Roman');


figure(4)

plot(0:N-1, epslon1(1,:))
hold on;
plot(0:N-1, epslon2(1,:))
hold on;
plot(0:N-1, epslon3(1,:))
hold on;
plot(0:N-1, epslon4(1,:))
hold on;
plot(0:N-1, epslon5(1,:))

figure(5)
plot3(x1(1,:),x1(2,:),0:N-1)
hold on;
plot3(x2(1,:),x2(2,:),0:N-1)
hold on;
plot3(x3(1,:),x3(2,:),0:N-1)
hold on;
plot3(x4(1,:),x4(2,:),0:N-1)
hold on;
plot3(x5(1,:),x5(2,:),0:N-1)


xlabel('$x_{1}$','interpreter','latex'),ylabel('$x_{2}$','interpreter','latex'),zlabel('Time Step');

figure(6)

plot(0:N-1, epslon1(2,:))
hold on;
plot(0:N-1, epslon2(2,:))
hold on;
plot(0:N-1, epslon3(2,:))
hold on;
plot(0:N-1, epslon4(2,:))
hold on;
plot(0:N-1, epslon5(2,:))