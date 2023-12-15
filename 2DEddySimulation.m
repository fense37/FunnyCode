clc,clear,close all
% 数值模拟
% 模拟一个反气旋的运动轨迹
%% 基本参数
g = 9.8; % 重力加速度
R = 6371000; % 地球半径 
H = 1000; % 水深
dx = 10000; % 10km网格x
dy = 10000; % 10km网格y
dt = 0.1; % 时间步长
Nx = 50;  % x轴大小
Ny = 50;  % y轴大小
Nt = 20000000; % 时间跨度，1year
omega = 2*pi/24/60/60;  %  地球自转角速度

%% 三维时空运动场
Hb = H*ones(Nx, Ny); % 水深场
u = zeros(Nx, Ny);
v = u; h = u;  % 初始化速度和海面高度场
x = (1:Nx)*dx; % x坐标轴
y = (1:Ny)*dy; % y坐标轴
[mx, my] = meshgrid(x ,y);
f = 2*omega*sin((my+500*dy)/R);  % 科氏力参数
% 初始值设置
h(:, :) = -0.01*exp((-(mx/max(x)-1/2).^2-(my/max(y)-1/2).^2)*100);
for i = 2:Nx-1
    for j = 2:Ny-1
        u(i,j) = -g/2/f(i,j)/dy*(h(i+1,j)-h(i-1,j));
        v(i,j) = g/2/f(i,j)/dx*(h(i,j+1)-h(i,j-1));
    end
end


%% 时间演化
% for k = 1: Nt-1
%     H = Hb + h(:, :, k);
%     for i = 2: Nx-1
%         for j = 2:Ny-1
%             h(i,j,k+1) = h(i,j,k) + dt*(-H(i,j)*((u(i, j+1, k)-u(i, j-1, k))/2/dx+(v(i+1, j, k)-v(i-1, j, k))/2/dy)...;
%                 -u(i,j,k)/2/dx*(h(i,j+1,k)-h(i,j-1,k))-v(i,j,k)/2/dy*(h(i+1,j,k)-h(i-1,j,k)));
%             u(i,j,k+1) = u(i,j,k) + dt*(f(i,j)*v(i,j,k)-g/2/dx*(h(i,j+1,k)-h(i,j-1,k))...
%                 -u(i,j,k)/2/dx*(u(i,j+1,k)-u(i,j-1,k))-v(i,j,k)/2/dy*(u(i+1,j,k)-u(i-1,j,k)));
%             v(i,j,k+1) = v(i,j,k) + dt*(-f(i,j)*u(i,j,k)-g/2/dy*(h(i+1,j,k)-h(i-1,j,k))...
%                 -u(i,j,k)/2/dx*(v(i,j+1,k)-v(i,j-1,k))-v(i,j,k)/2/dy*(v(i+1,j,k)-v(i-1,j,k)));
%         end
%     end
% end


% % 创建一个示例三维矩阵（假设为 10x10xN 大小的矩阵，其中 N 是时间步数）
% N = Nt;
% data = h;
% 
% % 提取矩阵的维度
% [x, y, z] = size(data);
% 
% % 创建一个 VideoWriter 对象
% videoFile = 'matrix_animation.mp4';
% videoObj = VideoWriter(videoFile, 'MPEG-4');
% open(videoObj);
% 
% % 创建一个图形窗口
% figure;
% 
% % 循环遍历时间步，绘制每一帧并将其写入视频
% for t = 1:1000:N
%     % 获取当前时间步的切片
%     slice = data(:, :, t);
%     sliceu = u(:,:,t);
%     slicev = v(:,:,t);
% 
%     % 绘制三维图形
%     contourf(mx,my,slice);
%     hold on;
%     quiver(mx,my,sliceu,slicev)
%     xlabel('X轴');
%     ylabel('Y轴');
%     zlabel('数值');
%     title(['三维矩阵随时间变化的动画 - 时间步 ', num2str(t*dt) 's']);
% 
%     % 添加颜色映射（可选）
%     colormap('jet');
%     colorbar;
%     caxis([-0.01 0.01])
% 
%     % 暂停一小段时间以便观察
%     pause(0.1);
% 
%     % 将当前帧写入视频
%     frame = getframe(gcf);
%     writeVideo(videoObj, frame);
% 
%     % 清空图形窗口，准备下一帧
%     clf;
% end
% 
% % 关闭视频文件
% close(videoObj);
 h1 = h; u1 = u; v1 = v;
% 创建一个 VideoWriter 对象
videoFile = 'matrix_animation.mp4';
videoObj = VideoWriter(videoFile, 'MPEG-4');
open(videoObj);
 
 
%% 减少内存消耗
for k = 1: Nt-1
    H = Hb + h;
    for i = 2: Nx-1
        for j = 2:Ny-1
            h1(i,j) = h(i,j) + dt*(-H(i,j)*((u(i, j+1)-u(i, j-1))/2/dx+(v(i+1, j)-v(i-1, j))/2/dy)...;
                -u(i,j)/2/dx*(h(i,j+1)-h(i,j-1))-v(i,j)/2/dy*(h(i+1,j)-h(i-1,j)));
            u1(i,j) = u(i,j) + dt*(f(i,j)*v(i,j)-g/2/dx*(h(i,j+1)-h(i,j-1))...
                -u(i,j)/2/dx*(u(i,j+1)-u(i,j-1))-v(i,j)/2/dy*(u(i+1,j)-u(i-1,j)));
            v1(i,j) = v(i,j) + dt*(-f(i,j)*u(i,j)-g/2/dy*(h(i+1,j)-h(i-1,j))...
                -u(i,j)/2/dx*(v(i,j+1)-v(i,j-1))-v(i,j)/2/dy*(v(i+1,j)-v(i-1,j)));
        end
    end
    h = h1;
    u = u1;
    v = v1;
    if mod(k,86400/dt/20)==0
    contourf(mx,my,h);
    hold on;
    quiver(mx,my,u,v)
    xlabel('X轴');
    ylabel('Y轴');
    zlabel('数值');
    title([' 时间步 ', num2str(k*dt/86400) 'day']);
    colormap(othercolor('BuDRd_12'));
    colorbar;
    caxis([-0.01 0.01])
    % 暂停一小段时间以便观察
    pause(0.1);

    % 将当前帧写入视频
    frame = getframe(gcf);
    writeVideo(videoObj, frame);

    % 清空图形窗口，准备下一帧
    clf;
    end
    
end
close(videoObj);