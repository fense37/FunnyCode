function psi = calc_stream_func(lon, lat, u, v)
%CALC_STREAM_FUNC Calculate stream function from velocity field
%   psi = calc_stream_func(lon, lat, u, v)
% inputs:
%   lon - longitude grid [n x 1]
%   lat - latitude grid [1 x m]
%   u   - zonal velocity [n x m]
%   v   - meridional velocity [n x m]
% outputs:
%   psi - stream function [n x m]
%----------------------------------------------------------------%

    % Ensure inputs are consistent
    if size(u, 1) ~= length(lon) || size(u, 2) ~= length(lat)
        error('Dimensions of u do not match lon and lat grids.');
    end
    if size(v, 1) ~= length(lon) || size(v, 2) ~= length(lat)
        error('Dimensions of v do not match lon and lat grids.');
    end
    psi = calc_stream_func_relax(lon, lat, u, v);

end

function psi = calc_stream_func_relax(lon, lat, u, v)
% calc stream function based on Richardson method (Richardson, 1920)
%   psi = calc_stream_func_rich(lon, lat, u, v)
% method:
%   1. fix boundary conditions
%   2. calculate the curl and divergence of the velocity field
%   3. solve the Poisson equation using relaxation method
%   4. return the stream function
    eps_psi = 50; % convergence criterion
    [dx,dy]=dxdy(lat,lon);
    dx = mean(mean(dx));
    dy = mean(mean(dy));
    % out of boundary is negative
    oiintVn = sum(u(1,:))*dy + sum(v(:,1))*dx + sum(u(end,:))*dy + sum(v(:,end))*dx;
    oiintaVn= sum(abs(u(1,:)))*dy + sum(abs(v(:,1)))*dx + sum(abs(u(end,:)))*dy + sum(abs(v(:,end)))*dx;
    eps = -oiintVn/oiintaVn;
    % fix boundary speed
    u(1,:) = u(1,:) + eps*abs(u(1,:));
    u(end,:) = u(end,:) + eps*abs(u(end,:));
    v(:,1) = v(:,1) + eps*abs(v(:,1));
    v(:,end) = v(:,end) + eps*abs(v(:,end));
    % cal_stream_fun
    [~, uy] = calGrad2d(lon, lat, u);
    [vx, ~] = calGrad2d(lon, lat, v);
    zeta = vx - uy;
    % residual
    res = -zeta;
    denominator = 2*(1/dx^2 + 1/dy^2);
    psi_n = zeros(size(u));
    psi_np1 = ones(size(u))*1e4;
    it_no = 0;
    eri = nanmean(abs(psi_np1 - psi_n),'all');
    while eri > eps_psi
        psi_n = psi_np1;
        psi_np1 = psi_n + (1/denominator) * res;
        psi_im1 = [psi_n(:,2:end)  psi_n(:,end)+v(:,end)*dx];
        psi_ip1 = [psi_n(:,1)-v(:,1)*dx    psi_n(:,1:end-1)];
        psi_jm1 = [psi_n(2:end,:); psi_n(end,:)-u(end,:)*dy];
        psi_jp1 = [psi_n(1,:)+u(1,:)*dy;   psi_n(1:end-1,:)];
        res = (psi_im1 + psi_ip1)/dx^2 + (psi_jm1 + psi_jp1)/dy^2 - denominator*psi_n - zeta;
        it_no = it_no + 1;
        eri = nanmean(abs(psi_np1 - psi_n),'all');
    end
    psi = psi_np1;
end

function [dx,dy]=dxdy(lat,lon)
    R0=6.371*1.0e+6;
    [lat1,lon1]=meshgrid(lat,lon);
    [~,dy]=gradient(lon1);
    dy=dy.*2*pi*R0/360;
    [dx,~]=gradient(lat1);
    dx=dx.*2.*pi.*R0.*cosd(lat1)/360;
end
