function s = contourdata(c)
% CONTOURDATA 根据给定的等高线信息返回等高线数据
% level 等高线等级
% number 等高线像素数
% x 等高线横坐标
% y 等高线纵坐标
% isopen 等高线是否闭合
if nargin < 1 || ~isfloat(c) || size(c, 1) ~= 2 || size(c, 2) < 4
    error('CONTOURDATA:rhs',...
    'Input Must be the 2-by-M matrix C.')
end
tol = 1e-12;
k = 1;
col = 1;
while col < size(c, 2)
    s(k).level = c(1, col);
    s(k).number = c(2, col);
    idx = col + 1: col + c(2, col);
    s(k).x = c(1, idx);
    s(k).y = c(2, idx);
%     s(k).isopen = abs( sum(diff( c(1, idx(1:end)) )) ) > tol ||...
%         abs( sum(diff( c(2, idx(1: end))))) > tol;
    s(k).isopen = abs( diff( c(1, idx([1 end]) )) ) > tol ||...
        abs( diff( c(2, idx([1 end])))) > tol;
    k = k + 1;
    col = col + c(2, col) + 1;
end

