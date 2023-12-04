function s = contourdata(c)
% CONTOURDATA ���ݸ����ĵȸ�����Ϣ���صȸ�������
% level �ȸ��ߵȼ�
% number �ȸ���������
% x �ȸ��ߺ�����
% y �ȸ���������
% isopen �ȸ����Ƿ�պ�
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

