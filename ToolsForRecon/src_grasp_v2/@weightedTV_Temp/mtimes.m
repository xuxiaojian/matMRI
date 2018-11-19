function res = mtimes(a,b)

if a.adjoint
    res = adjDz(b,a.w); % wTV
else
    % res = b(:,:,[2:end,end]) - b;
    % wTV block begin
    res = zeros(size(b));
    for k = 1:size(b,3)-1
        res(:,:,k) = a.w(k)*(b(:,:,k+1)-b(:,:,k));
    end
    % wTV block end
end

function y = adjDz(x,w) % wTV
% y= x(:,:,[1,1:end-1]) - x;
% y(:,:,1) = -x(:,:,1);
% y(:,:,end) = x(:,:,end-1);

% wTV block begin
y = zeros(size(x));
y(:,:,1) = -w(1)*x(:,:,1);
y(:,:,end) = w(end)*x(:,:,end-1);
for k = 2:size(x,3)-1
    y(:,:,k) = w(k-1)*x(:,:,k-1)-w(k)*x(:,:,k);
end
% wTV block end