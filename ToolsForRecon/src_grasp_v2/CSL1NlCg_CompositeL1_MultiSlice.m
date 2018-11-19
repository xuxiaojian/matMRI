function x = CSL1NlCg_CompositeL1_MultiSlice(x0,param,param2,param3) % CompositeL1
%
% res = CSL1NlCg(param)
%
% Compressed sensing reconstruction of undersampled k-space MRI data
%
% L1-norm minimization using non linear conjugate gradient iterations
%
% Given the acquisition model y = E*x, and the sparsifying transform W,
% the pogram finds the x that minimizes the following objective function:
%
% f(x) = ||E*x - y||^2 + lambda * ||W*x||_1
%
% Based on the paper: Sparse MRI: The application of compressed sensing for rapid MR imaging.
% Lustig M, Donoho D, Pauly JM. Magn Reson Med. 2007 Dec;58(6):1182-95.
%
% Ricardo Otazo, NYU 2008
%

fprintf('\n Non-linear conjugate gradient algorithm')
fprintf('\n ---------------------------------------------\n')

% starting point
x=x0;

% line search parameters
maxlsiter = 150 ;
gradToll = 1e-3 ;
param.l1Smooth = 1e-15;
param2.l1Smooth = 1e-15;	% CompositeL1
param3.l1Smooth = 1e-15;	% Multislice
alpha = 0.01;
beta = 0.6;
t0 = 1 ;
k = 0;

% compute g0  = grad(f(x))
g0 = grad(x,param,param2,param3); % CompositeL1 % Multislice
dx = -g0;

% iterations
while(1)
    
    % Cihat's correction begins here...
    
    % % %     % backtracking line-search
    % % % 	f0 = objective(x,dx,0,param,param2); % CompositeL1
    % % % 	t = t0;
    % % %     f1 = objective(x,dx,t,param,param2); % CompositeL1
    % % % 	lsiter = 0;
    % % %
    % % % %     fprintf('lsiter = %d, f0 = %1.15f, f1 = %1.15f, t = %1.15f\n',lsiter,f0,f1,t);
    % % % 	while (f1 > f0 - alpha*t*abs(g0(:)'*dx(:))) && (lsiter<maxlsiter)
    % % % 		lsiter = lsiter + 1;
    % % % 		t = t * beta;
    % % % 		f1 = objective(x,dx,t,param,param2); % CompositeL1
    % % %         % disp(['lsiter = ' num2str(lsiter) ', f0 = ' num2str(f0) ', f1 = ' num2str(f1) ', t = ' num2str(t)])
    % % % %         fprintf('lsiter = %d, f0 = %1.15f, f1 = %1.15f, t = %1.15f\n',lsiter,f0,f1,t);
    % % % 	end
    % % %
    % % % 	if lsiter == maxlsiter
    % % % 		disp('Error - line search ...');
    % % % 		return;
    % % % 	end
    % % %
    % % % 	% control the number of line searches by adapting the initial step search
    % % % 	if lsiter > 2, t0 = t0 * beta;end
    % % % 	if lsiter<1, t0 = t0 / beta; end
    % % %     % disp(['t0 = ' num2str(t0)])
    % % % %     fprintf('t0 = %1.15f\n',t0);
    % % %
    % % %     % update x
    % % % 	x = (x + t*dx);
    % % %
    % % % 	% print some numbers
    % % %     if param.display,
    % % %         fprintf(' ite = %d, cost = %1.15f, t = %1.5f, t0 = %1.5f, mean(abs(x(:,:,1))) = %1.5e\n',k,f1,t,t0,mean(mean(abs(x(:,:,1)))));
    % % %     end
    % % %
    % % %     %conjugate gradient calculation
    % % % 	g1 = grad(x,param,param2); % CompositeL1
    % % % 	bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps);
    % % % 	g0 = g1;
    % % % 	dx =  - g1 + bk* dx;
    % % % 	k = k + 1;
    % % %
    % % %
    % % %     % disp(['bk = ' num2str(bk) ', norm(dx(:)) = ' num2str(norm(dx(:)))])
    % % % %     fprintf('bk = %1.15f, norm(dx(:)) = %1.15f\n',bk,norm(dx(:)));
    % % %
    % % % 	% stopping criteria (to be improved)
    % % % 	if (k > param.nite) || (norm(dx(:)) < gradToll), break;end
    
    if k<3
        t0 = beta*t0;
    end
    
    t = t0;
        
    % update x
    x = (x + t*dx);
    
    % print some numbers
    if param.display,
        fprintf('CS iteration = %d\n',k);
    end
    
    %conjugate gradient calculation
    g1 = grad(x,param,param2,param3); % CompositeL1 % Multislice
    bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps);
    g0 = g1;
    dx =  - g1 + bk* dx;
    k = k + 1;
    
    % stopping criteria (to be improved)
    if k > param.nite
        break;
    end
    
    % Cihat's correction ends here.
    
end
return;

function res = objective(x,dx,t,param,param2) %********************************** % CompositeL1

% [Ny,Nx,Nz,Nt] = size(x);

% L2-norm part
w=param.E*(x+t*dx)-param.y;
L2Obj=w(:)'*w(:);

% L1-norm part
if param.lambda
    N = size(x,1); % Assuming that x is square
    M = 2^floor(log2(N));
    D = (N-M)/2;
    w = param.W*(x(D+1:end-D,D+1:end-D,:,:)+t*dx(D+1:end-D,D+1:end-D,:,:));
    L1Obj = sum((conj(w(:)).*w(:)+param.l1Smooth).^(1/2));
else
    L1Obj = 0;
end

% CompositeL1 block begin
if param2.lambda
    w = param2.W*(x+t*dx);
    L1Obj2 = sum((conj(w(:)).*w(:)+param2.l1Smooth).^(1/2));
else
    L1Obj2 = 0;
end

if param3.lambda % Multislice
    w = param3.W*(x+t*dx);
    L1Obj3 = sum((conj(w(:)).*w(:)+param3.l1Smooth).^(1/2));
else
    L1Obj3 = 0;
end

% CompositeL1 block end

% objective function
res=L2Obj+param.lambda*L1Obj+param2.lambda*L1Obj2+param3.lambda*L1Obj3; % CompositeL1

% disp(['L2Obj = ' num2str(L2Obj)])
% disp(['L1Obj1 = ' num2str(L1Obj)])
% disp(['L1Obj2 = ' num2str(L1Obj2)])
% disp(['param.lambda = ' num2str(param.lambda)])
% disp(['param2.lambda = ' num2str(param2.lambda)])
% disp(['param.lambda*L1Obj1 = ' num2str(param.lambda*L1Obj)])
% disp(['param2.lambda*L1Obj2 = ' num2str(param2.lambda*L1Obj2)])

function g = grad(x,param,param2, param3)%*********************************************** % CompositeL1 % Multislice

% L2-norm part
L2Grad = 2.*(param.E'*(param.E*x-param.y));

% L1-norm part
if param.lambda
    N = size(x,1); % Assuming that x is square
    M = 2^floor(log2(N));
    D = (N-M)/2;
    w = param.W*x(D+1:end-D,D+1:end-D,:,:);
    L1Grad = param.W'*(w.*(w.*conj(w)+param.l1Smooth).^(-0.5));
    L1Grad = padarray(L1Grad,[D D 0 0]);
else
    L1Grad=0;
end

% CompositeL1 block begin
if param2.lambda
    w = param2.W*x;
    L1Grad2 = param2.W'*(w.*(w.*conj(w)+param2.l1Smooth).^(-0.5));
else
    L1Grad2=0;
end

if param3.lambda
    w = param3.W*x;
    L1Grad3 = param3.W'*(w.*(w.*conj(w)+param3.l1Smooth).^(-0.5));
else
    L1Grad3=0;
end
% CompositeL1 block end


% composite gradient
g=L2Grad+param.lambda*L1Grad+param2.lambda*L1Grad2+param3.lambda*L1Grad3; % CompositeL1

