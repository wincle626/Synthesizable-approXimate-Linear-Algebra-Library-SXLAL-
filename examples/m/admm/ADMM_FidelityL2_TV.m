function   xout=ADMM_FidelityL2_TV(S,b,Init)


% beta*||Hx||_{1}  +\frac{1}{2}||Sx-b||_{2}^{2}
%b Lidar data points
%S measurement matrix
%x depth image to recover
%u side information (either radar or lidar) transformed data points


[rows cols] = size(b);
[D, Dt]     = defDDt;


% initiate intermediate variables - r, h, v
x                   = Init;
r                 = zeros(rows,cols);
v                 = D(r);

% fixed internal parameters rho1, mu and gamma
rho = 0.1;
% rho2=0.1;


rho2=rho;

% pre-calculated fourier term
eigDtD = abs(fftn([1 -1],  [rows cols])).^2 + abs(fftn([1 -1]', [rows cols])).^2;


% parameter settings
beta    = 0.001;
%beta       =0.1;
max_itr  = 4000;
tol      = 1e-5;


% initiate intermediate variables - y1, y2, y3
y2 = zeros(size(r));
y3 = zeros(size(v));



%u=D(u);
for itr=1:max_itr
    xold = x;
    
    % solving the x-subproblem
    rhs   = S.*b  +rho.*r - (y2);
    x      = (1./(rho+ (S.^2))).*rhs;
   % x=max(x,0);
    
    %subproblem r
    rhs   = y2 +rho.*x  -Dt(y3- rho.*v);
    r= real( ifft2( fft2(rhs)./(rho+rho.*eigDtD)) );
   

    %subproblem v
    tempv=D(r);
    tempv1=tempv+(y3./rho2);
    v=sign(tempv1).*max(abs(tempv1)-(beta./rho2),0);
    
    
    % update multipliers
   % y1          = y1 + rho*(S.*x - b);
    y2          = y2 + rho2*(x - r);
    y3          = y3 + rho2*(tempv-v);

    
    
%     relchg = norm(x - xold)/norm(x);
%     
%     fprintf('%3g \t %3.5e \n', itr, relchg);
    
    
    %     if (relchg<=tol)&&(itr>1)
    %   %      break;
    %     end
    %
    
end
xout=x;