% function   xout=ADMM_FidelityL2_TV(S,b,Init)
function   xout=ADMM_FidelityL2_TV(~,~,~)


% beta*||Hx||_{1}  +\frac{1}{2}||Sx-b||_{2}^{2}
%b Lidar data points
%S measurement matrix
%x depth image to recover
%u side information (either radar or lidar) transformed data points


% csvwrite('./input/S.csv', S);
% csvwrite('./input/b.csv', b);
% csvwrite('./input/Init.csv', Init);
% csvwrite('./input/eigDtD.csv', eigDtD);
S = csvread('./input/S.csv');
b = csvread('./input/b.csv');
Init = csvread('./input/Init.csv');
eigDtD = csvread('./input/eigDtD.csv');

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
% eigDtD = abs(fftn([1 -1],  [rows cols])).^2 + abs(fftn([1 -1]', [rows cols])).^2;


% parameter settings
beta    = 0.001;
%beta       =0.1;
max_itr  = 4000;
tol      = 1e-5;


% initiate intermediate variables - y1, y2, y3
y2 = zeros(size(r));
y3 = zeros(size(v));


%u=D(u);
% for itr=1:max_itr
for itr=1:max_itr
    fprintf('Iteration-%d\n',itr);
    
    xold = x;
%     csvwrite('xold.csv', xold);
    
    % solving the x-subproblem
    rhs   = S.*b  +rho.*r - (y2);
%     csvwrite('Sb.csv', S.*b);
%     csvwrite('rhor.csv', rho.*r);
%     csvwrite('rhs1.csv', rhs);
    x      = (1./(rho+ (S.^2))).*rhs;
%     csvwrite('x.csv', x);
   % x=max(x,0);
    
    %subproblem r
%     csvwrite('y3x1.csv', y3(:,:,1));
%     csvwrite('y3y2.csv', y3(:,:,2));
%     temp = rho.*v;
%     csvwrite('rho.xvx.csv', temp(:,:,1));
%     csvwrite('rho.xvy.csv', temp(:,:,2));
%     temp = y3- rho.*v;
%     csvwrite('y3x_rho.xvx.csv', temp(:,:,1));
%     csvwrite('y3y_rho.xvy.csv', temp(:,:,2));
    rhs   = y2 +rho.*x  -Dt(y3- rho.*v);
%     csvwrite('Dt.csv', Dt(y3- rho.*v));
%     csvwrite('rhs2.csv', rhs);
    r= real( ifft2( fft2(rhs)./(rho+rho.*eigDtD)) );
%     csvwrite('r.csv', r);
   

    %subproblem v
    tempv=D(r);
%     csvwrite('Dr1.csv', tempv(:,:,1));
%     csvwrite('Dr2.csv', tempv(:,:,2));
    tempv1=tempv+(y3./rho2);
%     csvwrite('abstempv1.csv', abs(tempv1(:,:,1)));
%     csvwrite('abstempv2.csv', abs(tempv1(:,:,2)));
%     csvwrite('abstempv1betarho2.csv', abs(tempv1(:,:,1))-(beta./rho2));
%     csvwrite('abstempv2betarho2.csv', abs(tempv1(:,:,2))-(beta./rho2));
%     tmp =  max(abs(tempv1)-(beta./rho2),0);
%     csvwrite('maxabstempv1betarho2.csv', tmp(:,:,1));
%     csvwrite('maxabstempv2betarho2.csv', tmp(:,:,2));
    v=sign(tempv1).*max(abs(tempv1)-(beta./rho2),0);
%     csvwrite('vx.csv', v(:,:,1));
%     csvwrite('vy.csv', v(:,:,2));
    
    
    % update multipliers
   % y1          = y1 + rho*(S.*x - b);
    y2          = y2 + rho2*(x - r);
%     csvwrite('xr.csv', x - r);
%     csvwrite('rho2xr.csv', rho2*(x - r));
%     csvwrite('y2.csv', y2);
    y3          = y3 + rho2*(tempv-v);
%     tmp = tempv-v;
%     csvwrite('tempvvx.csv', tmp(:,:,1));
%     csvwrite('tempvvy.csv', tmp(:,:,2));
%     csvwrite('y3x.csv', y3(:,:,1));
%     csvwrite('y3y.csv', y3(:,:,2));

    
    
%     relchg = norm(x - xold)/norm(x);
%     
%     fprintf('%3g \t %3.5e \n', itr, relchg);
    
    
    %     if (relchg<=tol)&&(itr>1)
    %   %      break;
    %     end
    %
%     csvwrite('x.csv', x);
end
xout=x;