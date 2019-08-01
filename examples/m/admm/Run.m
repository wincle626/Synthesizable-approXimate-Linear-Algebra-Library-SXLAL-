clear
clc


%co ordinates (X1,Y1) and (X2,Y2). Select region of interest in the depth image
X1=550;
X2=650;

Y1=X1+100;
Y2=X2+100;
%downsampling parameters
downsample_method = 'bicubic';
uf = 2;

% read high-res rgb and depth
disp_gt = double(imread('depth0000.png'));
gray_gt_ = im2double((imread('intensity0000.png')));


% figure(1)
% imagesc(disp_gt)
% title('Original Depth')

%Downsample and normalise the data
smin = min(disp_gt(:));
smax = max(disp_gt(:));
[M, N] = size(disp_gt);
[Mmb, Nmb] = size(gray_gt_);
dd = [Mmb - M; Nmb - N];
gray_gt = gray_gt_(dd(1)/2+1:end-dd(1)/2, dd(2)/2+1:end-dd(2)/2);
disp_res_n = imresize(disp_gt, 1/(2^uf), downsample_method);
ours = zeros(M,N);
ours(2^uf-(2^uf/2):2^uf:end, 2^uf-(2^uf/2):2^uf:end) = disp_res_n;
d_min = min(ours(ours>0));
d_max = max(ours(ours>0));
disp_res_n_norm = (ours-d_min)/(d_max-d_min);


%Calculate the tensor


%Determine the non zero points
weights = zeros(M,N);
weights(ours > 0) = 1;




% figure(2)
% imagesc(disp_gt(X1:Y1,X2:Y2))
% title('True Depth Region to be processed')
% figure(3)
% imagesc(weights(X1:Y1,X2:Y2).*disp_res_n_norm(X1:Y1,X2:Y2))
% title('Sampled Depth Image')

%% Run the admm code
% csvwrite('S.csv', weights(X1:Y1,X2:Y2));
% csvwrite('b.csv', disp_res_n_norm(X1:Y1,X2:Y2));
% csvwrite('Init.csv', disp_res_n_norm(X1:Y1,X2:Y2));

xout=ADMM_FidelityL2_TV(weights(X1:Y1,X2:Y2),disp_res_n_norm(X1:Y1,X2:Y2),disp_res_n_norm(X1:Y1,X2:Y2));


disp('ADMM Optimisation Finished')

%%  Run cvx

R=100;  %Rows R+1
C=100;  %Columns C+1
N=(R+1)*(C+1); %Total number of points analysed
% Difference matrix - along the rows
Dy=toeplitz([-1,zeros(1,N-2),1],[-1,1,zeros(1,N-2)]);
% Difference matrix - along the columns
Row=[-1,zeros(1,N-1)];
Row(C+2)=1;
Col=[-1,zeros(1,N-1)];
Col(N-R)=1;
Dx=toeplitz(Col,Row);

%Combined Difference matrix
D=[Dy;Dx];
b=vec(disp_res_n_norm(X1:Y1,X2:Y2));                %Lidar measurements
A=vec(weights(X1:Y1,X2:Y2)).*eye(N);          %Sensing Matrix



%Construct the Anistropic tensor

beta    = 0.001;

%TV with anistropic weighting + HTV
cvx_begin quiet
variable xWTV(N,1)
minimize(0.5*sum_square(A*xWTV - A*b(:)) + beta*norm(D*xWTV,1)  )
cvx_end




disp('CVX Optimisation Finished')

%%  Print Results


fprintf('Relative error (ADMM/CVX):           %f\n', norm(xWTV-vec(xout),2)/norm(xWTV,2))


