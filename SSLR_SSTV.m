function [output_image B S J W] = SSLR_SSTV(oriData3_noise, opts)

% -------------------------------------------------------------------------------------------------------------
blocksize = opts.blocksize;
stepsize  = opts.stepsize;
lambda    = opts.lambda;
beta1      =opts.beta1;
tau       = opts.tau;
r1        = opts.r1;
r2        = opts.r2;
maxIter   = opts.maxIter;
tol       = opts.tol;
% the parameters in this part is the link between global and local iamge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
[M N p] = size(oriData3_noise);

R         =   M-blocksize+1;
C         =   N-blocksize+1;
rr        =   [1:stepsize:R];
% rr        =   [rr rr(end)+1:R];
if rr(end)<R;rr = [rr R];else end;
cc        =   [1:stepsize:C];
% cc        =   [cc cc(end)+1:C];
if cc(end)<R;cc = [cc C];else end;
row       =   length(rr);
column    =   length(cc);
% Index image
Idxx     =   (1:row*column);
Idxx     =   reshape(Idxx, row, column);

opts.rr   = rr;
opts.cc   = cc;
opts.row  = row;
opts.column = column;
opts.Idxx   = Idxx;
opts.imagesize = [M,N,p];
%% Initializing optimization variables
Yp = imageTopatch(oriData3_noise,opts); % store 3-D image patch by patch
Y=oriData3_noise;
weight = calweight(opts);
[M1,N1,p1] = size(Yp);
Xp = Yp; 
Jp = Xp; Sp = zeros(M1,N1,p1);S=zeros(M,N,p);Bp = zeros(M1,N1,p1); 
W = patchToimage(Jp,opts);
YY = zeros(M1,N1,p1); YX = zeros(M1,N1,p1); YS = zeros(M,N,p);
 YJ = zeros(M,N,p); 

u1 = zeros(M,N,p); u2 = zeros(M,N,p); u3 = zeros(M,N,p);
y1 = zeros(M,N,p); y2 = zeros(M,N,p); y3 = zeros(M,N,p);

% for the 3D-TV norm
% define operators
% in this problem H=1;beta = [1 1 0.5];
H =1;beta = [1 1 0.5];  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eigHtH      = abs(fftn(H, [M N p])).^2;
eigDtD      = abs(beta(1)*fftn([1 -1],  [M N p])).^2 + abs(beta(2)*fftn([1 -1]', [M N p])).^2;
if p>1
    d_tmp(1,1,1)= 1; d_tmp(1,1,2)= -1;
    eigEtE  = abs(beta(3)*fftn(d_tmp, [M N p])).^2;
else
    eigEtE = 0;
end
% Htg         = imfilter(g, H, 'circular');
[D,Dt]      = defDDt(beta);
%% 
rho = 1.5;
max_mu1 = 1e6;
mu = 1e-1;
sv =10;
%% 
out_value = [];
out_value.SSIM =[];
out_value.PSNR =[];
out_value.obj= [];

%% main loop
iter = 0;

while iter<maxIter
    iter = iter + 1;
   
%%     Update X & Xp
temp = (Yp - Sp - Bp + Jp)/2 + (YY - YX) /(2*mu);
for i =1:p1
%     Update Xp
  temp1 = temp(:,:,i);
  if  choosvd(p,sv) ==1
      [U, sigma, V] = lansvd(temp1, sv, 'L');
  else
      [U,sigma,V] = svd(temp1,'econ');
  end
    sigma = diag(sigma);
    svp = min(length(find(sigma>1/(2*mu))),r1);
    if svp<sv
        sv = min(svp + 1, p);
    else
       sv = min(svp + round(0.05*p), p);  
    end
  Xp(:,:,i) = U(:, 1:svp) * diag(sigma(1:svp) - 1/(2*mu)) * V(:, 1:svp)'; 
end

    %%  update S
 X=patchToimage(Xp,opts);
 B=patchToimage(Bp,opts); 
  temp_S=zeros(M,N,p);
for i =1:p
    temp_S(:,:,i) = (Y(:,:,i) - X(:,:,i) - B(:,:,i)) + YS(:,:,i)/(mu);
    temp2 = temp_S(:,:,i);
  if  choosvd(p,sv) ==1
      [U, sigma, V] = lansvd(temp2, sv, 'L');
  else
      [U,sigma,V] = svd(temp2,'econ');
  end
    sigma = diag(sigma);
    svp = min(length(find(sigma>beta1/(mu))),r2);
    if svp<sv
        sv = min(svp + 1, p);
    else
       sv = min(svp + round(0.05*p), p);  
    end
    S(:,:,i) = U(:, 1:svp) * diag(sigma(1:svp) - beta1/(mu)) * V(:, 1:svp)'; 
end
Sp=imageTopatch(S,opts); 
% [~,Sp]=imageTopatch(S,opts); 
%%    Update B & Bp
  temp_B = Yp - Xp - Sp + YY/mu;
  Bp_hat = max(temp_B - lambda/mu, 0);
  Bp = Bp_hat+min(temp_B + lambda/mu, 0);
%%    Update J & Jp
 temp_C = Xp + YX/mu;
 temp_D = W - YJ/mu;
 J      = update_M(temp_C,temp_D,weight,opts);
%  [~,Jp]     = imageTopatch(J,opts);
Jp     = imageTopatch(J,opts);
%%       Update W
temp = J + YJ/mu;
Htg   = imfilter(temp, H, 'circular');
rhs   = fftn(Htg + Dt(u1+(1/mu)*y1,  u2+(1/mu)*y2, u3+(1/mu)*y3));
eigA  = (mu/mu)*eigHtH + eigDtD + eigEtE;
W     = real(ifftn(rhs./eigA));

%% update U = [u1 u2 u3]
 [Df1 Df2 Df3] = D(W);

  v1 = Df1-(1/mu)*y1;
  v2 = Df2-(1/mu)*y2;
  v3 = Df3-(1/mu)*y3;
  u1 = max(v1 - tau/mu, 0); u1 = u1 + min(v1 + tau/mu,0);
  u2 = max(v2 - tau/mu, 0); u2 = u2 + min(v2 + tau/mu,0);
  u3 = max(v3 - tau/mu, 0); u3 = u3 + min(v3 + tau/mu,0);
 
%%
  leq1 = Yp - Xp - Sp- Bp;
  leq2 = Xp - Jp;
  leq3 = J - W;
  leq4 = Y-X-S-B;
%   leq4 = Sp - Qp;
%% stop criterion          
%     stopC = max(max(max(abs(leq1))),max(max(abs(leq2))));
 out_value.obj1(iter) = sum(abs(leq1(:)));
 out_value.obj2(iter) = sum(abs(leq2(:)));
 out_value.obj3(iter) = sum(abs(leq3(:)));
  out_value.obj4(iter) = sum(abs(leq4(:)));
stopC1 = max(abs(leq1(:)));
stopC2 = max(abs(leq2(:)));
stopC3 = max(abs(leq3(:)));
stopC4 = max(abs(leq4(:)));
stopC = max(max(max(stopC1,stopC2),stopC3),stopC4);

% disp(['iter ' num2str(iter) ',mu=' num2str(mu,'%2.1e')  ...
%            ',stopALM=' num2str(stopC,'%2.3e')]);
 
 if stopC<tol
   break;
 else
   YY = YY + mu*leq1;
   YX = YX + mu*leq2;  
   YJ = YJ + mu*leq3;
   YS = YS + mu*leq4;
   
   y1   = y1 + mu*(u1 - Df1);
   y2   = y2 + mu*(u2 - Df2);
   y3   = y3 + mu*(u3 - Df3);
   mu = min(max_mu1,mu*rho);
   
 end
%% output SSIM and PSNR values of each step
output_image = J;
% output_image = S;
% figure,imshow(X(:,:,i),[])

end


end

function [D,Dt] = defDDt(beta)
D  = @(U) ForwardD(U, beta);
Dt = @(X,Y,Z) Dive(X,Y,Z, beta);
end

function [Dux,Duy,Duz] = ForwardD(U, beta)
frames = size(U, 3);
Dux = beta(1)*[diff(U,1,2), U(:,1,:) - U(:,end,:)];
Duy = beta(2)*[diff(U,1,1); U(1,:,:) - U(end,:,:)];
Duz(:,:,1:frames-1) = beta(3)*diff(U,1,3); 
Duz(:,:,frames)     = beta(3)*(U(:,:,1) - U(:,:,end));
end

function DtXYZ = Dive(X,Y,Z, beta)
frames = size(X, 3);
DtXYZ = [X(:,end,:) - X(:, 1,:), -diff(X,1,2)];
DtXYZ = beta(1)*DtXYZ + beta(2)*[Y(end,:,:) - Y(1, :,:); -diff(Y,1,1)];
Tmp(:,:,1) = Z(:,:,end) - Z(:,:,1);
Tmp(:,:,2:frames) = -diff(Z,1,3);
DtXYZ = DtXYZ + beta(3)*Tmp;
end

function [ weight] = calweight(par )
% this function is to calculate the weight from 3-D image to blockes
% by Wei He 
%%
% [M1 N1 p1] = size(outpatch);
blocksize  = par.blocksize;
imagesize  = par.imagesize;

rr        = par.rr;
cc        = par.cc;
row       = par.row;
column    = par.column;
Idxx      = par.Idxx;
Numofpatch= row * column; 

weight   = zeros(imagesize);
for idx = 1:Numofpatch
    [rowidx,columnidx] = find(Idxx==idx); % find the location in index image
    i = rr(rowidx); j = cc(columnidx);    % find the location in original image
    weight(i:i+blocksize-1,j:j+blocksize-1,:)= weight(i:i+blocksize-1,j:j+blocksize-1,:)+ 1;
end  
end