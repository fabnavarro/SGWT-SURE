%% Run wavelets on the Pittsburgh Synthetic Data
clear
%% 

addpath(genpath('DiffusionWavelets'))
addpath MTSG_Toolbox-master
addpath(genpath('DiffusionGeometry'))

Startup_DiffusionGeometry
Startup_DiffusionWavelets
MTSG_Path

%% get data


%% load pitt graph
load pittsburgh/pittsburgh.mat
clear ys
clear mu
%% load data from R
load pittsburgh/table_pitt_sig1
%load pittsburgh/table_pitt_sig2
%load pittsburgh/table_pitt_sig3

n=length(mu);
[D,edges]=preprocess_D(D,n);
base=norm(D,1);

y_val=mu;

%% Preprocessing
% James Sharpnack
load('JamesCode/PittsburghWavletBasis.mat') %got matrix B

% Mauro Maggioni
k=3;
[m,tmp]=size(edges);
    edges1=edges(:,1);edges2=edges(:,2);
    if tmp==3
        w=edges(:,3);
    else
        w=ones(m,1);
    end
    A=sparse(edges1,edges2,w,n,n);
    A=A+A';
    Deg=sum(A)';
    P=spdiags(1./sqrt(Deg),0,n,n);
    %symmetric diffusion operator
    T=P*A*P;
    Tree = DWPTree (T^k, 12, 1e-4, struct('Wavelets',true,'OpThreshold',1e-2,'GSOptions',struct('StopDensity',1,'Threshold',1e-2)));

% Jeff Irion

% no preprocessing needed 

%% Running experiments

sz=n;
rep=10;

xsh_list=zeros(n,sz,rep);% Sharpnack -  hard thresholding
xs_list=zeros(n,sz,rep);% Sharpnack -  soft thresholding
xm_list=zeros(n,sz,rep);% Maggioni
xi_list=zeros(n,sz,rep);% Irion
y_list=zeros(n,rep);

y_list=ys;


msesh_list=zeros(sz,rep);% Sharpnack -  hard thresholding
mses_list=zeros(sz,rep);% Sharpnack -  soft thresholding
msem_list=zeros(sz,rep);% Maggioni
msei_list=zeros(sz,rep);% Irion


%%
for j=1:rep
    y=y_list(:,j);
    
    [xm_list(:,:,j), ~] = maggioni_denoise(edges,y,k,Tree);
    %%
    [xi_list(:,:,j), ~] = irion_denoise(edges,y);
    %%
    
    
    msem_list(:,j)=sum(bsxfun(@minus,xm_list(:,:,j),y_val).^2)/n;
    msei_list(:,j)=sum(bsxfun(@minus,xi_list(:,:,j),y_val).^2)/n;
    
    
    yrot=B*y;
    lambdalist=sort(abs(yrot));
    for i=1:sz
        lambda=lambdalist(i);
        c=yrot;
        c(abs(c)<lambda)=0;%hard thresholding
        x=B'*c;
        xsh_list(:,i,j)=x;
        msesh_list(i,j)=norm(x-y_val)^2/n;
        snr_list_hard(i,j) = 20*log10(norm(y_val)/norm(y_val-x));
    end

    for i=1:sz
        c=yrot;
        lambda=lambdalist(i);
        c=max(0,c-lambda)...%soft thresholding
            + min(0,c+lambda); 
        %)(abs(yrot)<lambda*sqrt(2*log(n))); 
        x=B'*c;
        xs_list(:,i,j)=x;
        mses_list(i,j)=norm(x-y_val)^2/n;
        snr_list_soft(i,j) = 20*log10(norm(y_val)/norm(y_val-x));
    end

    
 tic   
    for i=1:sz

%n=length(y);
m=size(D,1);

O=construct_O(D,k-1);
scalefactor=base/norm(O,1);%mean(sum(abs(O),2);
lambda=scalefactor*lambdalist(i); rho=5*lambda;

lambda_list(i,j)=lambda;
tic
[ x ,history1] = gtf_admm_v2( y,edges,k-3,lambda,rho);
t_list(i,j)=toc;

% xm_list(:,i,j)=xm;
% xi_list(:,i,j)=xi;
% xs_list(:,i,j)=xs;


mse_list(i,j)=norm(x-y_val)^2/n;
snr_list(i,j) = 20*log10(norm(y_val)/norm(y_val-x));
reg_weight(i,j)=norm(O*x);
    end
    tftrend(j)  = toc;
end


nf = sqrt(mean(y_val.^2));

[aMSE,bMSE] = min(mse_list);
[cMSE,dMSE] = max(snr_list);
mean(aMSE)
mean(cMSE)

%msetosnr(nf,sqrt(mean(aMSE)))

[asharpSoft,bsharpSoft] = min(mses_list);
[csharpSoft,dsharpSoft] = max(snr_list_soft);
mean(asharpSoft)
mean(csharpSoft)
  
[asharpHard,bsharpHard] = min(msesh_list);
[csharpHard,dsharpHard] = max(snr_list_hard);
mean(asharpHard)
mean(csharpHard) 

[acoiff,bcoiff] = min(msei_list);
mean(acoiff)

%%

%%
hold off
errorbar(1:n,mean(msem_list,2),std(msem_list,[],2),'linewidth',1.5);
hold all
errorbar(1:n,mean(msei_list,2),std(msei_list,[],2),'linewidth',1.5);
errorbar(1:n,mean(mses_list,2),std(mses_list,[],2),'linewidth',1.5);
errorbar(1:n,mean(msesh_list,2),std(msesh_list,[],2),'linewidth',1.5);
%ylim([1.18^2/n,1.3^2/n])
%ylim([1.5e-4,2.5e-4])
xlabel('\lambda');
ylabel('MSE Validation error')
grid on
legend('Maggioni','Irion','Sharpnack','Sharpnack-Hard')


%k_list = 1:n;
%save('Wavelets_results_rep10.mat','k_list','xm_list','xi_list','xs_list','xsh_list','y_list');
