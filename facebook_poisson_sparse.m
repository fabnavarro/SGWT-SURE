load('facebook/facebook.txt')
addpath('RandomWalks')

facebook=facebook+1;%change to 1 based indexing
n=max(facebook(:));
[ D, edges ] = preprocess_D(facebook,n);

%% generate underlying signals
rng(15213)

headnodes=[0,107,1684,1912,3437,348,3980,414,686,698];
headnodes=headnodes+1;
%s=length(headnodes);
%idx=headnodes;

s=20;
idx=randperm(n);
idx=[idx(1:s), headnodes];
s=length(idx);

%0 based k
k=1;%O should be L
O=construct_O( D,k );

ss=sparse(n,1);
ss(idx)=randn(s,1)*100;
ss(idx)=ss(idx)-mean(ss(idx));

tol_pcg=1e-6;
pfun=cmg_sdd(O);
y=pcg(O,ss,tol_pcg,100,pfun);

%load Facebook_data_inhomogeneous;%Facebook_data_randomwalk.mat;
y=y/norm(y)*sqrt(n);%scale it up to order \sqrt(n)
%( Comment this out if it is the inhomogeneous data. but not homogeneous one)

%%

%% Get parameters of the experiments
rng(43651243)

SnR_list=[30,25,20,15,10,5,0,-5,-10];%in dB
% 10 log (signal_sigma^2/noise_sigma^2)
sigma_list=sqrt(10.^(-SnR_list/10));%/sqrt(n);
%sigma_list=[0.01,0.1,0.3,0.5,0.7,0.9];

lambdalist=-7:0.2:3;
lambdalist=exp(lambdalist);%exp(lambdalist);

sz=length(lambdalist);

%sigma_list=[0.01, 0.05];
szz=length(sigma_list);
mse_list=zeros(szz,1);

mse_listL=zeros(szz,1);
mse_listW=zeros(szz,1);
mse_list1=zeros(szz,1);
mse_list2=zeros(szz,1);



%filtered=zeros(n,szz,sz,5);


rep=5;

% the master result
mse_all=zeros(sz,szz,5,rep);
base=norm(D,1);

    mse_table=zeros(szz,rep); 
    mse_table1=zeros(szz,rep);
    mse_table2=zeros(szz,rep);
    mse_tableL=zeros(szz,rep);
    mse_tableW=zeros(szz,rep);

    %pre-generate all the noise
    Z=zeros(n,rep,szz);
    for i=1:length(sigma_list)
        Z(:,:,i)=randn(n,rep)*sigma_list(i);    
    end


%% k=0

for ii=1:rep
    for i=1:szz

        %% gtf k=0
        mse_list_inner=zeros(sz,1);
        reg_weight=zeros(sz,1);
        
        z=Z(:,ii,i);
        yhat=y+z;
        
        for j=1:sz % over all lambda
            
            
            
            m=size(D,1);
            scalefactor=base/norm(O,1);%mean(sum(abs(O),2);
            lambda=scalefactor*lambdalist(j); rho=10*lambda;

            [ x ,history1] = gtf_admm_v2( yhat,edges,k-1,lambda,rho);
            %[ x,history ] = gtf_proj_newton1( yhat,edges,k,lambda);

            mse_list_inner(j)=norm(x-y)^2/n;
            reg_weight(j)=norm(O*x);

        end
                loglog(lambdalist,mse_list_inner);
            hold on;
            mse_all(:,i,1,ii)=mse_list_inner;

        mse_list(i)=min(mse_list_inner);
               drawnow;
        %save('temp.mat');
    end
    %% To compare against the Laplacian smoother 
    for i=1:szz        
        
        
        z=Z(:,ii,i);
        yhat=y+z;
        
        mse_list_inner=zeros(sz,1);
        reg_weight=zeros(sz,1);

        for j=1:sz

            scalefactor=base/norm(O,1);%base/norm(O,'fro')^2;
            lambda=lambdalist(j)*scalefactor;
            %lambda=exp(j/2*log(lambdalist(j)));
            x=Laplacian_smoother(yhat,edges,lambda,0,0);

            mse_list_inner(j)=norm(x-y)^2/n;
        end
                    loglog(lambdalist,mse_list_inner);
            hold on;

                    mse_all(:,i,2,ii)=mse_list_inner;
        mse_listL(i)=min(mse_list_inner);
                   drawnow;
    end


%     %% To compare against James's Wavelets
% 
%     %load FacebookWavletBasis.mat %got matrix B
% 
% 
%     for i=1:szz
%         
%         z=Z(:,ii,i);
%         yhat=y+z;
%         
%         mse_list_inner=zeros(sz,1);
%         reg_weight=zeros(sz,1);
% 
%         for j=1:sz
% 
%             scalefactor=base/norm(O,1);%base/norm(O,'fro')^2;
%             lambda=lambdalist(j)*scalefactor;
% 
%                  yrot=B*yhat;
%          yrot=max(0,yrot-lambda*sqrt(2*log(n)))...%soft thresholding
%              + min(0,yrot+lambda*sqrt(2*log(n))); 
%          %)(abs(yrot)<lambda*sqrt(2*log(n))); 
%          x=B'*yrot;
% 
% 
%  
%             mse_list_inner(j)=norm(x-y)^2/n;
%         end
%                     loglog(lambdalist,mse_list_inner);
% 
% 
%             hold on;
% 
%             mse_all(:,i,3,ii)=mse_list_inner;
%         mse_listW(i)=min(mse_list_inner);
%                    drawnow;
%     end

    %% k=1
    k=1;
    for i=1:szz

        
        z=Z(:,ii,i);
        yhat=y+z;

        mse_list_inner=zeros(sz,1);
        reg_weight=zeros(sz,1);
        for j=1:sz
            m=size(D,1);
            scalefactor=base/norm(O,1);%mean(sum(abs(O),2);
            lambda=scalefactor*lambdalist(j); rho=10*lambda;

            [ x ,history1] = gtf( yhat,edges,k,lambda,rho);
            %[ x,history ] = gtf_proj_newton1( yhat,edges,k,lambda);

            mse_list_inner(j)=norm(x-y)^2/n;
            reg_weight(j)=norm(O*x);



        end
                loglog(lambdalist,mse_list_inner);
            hold on;

            mse_all(:,i,4,ii)=mse_list_inner;
        mse_list1(i)=min(mse_list_inner);
               drawnow;
        %save('temp.mat');
    end

    %% k=2
     k=2;
    for i=1:szz

                
        z=Z(:,ii,i);
        yhat=y+z;

        mse_list_inner=100*ones(sz,1);
        reg_weight=zeros(sz,1);
        for j=1:sz
            m=size(D,1);
            scalefactor=base/norm(O,1);%mean(sum(abs(O),2);
            lambda=scalefactor*lambdalist(j); rho=0.5*lambda;

            [ x ,history1] = gtf( yhat,edges,k,lambda,rho);
            %[ x,history ] = gtf_proj_newton1( yhat,edges,k,lambda);

            mse_list_inner(j)=norm(x-y)^2/n;
            reg_weight(j)=norm(O*x);




        end
                loglog(lambdalist,mse_list_inner);
            hold on;
            mse_all(:,i,5,ii)=mse_list_inner;
        mse_list2(i)=min(mse_list_inner);
            drawnow;
        
    end
 
    %%
    mse_table(:,ii) =mse_list;
    mse_table1(:,ii)= mse_list1;
    mse_table2(:,ii)=mse_list2;
    mse_tableL(:,ii)=mse_listL;
    mse_tableW(:,ii)=mse_listW;
    save('temp.mat');
end % end rep

%%

%save('Exp_Facebook_RandomWalk_last.mat')

%save('Exp_Facebook_results.mat','y','sigma_list','mse_list','mse_list2','mse_list3')
%%
h=figure;

% plot the relative denoising performance in dB, the maller the better of
% course.

% plot(-SnR_list,10*log(mse_list'./sigma_list.^2)/log(10)','ro-',...
%     -SnR_list,10*log(mse_list1'./sigma_list.^2)/log(10)','r*-',...
%     -SnR_list,10*log(mse_list2'./sigma_list.^2)/log(10)','r^-',...
%     -SnR_list,10*log(mse_listL'./sigma_list.^2)/log(10)','kx--',...
%     -SnR_list,10*log(mse_listW'./sigma_list.^2)/log(10)','bs-.','linewidth',2,'markersize',12);

plot(-SnR_list,10*log(mse_list')/log(10)','ro-',...
    -SnR_list,10*log(mse_list1')/log(10)','r*-',...
    -SnR_list,10*log(mse_list2')/log(10)','r^-',...
    -SnR_list,10*log(mse_listL')/log(10)','kx--',...%'linewidth',2,'markersize',12);
    linewidth',2,'markersize',12);

xlabel('Noise Level: -SnR in (dB)','fontsize',14)
ylabel('Denoising Performance in dB','fontsize',14)
ylabel('Denoised -SnR in dB','fontsize',14)
lg=legend('Trend filtering k=0','Trend filtering k=1',...
  'Trend filtering k=2','Laplacian smoothing');
set(lg,'fontsize',14,'location','best')
grid on;



