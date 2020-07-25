%% Simulating graph signal on real geographic graph
% Load a road network to use for statistical computations
clear
load_gaimc_graph('minnesota');
%%

%% solve the problem
[edges1,edges2,val]=find(triu(A));

%% Load R data corresponding to table 1
% (each file enables to reproduce one column of the table)
load minnesota/table_minesota_f1_sig1
%load minnesota/table_minesota_f1_sig2
%load minnesota/table_minesota_f1_sig3
%load minnesota/table_minesota_f2_sig1
%load minnesota/table_minesota_f2_sig2
%load minnesota/table_minesota_f2_sig3

y_val = y;
rep = 10;

tic
k = 2;
for j=1:rep
  y = y1(:,j);
  lambda = sort(abs(y));
  %lambda = linspace(7.5e-04,9.2e-04,2b0);
  for i=1:length(lambda)
    n=length(y);
    rho = 5*lambda(i);
    [ x ,history1] = gtf_admm_v2(y,[edges1,edges2],k,lambda(i),rho);
    mse_list(i,j) = norm(x-y_val)^2/n;
    snr_list(i,j) = 20*log10(norm(y_val)/norm(y_val-x));
  end
end
[aMSE,bMSE] = min(mse_list);
[cMSE,dMSE] = max(snr_list);
mean(aMSE)
mean(cMSE)
std(cMSE)
toc
