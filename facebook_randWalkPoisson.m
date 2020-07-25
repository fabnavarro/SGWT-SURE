close all
clear
load('facebook/res_inhomogene.mat')
load('facebook/inhomoRandWalkR.mat')

h1 = figure;
plot(-SnR_list,10*log(mse_list')/log(10)','ro-',...
    -SnR_list,10*log(mse_list1')/log(10)','r*-',...
    -SnR_list,10*log(mse_list2')/log(10)','r^-',...
    -SnR_list,10*log(mse_listL')/log(10)','kx--',...%'linewidth',2,'markersize',12);% -SnR_list,10*log(mse_b2')/log(10)','b*-',...
    -SnR_list,10*log(mse')/log(10)','bs-','linewidth',2,'markersize',12);

xlabel('Noise Level: -SnR in (dB)','fontsize',25)
ylabel('Denoising Performance in dB','fontsize',25)
ylabel('Denoised -SnR in dB','fontsize',25)
lg=legend('Tf k=0','Tf k=1',...
  'Tf k=2','Laplacian',...%'Level dependant \beta=2 (J=7)',...
  'LD \beta=2');
set(lg,'fontsize',25,'location','northwest')
grid on;
ax = gca; % current axes
ax.FontSize = 20;
%ax.TickDir = 'out';
%ax.TickLength = [0.02 0.02];
%axis equal 
ylim([-45 2])
xlim([-30 10])
%axis tight
%saveTightFigure(h1,'facebook_inhomo.pdf')


clear
load('facebook/res_homogene.mat')
load('facebook/homoRandWalkR.mat')

h2=figure;
plot(-SnR_list,10*log(mse_list')/log(10)','ro-',...
    -SnR_list,10*log(mse_list1')/log(10)','r*-',...
    -SnR_list,10*log(mse_list2')/log(10)','r^-',...
    -SnR_list,10*log(mse_listL')/log(10)','kx--',...%'linewidth',2,'markersize',12);% -SnR_list,10*log(mse_b2')/log(10)','b*-',...
    -SnR_list,10*log(mse')/log(10)','bs-','linewidth',2,'markersize',12);

xlabel('Noise Level: -SnR in (dB)','fontsize',25)
ylabel('Denoising Performance in dB','fontsize',25)
ylabel('Denoised -SnR in dB','fontsize',25)
l=legend('Tf k=0','Tf k=1',...
  'Tf k=2','Laplacian',...%'Level dependant \beta=2 (J=7)',...
  'LD \beta=2');
set(l,'fontsize',25,'location','northwest')
grid on;
ax = gca; % current axes
ax.FontSize = 20;
%ax.TickDir = 'out';
%ax.TickLength = [0.02 0.02];
%axis equal 
ylim([-45 2])
xlim([-30 10])
%axis tight
%saveTightFigure(h2,'facebook_homo.pdf')

clear
load('facebook/res_poisson_dense.mat')
load('facebook/densePoissonR.mat')

h3=figure;
plot(-SnR_list,10*log(mse_list')/log(10)','ro-',...
    -SnR_list,10*log(mse_list1')/log(10)','r*-',...
    -SnR_list,10*log(mse_list2')/log(10)','r^-',...
    -SnR_list,10*log(mse_listL')/log(10)','kx--',...%'linewidth',2,'markersize',12);% -SnR_list,10*log(mse_b2')/log(10)','b*-',...
    -SnR_list,10*log(mse')/log(10)','bs-','linewidth',2,'markersize',12);

xlabel('Noise Level: -SnR in (dB)','fontsize',25)
ylabel('Denoising Performance in dB','fontsize',25)
ylabel('Denoised -SnR in dB','fontsize',25)
g=legend('Tf k=0','Tf k=1',...
  'Tf k=2','Laplacian',...%'Level dependant \beta=2 (J=7)',...
  'LD \beta=2');
set(g,'fontsize',25,'location','northwest')
grid on;
ax = gca; % current axes
ax.FontSize = 20;
%ax.TickDir = 'out';
%ax.TickLength = [0.02 0.02];
%axis equal 
ylim([-45 2])
xlim([-30 10])
%axis tight
%saveTightFigure(h3,'facebook_dense.pdf')

clear
load('facebook/res_poisson_sparse.mat')
load('facebook/sparsePoissonR.mat')

h4 = figure;
plot(-SnR_list,10*log(mse_list')/log(10)','ro-',...
    -SnR_list,10*log(mse_list1')/log(10)','r*-',...
    -SnR_list,10*log(mse_list2')/log(10)','r^-',...
    -SnR_list,10*log(mse_listL')/log(10)','kx--',...%'linewidth',2,'markersize',12);% -SnR_list,10*log(mse_b2')/log(10)','b*-',...
    -SnR_list,10*log(mse')/log(10)','bs-','linewidth',2,'markersize',12);

xlabel('Noise Level: -SnR in (dB)','fontsize',25)
ylabel('Denoising Performance in dB','fontsize',25)
ylabel('Denoised -SnR in dB','fontsize',25)
lg=legend('Tf k=0','Tf k=1',...
  'Tf k=2','Laplacian',...%'Level dependant \beta=2 (J=7)',...
  'LD \beta=2');
set(lg,'fontsize',25,'location','northwest')
grid on;
ax = gca; % current axes
ax.FontSize = 20;
%ax.TickDir = 'out';
%ax.TickLength = [0.02 0.02];
%axis equal 
xlim([-30 10])
ylim([-45 2])
%axis tight
%saveTightFigure(h4,'facebook_sparse.pdf')

tilefigs
