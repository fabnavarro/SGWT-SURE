clear
load('res_homogene.mat')
load('homob3.mat')
mseb3 = mse;
clear mse
load('homob4.mat')
mseb4 = mse;
clear mse
load('homob5.mat')
mseb5 = mse;
clear mse

fontx  = 15;
h2=figure;
plot(-SnR_list,10*log(mse_list')/log(10)','ro-',...
    -SnR_list,10*log(mse_list1')/log(10)','r*-',...
    -SnR_list,10*log(mse_list2')/log(10)','r^-',...
    -SnR_list,10*log(mse_listL')/log(10)','kx--',...%'linewidth',2,'markersize',12);% -SnR_list,10*log(mse_b2')/log(10)','b*-',...
    -SnR_list,10*log(mseb3')/log(10)','bs-',...
    -SnR_list,10*log(mseb4')/log(10)','bd-',...
    -SnR_list,10*log(mseb5')/log(10)','bp-',...
    'linewidth',2,'markersize',12);

xlabel('Noise Level: -SnR in (dB)','fontsize',15)
ylabel('Denoising Performance in dB','fontsize',15)
ylabel('Denoised -SnR in dB','fontsize',15)
g=legend('$k=0$','$k=1$',...
  '$k=2$','Laplacian',...%'Level dependant \beta=2 (J=7)',...
  'LD $b=3$', 'LD $b=4$', 'LD $b=5$', 'LD $b=6$');
set(g,'interpreter', 'latex','fontsize',15,'location','northwest')
grid on;
ax = gca; % current axes
ax.FontSize = fontx;
axis equal 
ylim([-40 0])
xlim([-30 10])
%axis tight

