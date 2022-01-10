% clear; close all
clear
%% Generate data

PI=pi;

colorsgt={[ .5 0.447058823529412 0.741176470588235 ],...
        [0.850980392156863   0.825490196078431   0.098039215686275],...
        [0.929411764705882   0.694117647058824   0.625490196078431],...
        [0.929411764705882   0.394117647058824   0.325490196078431],...
         };

rgb=@(x,y,z) [x, y,z]/255;
colors={ rgb(84, 177, 159),
         rgb(217, 108, 25),
         rgb(211, 101, 159),
         rgb(128, 63, 189),
       };

% rng(2)
N=2;
% limits = [ -fix(abs(10*rand)) fix(abs(10*rand));
%            -fix(abs(10*rand)) fix(abs(10*rand));
%            -fix(abs(10*rand)) fix(abs(10*rand));
%          ];
% for i=1:N
%     part{i}=[limits(i,1)  limits(i,1)+(limits(i,2)-limits(i,1))*rand() limits(i,2)];
% end
% part{:}

% for i=1:2^N
%     zone(:,:,i)
% end
% [ v{1:N} ] = ndgrid(part{:});
% cell2mat(cellfun(@(x) reshape(x,[],1),v,'UniformOutput',0))

% -5 2 4
% -8 0 7


par(:,:,1)=[-5 2;
            -8 0];
par(:,:,2)=[-5 2;
            0 7];
par(:,:,3)=[2 4;
            -8 0];
par(:,:,4)=[2 4;
            0 7];
x1 = linspace(min(par(1,:)),max(par(1,:))-0.1,20);
x2 = linspace(min(par(2,:)),max(par(2,:))-0.1,20);
v={x1,x2};


[X{1:N}] = ndgrid(v{:});
x=cell2mat(cellfun(@(x) reshape(x,[],1),X,'UniformOutput',0))';
% x_orig=[X{1}(:)';X{2}(:)'];
% sum(sum(x==x_orig,2)==size(x,2))==n

m(:,:,1)=[1;
          0];
m(:,:,2)=[5;
          3];
m(:,:,3)=[15;
          4];
m(:,:,4)=[2;
          2];

n(:,:,1)=[0];
n(:,:,2)=[0];
n(:,:,3)=[-40];
n(:,:,4)=[0];

theta=vertcat(m, n);

y = pwa(par,theta,x);

figure(1)
subplot(2,1,1)
cla
plot_pwa(par,x,y,colorsgt)
view(7,17)
title('Ground truth');
xlabel('x1')
ylabel('x2')

modes=size(theta,3);

%% Initialize estimated parameters

pi_hat=repmat(1/modes,1,modes);
maxErr=1e-4;
emMaxIter=200;

%% EM Algo
[C,d,responsabilities,pi,Sigma] = emgm_estimate(x,y,modes,emMaxIter,maxErr,m,n,colors);

sort(reshape(m,size(m,1),size(m,3)))
sort(reshape(C,size(C,1),size(C,3)))

return
%%

figure(1)
subplot(2,1,1)
cla
plot_pwa(par,x,y,colorsgt)
view(7,17)
title('Ground truth');
xlabel('x1')
ylabel('x2')


figure(1)
subplot(2,1,2)
cla
plot_responsibles(x, y, responsabilities, C, d, colors);
view(7,17)
title(['EM Gaussian Mixture iter=' num2str(emInd) ])
zlim([min(y) max(y)])
xlim([min(x(1,:)) max(x(1,:))])
ylim([min(x(2,:)) max(x(2,:))])
xlabel('x1')
ylabel('x2')
