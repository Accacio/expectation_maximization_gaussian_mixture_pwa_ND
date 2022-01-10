% clear; close all
clear
%% Generate data

PI=pi;
emMaxIter=200;

colorsgt={[ .5 0.447058823529412 0.741176470588235 ],...
        [0.850980392156863   0.825490196078431   0.098039215686275],...
        [0.929411764705882   0.694117647058824   0.625490196078431],...
        [0.929411764705882   0.394117647058824   0.325490196078431],...
         };

colors={[ .5 0.247058823529412 0.741176470588235 ], ...
        [0.850980392156863   0.425490196078431   0.098039215686275], ...
        [0.829411764705882   0.394117647058824   0.625490196078431], ...
        [0.329411764705882   0.694117647058824   0.625490196078431], ...
        };

par(:,:,1)=[-5 2;
            -8 0];
par(:,:,2)=[-5 2;
            0 7];
par(:,:,3)=[2 4;
            -8 0];
par(:,:,4)=[2 4;
            0 7];
n=2;
x1 = linspace(min(par(1,:)),max(par(1,:))-0.1,10);
x2 = linspace(min(par(2,:)),max(par(2,:))-0.1,10);
v={x1,x2};


[X{1:n}] = ndgrid(v{:});
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

% figure(1)
% subplot(2,1,1)
% cla
% plot_pwa(par,x,y,colorsgt)
% view(7,17)
% title('Ground truth');
% xlabel('x1')
% ylabel('x2')

M=size(theta,3);

pi=repmat(1/M,1,M);
eps=10;
Sigma(:,:,1:M)=repmat(eps*eye(1),1,1,M);

%% Initialize estimated parameters

pi_hat=repmat(1/M,1,M);

% C=m;
% d=n;
C=m+2*rand(size(m));
d=n+2*rand(size(n));
% C=10*rand(size(m));
% d=10*rand(size(n));

maxErr=1e-4;

%% EM Algo
OldclusterSize=zeros(1,M);
for emInd=1:emMaxIter

    responsabilities=calculate_responsabilities(x,y,C,d,Sigma,pi_hat);


    % figure(1)
    % subplot(2,1,2)
    % cla
    % plot_responsibles(x, y, responsabilities, C, d, colors);
    % view(7,17)
    % title(['EM Gaussian Mixture iter=' num2str(emInd) ])
    % zlim([min(y) max(y)])
    % xlim([min(x(1,:)) max(x(1,:))])
    % ylim([min(x(2,:)) max(x(2,:))])
    % xlabel('x1')
    % ylabel('x2')


    [C, d, pi_hat] = update_parameters(x, y, responsabilities);
    [~,z_hat]=max(responsabilities,[],1);
    for i=1:M
        z_i=find(z_hat==i);
        clusterSize(i)=size(z_i,2);
        if OldclusterSize(i)==clusterSize(i)
            Sigma(:,:,i)=Sigma(:,:,i)*.9;
        else
            Sigma(:,:,i)=Sigma(:,:,i)*1.;
        end
    end
    if sum(Sigma < maxErr)==M
        break;
    end
    OldclusterSize=clusterSize;
end

reshape(m,size(m,1),size(m,3))
reshape(C,size(C,1),size(C,3))

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
