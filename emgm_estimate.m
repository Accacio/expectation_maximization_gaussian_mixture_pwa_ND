function [C,d,responsabilities,pi,Sigma] = emgm_estimate(x,y,modes,emMaxIter,maxErr,m,n,colors)
% EMGM_ESTIMATE -

    pi=repmat(1/modes,1,modes);
    OldclusterSize=zeros(1,modes);

    % TODO(accacio): find way to initialize it better (use gradient)
    % Something like:
    % indexpts1=randi([1 N-1],1,3);
    % C=(y(:,indexpts1+1)-y(:,indexpts1))./(x(:,indexpts1+1)-x(:,indexpts1));
    % d=y(:,indexpts1)-C.*x(:,indexpts1);
    % C=C+rand(1,M);
    % d=d+rand(1,M);

    % C=5*rand(size(x,1),1,modes);
    % d=5*rand(1,1,modes);
    C=m+2*rand(size(m));
    d=n+2*rand(size(n));
    % C=10*rand(size(m));
    % d=10*rand(size(n));

    eps=10;
    Sigma(:,:,1:modes)=repmat(eps*eye(1),1,1,modes);

    for emInd=1:emMaxIter

        responsabilities=calculate_responsabilities(x,y,C,d,Sigma,pi);


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


        [C, d, pi] = update_parameters(x, y, responsabilities);
        [~,z_hat]=max(responsabilities,[],1);
        for i=1:modes
            z_i=find(z_hat==i);
            clusterSize(i)=size(z_i,2);
            if OldclusterSize(i)==clusterSize(i)
                Sigma(:,:,i)=Sigma(:,:,i)*.9;
            else
                Sigma(:,:,i)=Sigma(:,:,i)*1.;
            end
        end
        if sum(Sigma < maxErr)==modes
            break;
        end
        OldclusterSize=clusterSize;
    end

end
