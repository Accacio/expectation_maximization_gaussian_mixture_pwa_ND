function y = pwa(par,theta,p)
    y=zeros(size(theta,2),size(p,2));
    for i=1:size(par,3)
        y=y+(sum((p>=par(:,1,i)&p<par(:,2,i)))==size(p,1)).*(theta(:,:,i)'*[p ;ones(1,size(p,2))]);
    end
end
