function plot_pwa(par,x,y,colors)
  hold on
  for sidx=1:length(par)
    mask=sum(x>=par(:,1,sidx)&x<par(:,2,sidx))==size(x,1);
    scatter3(x(1,mask),x(2,mask),y(mask),10,colors{sidx})
  end
  hold off
end
