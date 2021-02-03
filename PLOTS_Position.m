function PLOTS_Position(q,m)
i=1;
for k=1:5:496
    x(i,1)=q(k);
    y(i,1)=q(k+1);
    theta(i,1)=q(i+2);
    i=i+1;
    
end

figure(m)
 plot(x(:,1),y(:,1),'->');



end