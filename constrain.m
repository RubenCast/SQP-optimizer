function [c,ceq]=constrain(x,X0,Ts)

c=[];
i=1;
q(1:3)=X0;
for k=1:5:496
if k==1
  
ceq(i,1)=X0(1)-x(k,1)+Ts*x(k+3,1)*cos(x(k+2,1));
ceq(i+1,1)=X0(2)-x(k+1,1)+Ts*x(k+3,1)*sin(x(k+2,1));
ceq(i+2,1)=X0(3)-x(k+2,1)+Ts*x(k+4,1);
i=i+3;
else

ceq(i,1)=x(k-5,1)-x(k,1)+Ts*x(k+3,1)*cos(x(k+2,1));
ceq(i+1,1)=x(k-4,1)-x(k+1,1)+Ts*x(k+3,1)*sin(x(k+2,1));
ceq(i+2,1)=x(k-3,1)-x(k+2,1)+Ts*x(k+4,1);
i=i+3;

end

end
end