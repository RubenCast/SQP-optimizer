%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Rubén Castro González      %%%
%%%                               %%%
%%% SECUENTIAL QUADRATIC PROGRAM  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%   qd=[4,9,pi]
%   q0=[0 0 0]
%
%
%
%
%

clc
clear all;
etha=.5;              %% Gain of the directional derivative
tao=.5;               %% Bisection scalar
qd=[4,9,pi];           %% Desired point
q0=[0;0;0];           %% Initial position
Ts=0.1;               %% Increment of time
B=eye(500,500);       %% Initial Hessian of the Lagrangian
tryk=1;               %% Tries
maxtry=200;           %% Maximmum of tries
rho=.5;               %% Derivative Directional Parametter
mu=0;                 %% Penalty parameter 
Tol=0.0001;           %% Tolerance in the diference of the two last function cost values
 error=1;             %% error beetwin the last function cost values
for i=1:500             
   q(i,1)=rand(1);    %%Initial values of trayectories  
end


fprintf(1, ' try   Costfunction   ||constraints||1    Merit    Merit_ap     alpha      Cost_error \n');
fprintf(1, '-------------------------------------------------------------------------------------- \n');



%%%% Evaluating: Cost function(fi), Gradient of the Cost
%%%% Function(deltafi),constraints (ceq), Gradient of the constraints(A)
     [A]=build_A(q,Ts);
     [c,ceq]=constrain(q,q0,Ts); 
     [fi,deltafi]=Function_Cost(qd,q);
   
%%%%

 
 errorminimo=.0001;

 a=1;

 while error>=Tol && tryk<maxtry
  %%% KKT SYSTEM
 Left_side=  [B             -A';
              A        zeros(300,300)];
 Right_side=[-deltafi;
              -ceq];
 PK=inv(Left_side)*Right_side; %%Solving for the step
 
 Lambda=PK(501:800);
 
 Lagrangian_gradient=deltafi-A'*Lambda;
 %%%%% Penalty parameter mu
 
 if PK(1:500,1)'*B*PK(1:500,1)>0
     sigma=1;
 else
     sigma=0;
 end
mutry=tryk+1;
 if mu>=(deltafi'*PK(1:500,1)+(sigma/2)*(PK(1:500,1)'*B*PK(1:500,1)))/((1-rho)*norm(ceq,1))
  mu= mu;
 else
     mu=(deltafi'*PK(1:500,1)+(sigma/2)*(PK(1:500,1)'*B*PK(1:500,1)))/((1-rho)*norm(ceq,1))+.1;
 end
 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 alphak=1;
%  Di2=-PK(1:500,1)'*B*PK(1:500,1)+PK(1:500,1)'*A*Lambda-mu*norm(ceq,1);
  Di=deltafi'*PK(1:500,1)-mu*norm(ceq,1);%-mu*norm(ceq,1)%%Directional derivative
 %%merit function
 fi_mux=fi + mu*norm(ceq,1); %%with qk-1
 
 %%//////with q=q+alpha*p
 q1=q+alphak*PK(1:500,1);
 
     [fi1,deltafi1]=Function_Cost(qd,q1);
     [c1,ceq1]=constrain(q1,q0,Ts);
     
     if mu>=(deltafi1'*PK(1:500,1)+(sigma/2)*(PK(1:500,1)'*B*PK(1:500,1)))/((1-rho)*norm(ceq1,1))
     mu1= mu;
     else
     mu1=(deltafi1'*PK(1:500,1)+(sigma/2)*(PK(1:500,1)'*B*PK(1:500,1)))/((1-rho)*norm(ceq1,1))+.1;
     end
     fi_alphax=fi1 + mu1*norm(ceq1,1);
     
     
 %%%%%%%%%%%%Bisection%%%%%%%%%%%%%%%%%%%%%%%%
 maxtrybisection=1;
 
 mu2=0;
 while fi_alphax>fi_mux+etha*alphak*Di && maxtrybisection<=30
 alphak=tao*alphak;
 
 q2=q+alphak*PK(1:500,1);
 
     
     [fi2,deltafi2]=Function_Cost(qd,q2);
     [c2,ceq2]=constrain(q2,q0,Ts);
     
     if mu>=(deltafi2'*PK(1:500,1)+(sigma/2)*(PK(1:500,1)'*B*PK(1:500,1)))/((1-rho)*norm(ceq2,1))
     mu2= mu;
     else
     mu2=(deltafi2'*PK(1:500,1)+(sigma/2)*(PK(1:500,1)'*B*PK(1:500,1)))/((1-rho)*norm(ceq2,1))+.1;
     end
     
%    mu1=norm(Lambda,Inf)+delta;
     fi_alphax=fi2 + mu2*norm(ceq2,1);
 
     
 maxtrybisection=maxtrybisection+1;
 end
  
  
  %%%Actualiza q con alpha
 q3=q+alphak*PK(1:500,1);
 
  sk=alphak*PK(1:500,1);
%  sk=q3-q;  %alphak*PK(1:500,1);
 q=q3;
 
 [A]=build_A(q,Ts);
 [c,ceq]=constrain(q,q0,Ts);
 [fi,deltafi]=Function_Cost(qd,q);

  Lagrangian_gradientk=deltafi-A'*Lambda;
  yk=Lagrangian_gradientk-Lagrangian_gradient;
  
  
 %%%Updating Hessian of the Lagrangian 
 theta=0.1;
 if sk'*yk>=theta*B*sk %%curvature condition
B=B;
else
    if sk'*yk >= 0.2*sk'*B*sk 
        theta=1;
    else
        theta=(0.8*sk'*B*sk)/(sk'*B*sk-sk'*yk);
    end

rk=theta*yk+(1-theta)*B*sk;

if sk'*rk>0
    B=B-((B*sk*sk'*B)*(1/(sk'*B*sk)))+(rk*rk')*(1/(sk'*rk));
   
end
 end
 mu=mu2;

 tryk=tryk+1;
 ficalc(tryk,1)=fi;
 error=abs(ficalc(tryk-1)-ficalc(tryk));

 fprintf(1,' %3i   %14.8e   %8.2e        %8.2e  %8.2e  %2i   %8.2e\n', ... 
               tryk-1, fi, norm(ceq,1), fi_mux, fi_alphax, alphak,   error);
 
 end
  PLOTS_Position(q,2);
 
 
 
 
 