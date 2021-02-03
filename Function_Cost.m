function [J,greshapeee]=Function_Cost(qd,q)


    Q=eye(3,3);
    R=eye(2,2);
    i1=1;
    J=0;    
    
     for k=1:5:496 %% Gradiente de la funcion de costo
     ex=(qd(1:3)'-q(k:k+2,1)); %% Define el error [x_deseado y_deseado theta_deseado]-[x y theta]
     u=(q(k+3:k+4,1));
       
     %%%Gradiente de J es dJ/dx=x dJ/dy=y..pero el error = -e con su
     %%%respectiva ganancia
     J=J+1/2*(ex'*Q*ex+u'*R*u);%%DEFINE LA FUNCION DE COSTO
     
     g(i1,:)=[-ex(1)*Q(1,1) -ex(2)*Q(2,2) -ex(3)*Q(3,3) u(1)*R(1,1) u(2)*R(2,2)];%% EVALUA EL GRADIENTE
     
     i1=i1+1;
     end
     g=g';
     greshapeee=reshape(g,500,1);



end