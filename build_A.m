%18/07/2018
%Funcion que dado un vector de estados y controles devuelve la matriz A

function [A]= build_A(q,Ts)


A=zeros(500,300);
% gk=[1 0 ts*sin(theta)*u_1 -ts*cos(theta) 0;
%     0 1 ts*cos(theta)*u_1 -ts*sin(theta) 0;
%     0 0       1                0         0];

I_1=[1 0 0;
     0 1 0;
     0 0 1;
     0 0 0;
     0 0 0];
 
 i=1;
     for k=1:5:495
     A(k:k+4,i:i+2) = [-1                          0                         0;
                        0                          -1                        0;
                     -Ts*sin(q(k+2))*(q(k+3))      Ts*cos(q(k+2))*(q(k+3))  -1;
                      Ts*cos(q(k+2))               Ts*sin(q(k+2))            0;
                       0                            0                        Ts];
                       
                       
                 
     A(k:k+4,i+3:i+5) = I_1;
     
     i=i+3; 
     end
     
     A(496:500,298:300) = [-1                          0                                0;
                            0                          -1                               0;
                          -Ts*sin(q(496+2))*q(496+3)   Ts*cos(q(496+2))*(q(496+3))     -1;
                           Ts*cos(q(496+2))            Ts*sin(q(496+2))                 0;
                            0                          0                               Ts];
  A=A';                      
end