clear;
clc;
close all;

%% question 1

x=[1/2 1/3]';
d=2;
y=task_1(x,d)

% y =
% 
%     1.0000
%     0.3333
%     0.1111
%     0.5000
%     0.1667
%     0.2500

%% Task 2

X=[1 0;0 1];
tau=0.1;
Xp=task_2(X,tau)

% Xp =
% 
%     0.9048         0
%    -0.0861    0.9048

%% Task 3

K = 10000;
d = 1;
tau = 1;
[A,Y,Yp]=task_3(K,d,tau);
A

% A =
% 
%     1.0000   -0.0000   -0.0000
%    -0.0770    0.3692   -0.0002
%     0.0000    0.0000    0.3679
% which is approximated to 
% A =
% 
%     1.0000         0         0
%    -0.0800    0.3700         0
%          0         0    0.3700

%% Task 4

K = 1000;
d = 10; 
tau= 0.001;
[A,Y,Yp]=task_3(K,d,tau);
norm(Yp-A*Y,"fro")/norm(Y',"fro")/tau
% linearizaztion error is a quite small value of 0.0062

[V, D] = eig(A');
exp_lambda=diag(D);
lambda=log(exp_lambda)/tau;% eigenvalues

figure()
scatter(real(lambda),imag(lambda));
title('The eigenvalues \lambda_n')
xlabel('Real');
ylabel('Img');

%eigenvalues in the range Re(λ) ≥ −5.5 lie at −5,−4, . . . , −1, 0

%% Task 5

c1=V(10);
c2=V(11);
norm(A'*c1-exp(lambda(1)*tau)*c1)
norm(A'*c2-exp(lambda(2)*tau)*c2)

% the norm of A'*c1-exp(lambda(1)*tau)*c1 and 
% norm of A'*c2-exp(lambda(2)*tau)*c2 are quite small with values of 
% 8.3319e-11 and 5.8565e-11 which could
% verify A'*c1=exp(lambda(1)*tau)*c1
% and A'*c2=exp(lambda(2)*tau)*c2

%% Task 6

x1=-4:0.001:4;
x2=x1.^2;

figure()
plot(zeros(size(x1)),x1);
hold on
plot(x1,x2);
axis([-4,4,-4,4]);
xlabel('x_1');
ylabel('x_2');
legend('Z1','Z2')
title('Zero-level sets Z1 and Z2')

% According to the plots of zero-level sets Z1 and Z2
% it can be observed that the zero-level sets Z1 and Z2 can only intersect
% at the origin (0,0)



%% Task 7

x01=[sqrt(2);2];
x02=[0;-2];
tspan=[0 20];
[t1,y1] = ode45(@(t,y) odefun1(t,y), tspan, x01);
[t2,y2] = ode45(@(t,y) odefun1(t,y), tspan, x02);

figure()
plot(zeros(size(x1)),x1,'c');
hold on
plot(x1,x2,'b');
hold on
plot(y1(:,1),y1(:,2),'ro')
hold on
plot(y2(:,1),y2(:,2),'m*')
axis([-4,4,-4,4]);
xlabel('x_1');
ylabel('x_2');
legend('Z1','Z2','trajectory from initial point x01','trajectory from initial point x02')
title('Zero-level sets and trajectories')

% It is verified that if initial conditions lie on the invariant sets
% the trajectories should not leave them and the trajectories all converge 
% towards the origin



%% function task_1

function y=task_1(x,d)
m=1;
 for i=1:(d+2)
   for j=1:(d+2-i)
        y(m,:)=x(1)^(i-1)*x(2)^(j-1);
        m=m+1;
    end
 end
end

%% function task_2

function Xp=task_2(X,tau)
m=size(X,2);
tspan=[0 tau];
for i=1:m
    y0=X(:,i);
    [t,y] = ode45(@(t,y) odefun1(t,y), tspan, y0);
    pp=size(y,1);
    Xp(:,i)=y(pp,:)';
end
end

function dydt = odefun1(t,y)
dydt = zeros(2,1);
dydt(1) = -y(1);
dydt(2) = -y(2)-y(1)^2;
end

%% function task_3

function [A,Y,Yp]=task_3(K,d,tau)
X=-1 + 2*rand(2,K);
% Y=task_1(X,d);
Xp=task_2(X,tau);
% Yp=task_1(Xp,d);
% A = lsqminnorm(Y,Yp);
Y=[];
for i=1:K
Y=[Y task_1(X(:,i),d)];
end
Yp=[];
for i=1:K
Yp=[Yp task_1(Xp(:,i),d)];
end
A=lsqminnorm(Y',Yp');
A=A';
end

% function dydt = odefun2(t,y)
% dydt = zeros(2,1);
% dydt(1) = y(1);
% dydt(2) = y(1)^2-y(2);
% end