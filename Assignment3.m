clear;
clc;
close all;

addpath('C:\Users\linxi\Documents\MATLAB\chebfun-master\chebfun-master');
addpath('D:\eigtool-master\eigtool-master');

%% Task1

d = [-25 25]; 
u = @(z)12*sech(z)*sech(z);
u = chebfun(u,d);
N1 = chebop(@(z) diff(z,2)+u*z, [-25 25], 'periodic');
[V, D] = eigs(N1,'LR',3);
lambda1=diag(D);
figure()
scatter(real(lambda1),imag(lambda1));
title('The eigenvalues \lambda_n')
xlabel('Real');
ylabel('Img');

%% Task2

phi1=abs(V{1});
lambda_1=lambda1(1);
lambda1_approx = @(z) exp(sqrt(lambda_1)*z)*phi1(z);
g1 =chebfun(lambda1_approx,[-25 25], 'periodic');
figure()
plot(g1)
title('The function values in the domain [-25,25]')
xlabel('Z');
% The estimate of \delta_1(0) is the right boundary of
% exp(sqrt(\lambda_1)z)\Phi^1(z), which is 25

%% Task3


dom = [-25 25];
x = chebfun('x',dom);
tmax = 0.5;
S = spinop(dom,[0 tmax]);
S.lin = @(u) - diff(u,3);
S.nonlin = @(u) -3*diff(u.^2);
u = @(z)12*sech(z)*sech(z);
S.init= chebfun(u,dom);

figure()
N = 800;   % numer of grid points
dt = 5e-6; % time-step
tic, u = spin(S,N,dt,'plot','off'); 
time_in_seconds = toc;
plot(S.init), hold on, plot(u), hold off
text(4.4,1300,'t = 0'), text(13.5,1300,'t = 0.5')
title('u(0.5,z) and u(0,z)')
xlabel('Z');
legend('u(0,z)','u(0.5,z)');

%% Task4

% Use the Tasks 1 and Task 2 to u(0.5, z)

N2 = chebop(@(z) diff(z,2)+u*z, [-25 25], 'periodic');
[V, D] = eigs(N2,'LR',3);
lambda2=diag(D);
figure()
scatter(real(lambda2),imag(lambda2));

% check the eigenvalues keep the same
t=0.5;
phi1=abs(V{1});
lambda_1=lambda2(1);
lambda1_approx = @(z) exp(sqrt(lambda_1)*z)*phi1(z);
g2 =chebfun(lambda1_approx,[-25 25], 'periodic');
figure()
plot(g2)
% delta0=25;
hold on
plot(exp(4*sqrt(lambda_1)*lambda_1*t)*g1)
title('u(0.5,z) and \delta_1(0.5)')
xlabel('Z');
legend('u(0.5,z)','\delta_1(0.5)');

%% Task5

% Verify the reconstruction formula (4.1) using the eigenvalues and -functions
% from Task 1

dom = [-25 25];
x = chebfun('x',dom);
tmax = 0.5;
S = spinop(dom,[0 tmax]);
S.lin = @(u) - diff(u,3);
S.nonlin = @(u) -3*diff(u.^2);
u = @(z)12*sech(z)*sech(z);
S.init= chebfun(u,dom);

% calculate u(0.5,z) by reconstruction from (4.1)

[V, D] = eigs(N2,'LR',3);%N2 corresponds to u(0.5, z)
lambda=diag(D);
phi1=V{1};
phi2=V{2};
phi3=V{3};
u05z=@(z) 4*sqrt(lambda(1))*phi1(z)^2+4*sqrt(lambda(2))*phi2(z)^2+4*sqrt(lambda(3))*phi3(z)^2;
g3 =chebfun(u05z,[-25 25], 'periodic');

% calculate u(0.5,z) by reconstruction from (4.1)

[V, D] = eigs(N1,'LR',3);
lambda=diag(D);
phi1=V{1};
phi2=V{2};
phi3=V{3};
u0z=@(z) 4*sqrt(lambda(1))*phi1(z)^2+4*sqrt(lambda(2))*phi2(z)^2+4*sqrt(lambda(3))*phi3(z)^2;
g0 =chebfun(u0z,[-25 25], 'periodic');


%########################In the case of u(0,z)#########################%
figure()
plot(S.init)
hold on
plot(g0)
title('u(0,z)')
xlabel('Z');
legend('using the formula from Task 1','reconstruction from (4.1)');
% The plot u(0.5, z) directly using the formula from Task 1 overlap with 
% the reconstruction from (4.1).

%########################In the case of u(0.5,z)#########################%
figure()
u = spin(S,N,dt,'plot','off'); 
plot(u)
hold on
plot(g3)
title('u(0.5,z)')
xlabel('Z');
legend('using the formula from Task 1','reconstruction from (4.1)');
% The plot u(0.5, z) directly using the formula from Task 1 overlap with 
% the reconstruction from (4.1).

%% Task6

z1=10;
la=1;
del1=sqrt(2);
dom2=[-z1,z1];

% solve differential equation 4.2
N=chebop(-10,10);
N.op=@(t,u) diff(u,2)+4*u^3-u;
N.rbc=[del1*exp(-10);-del1*exp(-10)];
M=N\0;

% u(t0,z)
u = @(z)2*sech(z)*sech(z);
compare= chebfun(u,dom2);

figure()
plot(compare)
hold on
plot(4*sqrt(la)*del1*exp(-sqrt(la)*z1)^2) %using the given solutions of 4.2
hold on
plot(4*M^2)
title('\delta=sqrt(2)')
xlabel('Z');
legend('u(t0,z)','using the given solutions of 4.2','solve differential equation 4.2');

%% Task7

% '\delta=sqrt(2)/2'
del2=sqrt(2)/2;
% solve differential equation 4.2
N=chebop(-10,10);
N.op=@(t,u) diff(u,2)+4*u^3-u;
N.rbc=[del2*exp(-10);-del2*exp(-10)];
M=N\0;

figure()
plot(compare)
hold on
plot(4*sqrt(la)*del2*exp(-sqrt(la)*z1)^2) %using the given solutions of 4.2
hold on
plot(4*M^2)
title('\delta=sqrt(2)/2')
xlabel('Z');
legend('u(t0,z)','using the given solutions of 4.2','solve differential equation 4.2');



%'\delta=sqrt(2)*2'
del3=sqrt(2)*2;
% solve differential equation 4.2
N=chebop(-10,10);
N.op=@(t,u) diff(u,2)+4*u^3-u;
N.rbc=[del3*exp(-10);-del3*exp(-10)];
M=N\0;

figure()
plot(compare)
hold on
plot(4*sqrt(la)*del3*exp(-sqrt(la)*z1)^2) %using the given solutions of 4.2
hold on
plot(4*M^2)
title('\delta=sqrt(2)*2')
xlabel('Z');
legend('u(t0,z)','using the given solutions of 4.2','solve differential equation 4.2');


%###############################analysis##################################%
% In the three cases of different delta using the given solutions of 4.2 and 
% solve differential equation 4.2 are the same
% When delat=sqrt(2), the solutions get from 4.2 overlap with u(t0,z)
% When delat=sqrt(2)/2, the solutions get from 4.2 diverge from u(t0,z)
% to the left at the top, the top point is (0.693,2)
% When delat=sqrt(2)*2, the solutions get from 4.2 diverge from u(t0,z)
% to the right at the top, the top point is (-0.693,2)
% delat=sqrt(2)/2 and sqrt(2)*2 were symmetrical and have the same amount
% of deviation from the top point (0,2)
