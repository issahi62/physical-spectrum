%***************************
%% INITIALIZING 
%**************************
close all 
clc 
clear all

%************************
%% DASHBOARD 
%*************************
B = 0.1;
I = 0.5; 
Omega =1; 
T = 2/sqrt(2); 
TT = [0, 1.5, 3, 4.5];
TT0=0; 
%Amplitdue section
Yahoo =[];%empty array to store the amplitudes
Storer_TT_1 =[];
Storer_TT_2 =[];
Storer_TT_3 =[];
Storer_TT_4 =[];
%G = zeros(300,300); 
%G1 =zeros(300, 300); 

%************************
%% INDEXING 
%*************************
% Number of cells
t = linspace(-5, 15, 300);
C_Omega = linspace(-3, 3, 300); 


%************************
%% PHYSICAL SPECTRUM 
%*************************

for i = 1:length(t) 
% Integral section
    Yahoo =[(Yahoo); sqrt(pi/2)*I*B^2*T*exp(1/2*B^2*T^2)*exp(-2*B*t(i))];
    for j = 1: length(C_Omega)
   
    f = @(s) exp((-1/8*Omega.^2)*s.^2).*cos(C_Omega(j).*s).*erfc((s-(2.*t(i))...
        - B.*T^2)/sqrt(2).*T); 
     Q = integral(f, 0, inf);
      f1 = @(s) exp((-1/8*Omega.^2)*s.^2).*erfc((s-(2.*t(i))...
        - B.*T^2)/sqrt(2).*T); 
     Q1 = integral(f1, 0, inf);
    G(i,j) = Q; 
    G1(i,j) = Q1; 
    end
      
    Final = Yahoo.*G;
    Final2 = Yahoo.*G1;
end


%************************
%% AT Different time values 
%*************************

for i = 1: length(C_Omega)
    %T =0
    KK2= sqrt(pi/2)*I*B^2*T*exp(1/2*B^2*T^2)*exp(-2*B*TT(1));
    f2 = @(s) exp((-1/8*Omega.^2)*s.^2).*cos(C_Omega(j).*s).*erfc((s-(2.*TT(1))...
        - B.*T^2)/sqrt(2).*T);
    Q2 = integral(f, 0, inf);
    Mul_TT = KK2*Q2;
    Storer_TT_1 =[(Storer_TT_1), Mul_TT];
    % T = 1.5
    KK3= sqrt(pi/2)*I*B^2*T*exp(1/2*B^2*T^2)*exp(-2*B*TT(2));
    f3 = @(s) exp((-1/8*Omega.^2)*s.^2).*cos(C_Omega(j).*s).*erfc((s-(2.*TT(2))...
        - B.*T^2)/sqrt(2).*T);
    Q3 = integral(f, 0, inf);
    Mul_TT = KK3*Q3;
    Storer_TT_2 =[(Storer_TT_2), Mul_TT];
    %T = 3
    KK4= sqrt(pi/2)*I*B^2*T*exp(1/2*B^2*T^2)*exp(-2*B*TT(3));
    f4 = @(s) exp((-1/8*Omega.^2)*s.^2).*cos(C_Omega(j).*s).*erfc((s-(2.*TT(3))...
        - B.*T^2)/sqrt(2).*T);
    Q4 = integral(f, 0, inf);
    Mul_TT = KK4*Q4;
    Storer_TT_3 =[(Storer_TT_3), Mul_TT];
    % T = 4.5
    KK5= sqrt(pi/2)*I*B^2*T*exp(1/2*B^2*T^2)*exp(-2*B*TT(4));
    f5 = @(s) exp((-1/8*Omega.^2)*s.^2).*cos(C_Omega(j).*s).*erfc((s-(2.*TT(4))...
        - B.*T^2)/sqrt(2).*T);
    Q5 = integral(f, 0, inf);
    Mul_TT = KK5*Q5;
    Storer_TT_4 =[(Storer_TT_4), Mul_TT];
end
%************************
%% PLOTTING OF FIGURES 
%*************************
figure(1);
pcolor(t, C_Omega, Final'); 
xlabel('t'); ylabel('$\Delta \omega$', 'Interpreter', 'Latex');
title('Physical Spectrum')
shading interp
colormap('jet')
colorbar
figure(2);
pcolor(t, C_Omega, Final');
xlabel('t'); ylabel('$\Delta \omega$', 'Interpreter', 'Latex');
title('Physical Spectrum $\Delta \omega$=0', 'Interpreter', 'Latex')
shading interp
colorbar
colormap('jet')

%************************
%% PLOTTING OF FIGURES ID PLOT of TT values 
%*************************
figure(3);
subplot(2,2,1)
plot(C_Omega, Storer_TT_1); 
xlabel('t'); ylabel('$\Delta \omega$', 'Interpreter', 'Latex');
title('Physical Spectrum t=0', 'Interpreter', 'Latex')
subplot(2,2,2)
plot(C_Omega, Storer_TT_2); 
xlabel('t'); ylabel('$\Delta \omega$', 'Interpreter', 'Latex');
title('Physical Spectrum t=1.5', 'Interpreter', 'Latex')
subplot(2,2,3)
plot(C_Omega, Storer_TT_3); 
xlabel('t'); ylabel('$\Delta \omega$', 'Interpreter', 'Latex');
title('Physical Spectrum t=3', 'Interpreter', 'Latex')
subplot(2,2,4)
plot(C_Omega, Storer_TT_4); 
xlabel('t'); ylabel('$\Delta \omega$', 'Interpreter', 'Latex');
title('Physical Spectrum t=4.5', 'Interpreter', 'Latex')



