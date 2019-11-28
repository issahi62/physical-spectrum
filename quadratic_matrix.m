%************************
%% INITIALIZE 
%************************
clc 
close all 
clear 


%************************
%% DASHBOARD
%************************
% Idea if to transform a plot into a quadratic form
xa = -100:100; 

%************************
%% MATRIX
%************************
A = [ 2, 2; 2 2]; % choose any matrix
NA = zeros(length(xa)); 
%************************
%% OPERATION
%************************
for xi = 1:length(xa)
    for xj= 1:length(xa)
        W = [xa(xi) xa(xj)]'; 
        Normalizer = W'*W; 
        NA(xi, xj) = (W'*A*W)/Normalizer; 
    end 
end
 
%************************
%% PLOTTING
%************************
figure; clf;
subplot(1,2,1)
surf(xa,xa, NA'); 
shading interp
xlabel('x_axis'); 
ylabel('y_axis');
colormap('jet');
rotate3d on; 
axis square; 

subplot(1,2,2); 
plot(A); 
axis square; 