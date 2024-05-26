function [result] = AniVortexFFTinterpLambdasFun(Lambda1, Lambda2, z)
%% useful parameters
% flux quanta in Wb not uT-um^2
Phi_0 = 2068;
% scan height in um
%z = 4.8;
% pickup loop effective radius in um
b = 3.57;
% Peral length in two directions  
%Lambda1 = 400;
%Lambda2 = 400;
%% defind the real and k space
% defind real space range in m
dx = .1;
dy = .1;
xRange = 400 + dx/2;
yRange = 400 + dy/2;
x = -xRange:dx:xRange;
y = -yRange:dy:yRange;
[X,Y] = meshgrid(x,y);

% defind the k-space
% dkx = 2pi/(range of x)
dkx = pi/xRange;
dky = pi/yRange;
% (range of kx) = 2pi/dx
kxRange = pi/dx;
kyRange = pi/dy;
%dkx = 2*pi/(length(x) * dx);
%dky = 2*pi/(length(y) * dy);
kx = linspace(-kxRange, kxRange-dkx, length(x));
ky = linspace(-kyRange, kyRange-dky, length(y));
[Kx, Ky] = meshgrid(kx,ky);
%% compute the expression in Fourirer space
K = sqrt(Kx.^2 + Ky.^2);
hzk = Phi_0 * exp(-K*z) ./ (1 + Kx.^2 ./K * Lambda2 + Ky.^2 ./K * Lambda1);
%hzk = ifftshift(hzk);

figure(135)
imagesc(kx,ky,hzk)
colorbar;
shading flat
axis image
xlabel('k_x (1/\mum)')
ylabel('k_y (1/\mum)')

hz = ifft2(fftshift(hzk));
hz = fftshift(hz)/(dx*dy);
total_flux = sum(sum(real(hz))) * dx * dy/Phi_0;
figure(136)
imagesc(x,y,abs(hz))
shading flat;
axis image;
colorbar;
xlabel('x(\mum)')
ylabel('y(\mum)')
title_string = sprintf('Total flux in image = % 4.2g Phi_0', total_flux);
title(title_string)

figure(235)
imagesc(x,y,real(hz))
shading flat;
axis image;
colorbar;
xlabel('x(\mum)')
ylabel('y(\mum)')
title('real part')

figure(236)
imagesc(x,y,imag(hz))
shading flat;
axis image;
colorbar;
xlabel('x(\mum)')
ylabel('y(\mum)')
title('imag part')

%% zoom in to center part
dxp = .1;
dyp = .1;
% plot range
xP = 30;
yP = 30;

% number of points to be plotted
Nx = round(xP/xRange * length(x)/2);
Ny = round(yP/yRange * length(x)/2);
% half of every vector
halfX = length(x)/2;
halfY = length(y)/2;
% range to interpolate to
%xq = x(halfX - Nx):dxp:x(halfX + Nx + 2);
%yq = y(halfX - Ny):dyp:y(halfY + Ny + 2);
%[Xq,Yq] = meshgrid(xq,yq);
%kx_useful = kx(halfX - Nx: halfX + Nx + 2);
%ky_useful = ky(halfY - Ny: halfY + Ny + 2);

x_useful = x(halfX - Nx:halfX + Nx + 2);
y_useful = y(halfY - Ny:halfY + Ny + 2);

%[X_useful, Y_useful] = meshgrid(x_useful, y_useful);
%xq = xq - x(halfX + 1) * ones(1, length(xq));
%yq = yq - y(halfY + 1) * ones(1, length(yq));

real_hz = real(hz);
real_hz_useful = real_hz(halfX - Nx: halfX + Nx + 2, halfY - Ny: halfY+ Ny + 2);

%real_hz_useful_interp = interp2(X_useful,Y_useful,real_hz_useful,Xq,Yq);
figure(346)
imagesc(x_useful, y_useful, real_hz_useful)
colormap jet
%colorbar 
c = colorbar;
c.Label.String = '\muT';
set(gca, 'FontSize',20)
xlabel('x (\mum)')
ylabel('y (\mum)')
set(gca, 'ydir', 'normal')
set(gcf,'Position',[100 100 500 427])

%% convolve with pickup loop
dPU = .1;
PU_range = 4;
xPU = -PU_range:dPU:PU_range ;
yPU = -PU_range:dPU:PU_range;
[XPU, YPU] = meshgrid(xPU, yPU);
RPU = sqrt((XPU).^2 + (YPU).^2);
PU = 0 * XPU;
PU(RPU - b <.001)=1;
figure(456)
imagesc(xPU,yPU,PU)
shading flat;
axis image;
colorbar;
colormap bone;
xlabel('x (\mum)')
ylabel('y (\mum)')
set(gcf,'Position',[100 100 500 427])
set(gca, 'FontSize',20)

flux = conv2(real_hz_useful,PU,'same')* dx * dy* 1e3/Phi_0;
figure(138)
pcolor(x_useful,y_useful,flux)
shading flat;
axis image;
colorbar;
xlabel('x(\mum)')
ylabel('y(\mum)')
title('m\Phi_0')

%% useful range for vortex fitting
% plot range
xP2 = 10;
yP2 = 10;

% number of points to be plotted
Nx = round(xP2/xP * length(x_useful)/2);
Ny = round(yP2/yP * length(y_useful)/2);
% half of every vector
halfX = round(length(x_useful)/2)-1;
halfY = round(length(y_useful)/2)-1;
% range to interpolate to
%xq = x(halfX - Nx):dxp:x(halfX + Nx + 2);
%yq = y(halfX - Ny):dyp:y(halfY + Ny + 2);
%[Xq,Yq] = meshgrid(xq,yq);
%kx_useful = kx(halfX - Nx: halfX + Nx + 2);
%ky_useful = ky(halfY - Ny: halfY + Ny + 2);

x_useful2 = x_useful(halfX - Nx + 1:halfX + Nx);
y_useful2 = y_useful(halfY - Ny + 1:halfY + Ny);

%[X_useful, Y_useful] = meshgrid(x_useful, y_useful);
%xq = xq - x(halfX + 1) * ones(1, length(xq));
%yq = yq - y(halfY + 1) * ones(1, length(yq));

flux_useful = flux(halfX - Nx + 1: halfX + Nx, halfY - Ny + 1: halfY+ Ny);
figure(138)
pcolor(x_useful2,y_useful2,flux_useful)
shading flat;
axis image;
colorbar;
xlabel('x(\mum)')
ylabel('y(\mum)')
title('m\Phi_0')
% % number of points to be plotted
% NPU = round(PU_range/xP * length(xq)/2) - 1;
% halfXq = round(length(xq)/2);
% halfYq = round(length(yq)/2);
% 
% hz_PU = real_hz_useful_interp(halfXq - NPU+1: halfXq + NPU-1, halfYq - NPU+1: halfYq + NPU-1);
% fluxPU = hz_PU .*PU * dPU * dPU;
% 
% figure(358)
% imagesc(xPU,yPU,hz_PU)
% shading flat;
% axis image;
% colorbar;
% colormap jet;
% title('unit um')
% 
% figure(357)
% imagesc(xPU,yPU,fluxPU)
% shading flat;
% axis image;
% colorbar;
% colormap jet;
% title('unit um')
% %%
% % convolve with a SQUID with inner radius 3.73 um
% % and outer radius 4.56 um
% %r_eff = 4.17;
% %r_eff = 3.73;
% PickupLoop = 0*X;
% b = 3.73;
% PickupLoop((X.^2+Y.^2)<b^2)=1;
% figure(137)
% imagesc(x,y,PickupLoop)
% shading flat;
% axis image;
% colorbar;
% colormap bone;
% 
% % convolve pickup loop with field
% flux = conv2(abs(hz),PickupLoop,'same')* dx * dy* 1e3/Phi_0;
% figure(138)
% pcolor(x,y,flux)
% shading flat;
% axis image;
% colorbar;
% xlabel('x(\mum)')
% ylabel('y(\mum)')
% title('m\Phi_0')

%save('trial.mat','flux', 'x_useful','y_useful', "hz");

% directly calculate mutual inductance
%pre_M = - mu_0* alpha^2 * a* b/(4 * pi * Phi_0);
%q = integral3(@M,-inf,inf,-inf,inf,-inf,inf, 'ArrayValued', 1);
%q_array = pre_M * arrayfun(@(xp) integral3(@(kx,ky,qx) M(kx,ky,qx,xp),-100000,100000,-100000,100000,-100000,100000), x);
%q_array = pre_M * arrayfun(@(xp) integral3(@(kx,ky,qx) M(kx,ky,qx,xp),-inf, inf, - inf, inf, -inf, inf, 'RelTol',1e-5 cv), x);
%fun = @(kx,ky,qx) exp(-k*z + 1i * kx .* X -1i * ky * x0+ 1i * qx * x0 - sqrt(qx^2 + ky^2)*z0 - 1i * qx * );
% function y = M(kx, ky, qx, xp)
% Q = sqrt(qx.^2 + ky.^2);
% k = sqrt(kx.^2 + ky.^2);
% % scan height in m
% z = 4.4 * 1e-6;
% % Pearl length in m
% Lambda = 350 * 1e-6;
% % defect distance to center in m
% x0 = 3 * 1e-6;
% % field coil effective radius in m
% a = 7.35 * 1e-6;
% % pickup loop effective radius in m
% b = 3.73 * 1e-6;
% y = exp(-k*z + 1i * kx .* xp -1i .* kx * x0 + 1i .* qx * x0 - Q*z - 1i * qx .* xp) .* besselj(1, k * b) .* besselj(1, Q * a) .* (qx .* kx + ky.^2)./ (k .* (1 + Lambda * k) .* Q .* (1 + Lambda * Q));
% end
result = flux_useful;
end