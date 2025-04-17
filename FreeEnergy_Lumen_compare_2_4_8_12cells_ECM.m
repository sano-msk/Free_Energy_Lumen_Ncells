% Free Energy for Lumen and N Cells (Cyst, Organoids)
% Calculation is based on the symmetric configuraltion of N=2,4,8,12 cells
% in 3D with the asumption of apical and lateral surface tension of cells, 
% elasiticity of cell volume, and elasticity of ECM under the lumen
% pressure \Delta p (see Lu et al. in Nature Comm.(2025) for details.)

global gl k kE R0 L L0;
ga = 100;   % Surface tension of cell apical surfaces [Pa*um]
gl = 100;   % Surface tension of cell lateral surfaces [Pa*um]
k = 1/50;   % Elasticity of cell volume
kE = 1/100; % Elasticity of ECM
L = 40;
L0 = 38.5;  % Related to the natural volume of ECM,V0_ecm=g_s R0^3
R0 = 20;    % Related to the natural volume of each cell, V0=g_s R0^3
dp = 100;   % [Pa] Pressure difference between lumen and outrside
N = 250;    % Number of points in R
Rmax = 25;  % Max value of the range of R
dR = Rmax/N;
% R is the radius of the organoid
% r is the radius of the lumen
% L is the radius of the ECM

deltaP = [0,30,100];  % List of \Delta P [Pa] to calculate

for ii = 1:3
    dp = deltaP(ii);
    
% 2 cells
for i = 1:N
    r(i) = dR*(i-1);
    x0 = 20;  % 30 initial value for R
    % Solv dF/dR = 0 for each r to find R_min (radius of organoid)
    fun = @(x)paramfun(x,r(i)); % fun now has the new r value
    R(i) = fsolve(fun,x0);
end

figure(1)
plot(r,R,'LineWidth',2)
hold on

% Free Energy for 2 cells
G2 = -dp*(4*pi/3)*r.^3 + 4*pi*ga*r.^2 + 2*gl*pi*(R.^2 - r.^2)  +k/2*(4*pi/3*R.^3-4*pi/3*r.^3-4*pi/3*R0.^3).^2 ...
     +kE*(4*pi/3*L^3-4*pi/3*R.^3-4*pi/3*L0^3).^2;
figure(ii+1)
plot(r, G2, 'LineWidth',2.5)
hold on
i2_min = find(G2==min(G2));
R2_min = dR*(i2_min - 1);
iR(ii,1) = 2;
r_min(ii,1) = R2_min;

% 4 cells
for i = 1:N
    r(i) = dR*(i-1);
    % Solv dF/dR = 0 for each r to find R_min (radius of organoid)
    fun4 = @(x)paramfun4(x,r(i)); % fun now has the new r value
    R(i) = fsolve(fun4,x0);
end

figure(1)
plot(r,R,'LineWidth',2)
hold on
% Free Energy for 4 cells
G4 = -dp*(4*pi/3)*r.^3 + 4*pi*ga*r.^2 + 4*gl*pi*(R.^2 - r.^2)  +k/4*(4*pi/3*R.^3-4*pi/3*r.^3-4*pi/3*R0.^3).^2 ...
     +kE*(4*pi/3*L^3-4*pi/3*R.^3-4*pi/3*L0^3).^2;
figure(ii+1)
plot(r, G4, 'LineWidth',2.5)
hold on
i4_min = find(G4==min(G4));
R4_min = dR*(i4_min - 1);
iR(ii,2) = 4;
r_min(ii,2) = R4_min;

% 8 cells
for i = 1:N
    r(i) = dR*(i-1);
    % Solv dF/dR = 0 for each r to find R_min (radius of organoid)
    fun8 = @(x)paramfun8(x,r(i)); % fun now has the new r value
    R(i) = fsolve(fun8,x0);
end

figure(1)
plot(r,R,'LineWidth',2.5)
hold on

figure(ii+1)
% Free Energy for 8 cells
G8 = -dp*(4*pi/3)*r.^3 + 4*pi*ga*r.^2 + 6*gl*pi*(R.^2 - r.^2)  +(1/8)*k*(4*pi/3*R.^3-4*pi/3*r.^3-4*pi/3*R0.^3).^2 ...
     + kE*(4*pi/3*L^3-4*pi/3*R.^3-4*pi/3*L0^3).^2;
plot(r, G8, 'LineWidth',2.5)
hold on
i8_min = find(G8==min(G8));
R8_min = dR*(i8_min - 1);
iR(ii,3) = 8;
r_min(ii,3) = R8_min;

% 12 cells (Regular dodecahedron)
%  R (solution x) corresponds to the radius of circumscribe sphere
%  RI corresponds to theradius of inscribe sphere
for i = 1:N
    r(i) = dR*(i-1);
    % Solv dF/dR = 0 for each r to find R_min (radius of organoid)
    fun12 = @(x)paramfun12(x,r(i)); % fun now has the new r value
    R(i) = fsolve(fun12,x0);
end

g = (1 + sqrt(5))/2;  % Golden Mean
%a = 2*R/sqrt(3)/g;
RI = (3*g + 2)/sqrt(4*g + 3)/2*2.*R/sqrt(3)/g;

figure(1)
plot(r,(RI+R)/2,'LineWidth',2.5)  % Plot the mean of R and RI
axis square
xlabel('r [\mu m]','FontSize',18)
ylabel('R_{min} [\mu m]','FontSize',18)
hold off

lambda = 2/sqrt(3)/g;
A1 = (7*g+4)/2;
G0 = 0.0;   % Offset for the free energy
figure(ii+1)
% Free Energy for 12 cells (Regular dodecahedron)
G12 = -dp*A1*lambda^3*r.^3 + 3*sqrt(5*(4*g+3))*ga*r.^2 + 15*gl*sqrt((3*g^2-1))*lambda^2*(R.^2 - r.^2)...
    +(1/12)*k*lambda^6*A1^2*(R.^3-r.^3-R0.^3).^2 ...
     + kE*(4*pi/3*L^3-A1*lambda^3*R.^3-4*pi/3*L0^3).^2 - G0;
plot(r, G12, 'LineWidth',2.5)
if ii==3
    ylim([-15e5,14e5]);
else
    ylim([2e5,14e5]);    
end
%axis square
xlabel('lumen radius r [\mum]','FontSize',18)
ylabel('F(r)','FontSize',18)
legend('N=2','N=4','N=8','N=12','Location','southeast','FontSize',12)
str = append('\Deltap = ',num2str(dp),'[Pa]');
title(str,'FontSize',14)
hold off
i12_min = find(G12==min(G12));
R12_min = dR*(i12_min - 1);
iR(ii,4) = 12;
r_min(ii,4) = R12_min;
axis square
box on

figure(5)
scatter(iR(ii,:), r_min(ii,:),'LineWidth',2)
hold on
xlabel('Number of Cells','FontSize',14)
ylabel('lumen radius [\mum]','FontSize',14)
axis square

end
figure(5)
legend('0 [Pa}','30 [Pa]','100 [Pa]')

function F = paramfun(x,r)
global gl k kE R0 L L0;
F = gl + k*4*pi/3*(x.^4-(r^3+R0^3).*x) + 8*kE*pi/3*(x.^4-(L^3-L0^3).*x);
end

function F = paramfun4(x,r)
global gl k kE R0 L L0;
F = gl + k*pi/3*(x.^4-(r^3+R0^3).*x) + 4*kE*pi/3*(x.^4+(L0^3-L^3).*x);
end

function F = paramfun8(x,r)
global gl k kE R0 L L0;
F = gl + k*pi/9*(x.^4-(r^3+R0^3).*x) + 8*kE*pi/9*(x.^4+(L0^3-L^3).*x);
end

function F = paramfun12(x,r)
global gl k kE R0 L L0;
%  R (solution x) corresponds to circumscribe sphere
g = (1+sqrt(5))/2;  % Golden mean
lambda = 2/sqrt(3)/g;
A1 = (7*g + 4)/2;
lambda = 2/sqrt(3)/g;
F = 30*sqrt(3*g^2-1)*lambda^2*gl*x + k/2*lambda^6*A1^2*(x.^5-(r^3+R0^3).*x.^2)...
    + 2*kE*(A1*lambda^3*x.^5+4*pi/3*(L0^3-L^3).*x.^2)*3*A1*lambda^3;
end
