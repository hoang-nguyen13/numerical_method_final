%SCATTERING 1D 

mystartdefaults

recipunit = 1.0E+10;                              % reciprocal space unit [1/m]
ekinscale = ((hbar*recipunit)^2/(2*elm))/qel;      % free electron kinetic energy 

datafile = 'scat1d.dat';   
pdf = true;
latex = true;

tol = 1e-12;
% Values for GaAs

mu = 7900E-4;              % mobility 
efm = 0.067;               % effective mass
relax = mu^efm*elm/qel;   % relaxation time

lifetime = 1e-9;           % lifetime = recombination time 
Gamma = hbar*(2*pi/lifetime/qel);  % Absolute value of imaginary part of energy 

%-----------------------------------------------------
% INCIDENT ENERGY RANGE [eV] - eq (4.19) IN x<0 DOMAIN
f1 = fopen ( datafile,'w'); % output file #####
spectrum = true; tune= false;

if (~spectrum && ~tune)       % Boolean variables defined in sourcing script
  E_0 = 0.15; 
  U_0 = 0;
  U_1 = 0.10; 
  x_min = -20;
  x_max = 100;
end 

if (spectrum && ~ tune)
  E_0 = 0.15;
  U_0 = 0;                     % reference level in x<0 domain
  Bias_min = -1.5;             % min ref level in x>0 domain
  Bias_max = 0.2;              % max ref level in x>0 domain 
  Bias_step = 0.005;          % resolution of bias tunning
  U_1 = [Bias_min:Bias_step:Bias_max];
end

if (spectrum && tune) 
  error('Simultaneous spectrum and tune mode are not implemented');
end 


%--------------------------
% MODIFIED ENERGY ARRAYS 
E_0 = E_0 - U_0 ;
E_1 = E_0-U_0-U_1 ;

% Real wavevectors
k_0 = sqrt(E_0/ekinscale); 
k_1 = sqrt(E_1/ekinscale);

% Imaginary parts of incident energy and wavefunction

ck_0 = sqrt((E_0+1i*Gamma)/ekinscale);
ck_1 = sqrt((E_1+1i*Gamma)/ekinscale);  

%------------------------------------
% Discription of localized perturbation such that xp>0

xp_min = 0; 
xp_max = 80; 
n = 80; 

step = (xp_max - xp_min)/n; 

xp = zeros(n,1);
for i = 1:n
  xp(i) = xp_min + step/2 + (i-1)*step;
end 

U = zeros(n,1)+U_0; % Initialization = reference level for x<0

for i=1:n 
  if (xp(i)>0 && xp(i) < 15)
    U(i) = U(i) + 0.2;
  end
  if (xp(i)>65 && xp(i)<80)
    U(i) = U(i) + 0.1;
  end
end 

if (abs(U_0-max(U_1)) > tol) 
  electric_field = -(max(U_1)-U_0)/(xp_max-xp_min);
  for i = 1:n
    u(i) = U(i) - electric_field*xp(i);
  end
end


V = u/ekinscale;
G0 = zeros(n,n);
length(ck_0), length(ck_1)
rb = (ck_0-ck_1) / (ck_0 + ck_1); % Reflection coefficient  - eq A 64 
tb = (2*ck_0) / (ck_0 + ck_1);        % Transmission coefficient - eq A 64 
   
Phi0p = tb * exp(1i*ck_1);                 % Incident plane 
G0 = step * Green1(xp,xp',ck_0,ck_1);     % Green's function matrix inside the perturbation
   
T = eye(n,n)-G0*diag(V);                      % matrix in rq (4.51)
   
Phip = T\Phi0p;

for i =1:m
  Phis(i) = 0;
  for j=1:n 
    Phis(i) = Phis(i) + step * Green1(x(i),xp(j),ck_0,ck_1)*V(j)*Phip(j);
  end
  if (x(i)>0)  % Adding incident field
    Phip(i) = tb*exp(1i*ck_0*x(i)) + rb*exp(-1i*ck_0*x(i)) + Phis(i) ;
  end
end
Ref = abs( (rb*exp(-1i*ck_0*x(1)) + Phis(1))/ (exp(1i*ck_0*x(1))))^2;
Tra = (ck_1/ck_0)*abs(Phi(2))^2;
Absor = 1-Ref-Tra;

