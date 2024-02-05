mystartdefaults
% Energy range and step size
E_min = 0; % in eV
E_max = 0.3; % in eV
deltaE = 0.0005; % in eV
step_size = 0.5;
x1_min = 0;
x1_max = 80;
tol = 1e-12;
tau = 1e-9; % 1 ns in seconds
gamma=(hbar*2*pi/tau)/qel; % damping factor
recipunit=1.0E+10;
ekinscale=(hbar*recipunit)^2/(2*elm)/qel;

%plotting conditions
unbiased = false;
prob = false;
roughly_biased = true;
biased_junction = false;
roughly_biased_01 = false;
U_minus_02 = false;
U_plus_02 = false;
varying_bias = false;

n = (x1_max - x1_min)/step_size;
U = zeros(1,n); % Potential barrier
x1 = zeros(1,n);

for i = 1:n
    x1(i) = x1_min + step_size/2 + (i - 1) * step_size;
    if x1(i) <= 15
        U(i) = 0.2;
    elseif x1(i) >= 65 && x1(i) <= 80
        U(i) = 0.2;
    else
        U(i) = 0;
    end
end
edivide = round((E_max-E_min)/deltaE);
for i = 1:edivide
    E0(i) = deltaE/2 + E_min+(i-1) * deltaE;
end

for i = 1:edivide
    [R_values(i),T_values(i),A_values(i)]=spectra(0,E0(i),gamma,x1,U,step_size,ekinscale) ; % the step is set to 0 to look at the constant background case
end

if unbiased == true
    plot(E0, R_values, E0, T_values, E0, A_values, 'LineWidth', 2); % Plotting
    xlabel('Energy(eV)');
    ylabel('Coefficient');
    title('Unbiased junction');
    legend('Reflection', 'Transmission', 'Absorption', "FontSize",10, "Location","northwest");
    legend_position = [0.25, 0.5, 0.1, 0.1];  % [left, bottom, width, height]
    set(gca, 'Position', get(gca, 'Position') + [0.05, 0, 0, 0]);  % Adjust plot position
    set(legend, 'Units', 'normalized', 'Position', legend_position);
    % saveas(gcf, 'unbiased_junction_spectra.png');  % You can specify other file formats like '.jpg', '.pdf', etc.
end

x_min = -20;
x_max = 100;
n = (x_max - x_min)/step_size;
n1 = (x1_max - x1_min)/step_size;
U_total = zeros(1,n); % Potential barriers
x = zeros(1,n);
for i = 1:n
    x(i) = x_min + step_size/2 + (i - 1) * step_size;
    if x(i) >= 0 && x(i) <= 15
        U_total(i) = 0.2;
    elseif x(i) >= 65 && x(i) <= 80
        U_total(i) = 0.2;
    else
        U_total(i) = 0;
    end
end
U = zeros(1,n1); % Potential barrier
x1 = zeros(1,n1);
for i = 1:n1
    x1(i) = x1_min + step_size/2 + (i - 1) * step_size;
    if x1(i) <= 15
        U(i) = 0.2;
    elseif x1(i) >= 65 && x1(i) <= 80
        U(i) = 0.2;
    else
        U(i) = 0;
    end
end

E_resonance=[0.01075,0.04325,0.09525,0.16325,0.1312];
[~,n]=size(E_resonance);
W = zeros(length(U),length(U));

wave_func = zeros(n-1, length(U_total));
for i = 1:n
    k1 = sqrt((E_resonance(i) + 1i * gamma)/ekinscale);
    for j = 1:length(U)
        phi0k(j) = exp(1i * k1 * x1(j));
        W(j,j) = step_size * U(j)/ekinscale;
        for k = 1:length(U)
            G0(j,k) = Green0(0,x1(j),x1(k),E_resonance(i),gamma,ekinscale);
        end
    end
    T = eye(length(U)) - G0 * W;
    phisol = T\(phi0k.');
    for j = 1:(length(U_total))
        if i < n
            wave_func(i,j) = exp(x(j) * 1i * k1) + outside_sol(x(j), x1, U, step_size, phisol, 0, E_resonance(i), gamma, ekinscale);
        else
            non_resonance_phi(j) = exp(x(j) * 1i * k1) + outside_sol(x(j),x1,U,step_size,phisol,0,E_resonance(i),gamma,ekinscale);
        end
    end
end
resonance_prob = abs(wave_func).^2;
non_resonance_prob = abs(non_resonance_phi).^2;
for i=1:length(U_total)    
  resonance_prob(1,i) = 0.0059 * resonance_prob(1,i) + E_resonance(1);
  resonance_prob(2,i) = 0.0019 * resonance_prob(2,i) + E_resonance(2);
  resonance_prob(3,i) = 0.0005 * resonance_prob(3,i) + E_resonance(3);
  resonance_prob(4,i) = 0.0030 * resonance_prob(4,i) + E_resonance(4); 
  non_resonance_prob(i) = 0.012 * non_resonance_prob(i) + E_resonance(n);
end

if prob == true
    plot(x,U_total,'LineWidth',3);
    title('Static Potential Barriers', 'FontSize',20);
    ylim([0, 0.25]);
    xlabel('x(Å)');
    ylabel('E(eV)');
    % saveas(gcf, 'Static_Potential_Barriers.png');
    % 
    plot(x, non_resonance_prob, x,resonance_prob, x, U_total, 'LineWidth',2);
    title('Probability densities against barriers','FontSize',20);
    ylim([0, 0.27]);
    xlabel('x(Å)');
    ylabel('E(eV)'); 
    % saveas(gcf, 'prob_density.png')
end
% Energy range and step size
E_min = 0; % in eV
E_max = 0.3; % in eV
deltaE = 0.0005; % in eV
step_size = 0.5;
x1_min = 0;
x1_max = 80;

n = (x1_max - x1_min)/step_size;

U = zeros(1, n); % Potential barrier
x1 = zeros(1, n);
for i = 1:n
    x1(i) = x1_min + step_size/2 + (i - 1) * step_size;
    if x1(i) <= 15
        U(i) = 0.2;
    elseif x1(i) >= 65 && x1(i) <= 80
        U(i) = 0.1;
    else
        U(i) = 0;
    end
end

edivide = round((E_max-E_min)/deltaE);
for i = 1:edivide
    E0(i) = deltaE/2 + E_min+(i-1) * deltaE;
end

for i = 1:edivide
    [R_values(i),T_values(i),A_values(i)]=spectra(0,E0(i),gamma,x1,U,step_size,ekinscale) ; % the step is set to 0 to look at the constant background case
end
if roughly_biased == true
    plot(E0, R_values, E0, T_values, E0, A_values, 'LineWidth', 2);
    xlabel('Energy (eV)');
    ylabel('Coeficient');
    title('Second barrier lowered (0.1 eV)');
    legend('Reflection', 'Transmission', 'Absorption', "FontSize",10, "Location","northwest");
    legend_position = [0.35, 0.5, 0.1, 0.1];  % [left, bottom, width, height]
    set(gca, 'Position', get(gca, 'Position') + [0.05, 0, 0, 0]);  % Adjust plot position
    set(legend, 'Units', 'normalized', 'Position', legend_position);
    % saveas(gcf, 'roughly_biased_junction_spectra.png');
end
% 
% Constants and parameters

% Parameters for the step function potential
x_min = -20; % Minimum x in Angstroms
x_max = 100; % Maximum x in Angstroms
U0 = 0; % Potential for x < 0
U1 = -0.1; % Step potential in eV
% 
% Discretize the space
n = (x_max - x_min)/step_size;
% Define the initial step function potential
U = zeros(n, 1); % Initialize potential array
% Define the region where the electric field is applied
x_prime_min = 0; % Starting point of the heterostructure in Angstroms
x_prime_max = 80; % Ending point of the heterostructure in Angstroms

% Modify the potential profile to include the electric field effect
W = zeros(n, 1); % Initialize modified potential array
U = zeros(n, 1); % Potential barrier
x = zeros(n, 1);

for i = 1:n
    x(i) = x_min + step_size/2 + (i - 1) * step_size;
    if x(i) <= 15 && x(i) > 0
        U(i) = 0.2;
    elseif x(i) >= 65 && x(i) <= 80
        U(i) = 0.2;
    elseif x(i) >= 80
        U(i) = -0.1;
    else
        U(i) = 0;
    end
end
%plot(x, U, 'LineWidth', 2);
%ylim([-0.1 0.25]);
%plot(x, V, 'LineWidth', 2);
%Calculate the magnitude of the electric field
E = - U1 / (x_prime_max - x_prime_min);

for i = 1:n
    if x(i) > x_prime_min && x(i) < x_prime_max
        W(i) = U(i) - E * x(i);
    elseif x(i) >= x_prime_max
        W(i) = -0.1;
    else
        W(i) = 0;
    end
end
% 
% Plot the modified potential profile
if biased_junction == true
    plot(x, W, 'LineWidth', 3);
    title('Modified Potential Profile with W(x\prime)');
    xlabel('x(Å)');
    ylabel('E(eV)');
    ylim([-0.1 0.25]);
    % saveas(gcf, 'biased_junction.png');
end
E_min = 0; % in eV
E_max = 0.3; % in eV
deltaE = 0.0005; % in eV
x1_min = 0;
x1_max = 80;
n = (x1_max - x1_min)/step_size;

xp = zeros(1,n); % done
for i = 1:n % done
  xp(i) = x1_min + step_size/2 ... 
      + (i - 1) * step_size;                       % done
end 
U = zeros(1,n) + U0;                               % done
x1 = zeros(1,n);
for i=1:n
  x1(i) = x1_min + step_size/2 + (i - 1) * step_size;
  if (xp(i) > 0 && xp(i) < 15)                     % done
    U(i) = U(i) + 0.2;
  end
  if (xp(i) > 65 && xp(i) < 80)                    % done
    U(i) = U(i) + 0.2;
  end
end
if (abs(U0 - max(U1)) > tol)                       % done
  E_field = -(max(U1) - U0)/ ...
  (x_prime_max - x_prime_min);                     % done
  for i = 1:n
    W_uu(i) = U(i) - E_field * xp(i);              % done
  end
end
W_uu = W_uu - U1;

edivide = round((E_max-E_min)/deltaE);
for i = 1:edivide
    E0(i) = deltaE/2 + E_min+(i-1) * deltaE;
end
for i = 1:edivide
    [R_values(i),T_values(i),A_values(i)]=spectra(U1,E0(i),gamma,x1,W_uu,step_size,ekinscale) ;% the step is set to 0 to look at the constant background case
end

if roughly_biased_01 == true
    plot(E0, R_values, E0, T_values, E0, A_values, 'LineWidth', 2);
    xlabel('Energy (eV)');
    ylabel('Coeficient');
    title('Roughly biased junction (-0.1eV)');
    legend('Reflection', 'Transmission', 'Absorption', "FontSize",10, "Location","northwest");
    legend_position = [0.25, 0.5, 0.1, 0.1];  [left, bottom, width, height]
    set(gca, 'Position', get(gca, 'Position') + [0.05, 0, 0, 0]);  Adjust plot position
    set(legend, 'Units', 'normalized', 'Position', legend_position);
    % saveas(gcf, 'roughly_biased_junction_spectra(-0.1eV).png');
end
% 
% Define the x array
x = linspace(x_min, x_max, n);
% Bias values
U1_values = [-0.2, 0.2];

% Initialize arrays to store the modified potentials
W_minus_02 = zeros(1, n);
W_plus_02 = zeros(1, n);

% Loop over the bias values
for idx = 1:length(U1_values)
    U1 = U1_values(idx);

    % Potential barrier
    U = zeros(1, n);
    U(x > 0 & x <= 15) = 0.2;
    U(x >= 65 & x <= 80) = 0.2;

    % Calculate the electric field and the modified potential
    E = -U1 / (x_prime_max - x_prime_min);
    W = U;
    W(x > x_prime_min & x < x_prime_max) = U(x > x_prime_min & x < x_prime_max) - E * x(x > x_prime_min & x < x_prime_max);

    % Store the modified potential for each bias
    if U1 == -0.2
        W_minus_02 = W;
    elseif U1 == 0.2
        W_plus_02 = W;
    end
end

for i = 1:length(W_minus_02)
    if x(i) >= 80
        W_minus_02(i) = -0.2;
    end
end

% Plot the modified potential profile for U1 = -0.2 eV
if U_minus_02 == true
    plot(x, W_minus_02, 'LineWidth', 3);
    title('U1 = -0.2 eV');
    xlabel('x(Å)');
    ylabel('E(eV)');
    ylim([-0.25 0.25]);
    %saveas(gcf, 'U1_-0.2.png');
end
for i = 1:length(W_plus_02)
    if x(i) >= 80
        W_plus_02(i) = 0.2;
    end
end
if U_plus_02 == true
    % Plot the modified potential profile for U1 = 0.2 eV
    plot(x, W_plus_02, 'LineWidth', 3);
    title('U1 = 0.2 eV');
    ylabel('E(eV)');
    xlabel('x(Å)');
    ylim([-0.05 0.45]);
    % saveas(gcf, 'U1_0.2.png');
end
E0 = 0.01;
U0 = 0;
U_min = -0.2;
U_max = 0.2;                                  
deltaU = 0.0005;                                  
U1 = U_min:deltaU:U_max;
x1_min = 0;
x1_max = 80;
tau=1.0E-9; % Lifetime

Udivide = round((U_max-U_min)/deltaU);
for i = 1:Udivide
  Eners(i) = deltaU/2 + U_min + (i-1) * deltaU; 
end

n = (x1_max - x1_min)/step_size;
for i=1:n
    x1(i) = x1_min + step_size/2 + (i - 1) * step_size;
end
for i = 1:length(U1)
    E(i) = -U1(i)/(x_prime_max - x_prime_min);
end
U = zeros(1,n);
xp = zeros(1,n);
for i=1:n
  xp(i) = x_prime_min + step_size/2 + (i - 1) * step_size;
  if (xp(i) > 0 && xp(i) < 15) 
    U(i) = 0.2;
  end
  if (xp(i) > 65 && xp(i) < 80)  
    U(i) = 0.2;
  end
end

for i = 1:Udivide
    for j = 1:n
        efield(j) = Eners(i)*x1(j)/(x1_max-x1_min) - Eners(i);
    end
    W_uu = U + efield;
    [R_values(i),T_values(i),A_values(i)]=spectra(Eners(i),E0,gamma,x1,W_uu,step_size,ekinscale) ; % the step is set to 0 to look at the constant background case
end
if varying_bias == true
    plot(-Eners, R_values, -Eners, T_values, -Eners, A_values); % Plotting
    xlabel('Energy (eV)');
    ylabel('Coeficient');
    title('Varying Bias');
    legend('Reflection', 'Transmission', 'Absorption', "FontSize",10, "Location","northwest");
    legend_position = [0.25, 0.5, 0.1, 0.1];  % [left, bottom, width, height]
    set(gca, 'Position', get(gca, 'Position') + [0.05, 0, 0, 0]);  % Adjust plot position
    set(legend, 'Units', 'normalized', 'Position', legend_position);
    % saveas(gcf, 'varying_bias.png');
end






function [Gf] = Green0(step,x,xp,E,gamma,ekinscale) 
    k0 = sqrt((E + 1i * gamma)/ekinscale);
    k1 = sqrt((E + 1i * gamma - step)/ekinscale);
    if(x >= 0 && xp >= 0)
        Gf = exp(1i * k1 * abs(x-xp))/(2i * k1) + exp(1i * k1 * (x + xp))*((k1 - k0)/(k1 + k0))/(2i * k1);
    end
    if (x < 0 && xp >= 0)
        Gf = exp(-1i * k0 * x + 1i * k1 * xp)/(1i * (k0 + k1));
    end
    if (xp < 0 && x >= 0)
        Gf = exp(-1i * k0 * xp + 1i * k1 * x)/(1i * (k0 + k1));
    end
end
function [extra]= outside_sol(x,pos,V,delta_x,phisol,step,E,gamma,ekinscale)
    extra=0;
    [~,n]=size(pos);
    for i = 1:n
        extra = extra + Green0(step,x,pos(i),E,gamma,ekinscale) * delta_x * (V(i)/ekinscale) * phisol(i);
    end
end
function [R_values,T_values,A_values] = spectra(step,E,gamma,pos,V,dx,ekinscale)
   [~,n]=size(pos);
   k0 = sqrt((E + 1i * gamma)/ekinscale);
   k1 = sqrt((E + 1i * gamma - step)/ekinscale); 
   tb = 2 * k0/(k0 + k1);
   rb =(k0 - k1)/(k1 + k0);
   W = zeros(n,n);
   for i=1:n
    W(i,i) = dx * V(i)/ekinscale;
    phi0k(i)=tb * exp(1i * k1 * pos(i));
    for j=1:n
        G0(i,j)= Green0(step,pos(i),pos(j),E,gamma,ekinscale);
    end
   end
   T = eye(n)- G0 * W;
   phisol = T \ (phi0k.');
   
   phirefl = rb * exp(-1i * k0 * (pos(1)-1)) + outside_sol(pos(1) - 1,pos,V,dx,phisol,step,E,gamma,ekinscale);
   phitrans = tb * exp(1i * k1 * (pos(n)+1)) + outside_sol(pos(n) + 1,pos,V,dx,phisol,step,E,gamma,ekinscale);
   
   T_values = (real(k1)/real(k0)) * abs(phitrans)^2;
   R_values = (abs(phirefl/(exp(1i * k0 * (pos(1)-1)))))^2;
   A_values = 1 - (R_values + T_values);
end
