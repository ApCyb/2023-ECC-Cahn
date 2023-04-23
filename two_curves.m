close all;
clc;
clear all;

tau_1 = 1;
tau_2 = 0.005;
K_vco = 120;
k = 2/pi;
omega_sep_formula(tau_1, tau_2, k, K_vco)
omega_e_free = 50;


mu = pi*k - 1;
xi = (k*tau_2*K_vco + 1)/(2*sqrt(K_vco*(tau_1 + tau_2)));
eta = (k*tau_2*K_vco - mu)/(2*sqrt(K_vco*(tau_1 + tau_2)));
rho = sqrt(abs(xi^2 - k));
kappa = sqrt(eta^2 + mu*k);


y1_ht = (kappa + eta)*(1/k + omega_e_free/(k*K_vco));
y1_max = 2*tau_1*sqrt(K_vco/(k*(tau_1 + tau_2)));
y0_ht = (kappa - eta)*(1/k - omega_e_free/(k*K_vco));
y0_max = 2*tau_1*sqrt(K_vco/(k*(tau_1 + tau_2)));


curve1 = @(y1, y0) ((y1 + (kappa - eta)*(1/k + omega_e_free/(k*K_vco)))/...
    (y0 - (kappa - eta)*(1/k - omega_e_free/(k*K_vco))))^...
    ((kappa - eta)/kappa) - ...
((y0 + (kappa + eta)*(1/k - omega_e_free/(k*K_vco)))/...
    (y1 - (kappa + eta)*(1/k + omega_e_free/(k*K_vco))))^...
    ((kappa + eta)/kappa);

    if(xi < sqrt(k))
        curve2 = @(y1, y0) ((y1^2 - 2*xi*(1/k + omega_e_free/(k*K_vco))*y1 + k*(1/k + omega_e_free/(k*K_vco))^2)/...
                 (y0^2 + 2*xi*(1/k - omega_e_free/(k*K_vco))*y0 + k*(1/k - omega_e_free/(k*K_vco))^2) - ...
        exp(2*xi/rho*(...
        atan(((1/k - omega_e_free/(k*K_vco))*rho)/(y0 + xi*(1/k - omega_e_free/(k*K_vco)))) - ...
        atan((y1 - xi*(1/k + omega_e_free/(k*K_vco)))/((1/k + omega_e_free/(k*K_vco))*rho)) + ...
        pi/2)));
    end
    if(xi == sqrt(k))
        curve2 = @(y1, y0) ((y1 - sqrt(k)*(1/k + omega_e_free/(k*K_vco)))/...
                 (y0 + sqrt(k)*(1/k - omega_e_free/(k*K_vco))) - ...
        exp(sqrt(k)*(1/k - omega_e_free/(k*K_vco))/(y0 + sqrt(k)*(1/k - omega_e_free/(k*K_vco))) + ...
            sqrt(k)*(1/k + omega_e_free/(k*K_vco))/(y1 - sqrt(k)*(1/k + omega_e_free/(k*K_vco)))));
    end
    if(xi > sqrt(k))
        curve2 = @(y1, y0) (((y1 + (rho - xi)*(1/k + omega_e_free/(k*K_vco)))/...
                  (y0 - (rho - xi)*(1/k - omega_e_free/(k*K_vco))))^((rho - xi)/rho) - ...
                 ((y0 + (rho + xi)*(1/k - omega_e_free/(k*K_vco)))/...
                  (y1 - (rho + xi)*(1/k + omega_e_free/(k*K_vco))))^((rho + xi)/rho));
    end

hold on;
fimplicit(curve1,[y1_ht y1_max y0_ht y0_max], 'color', [127/255 0 1], 'LineWidth', 2);
fimplicit(curve2,[y1_ht y1_max y0_ht y0_max], 'color', [1 165/255 0], 'LineWidth', 2);

plot((kappa + eta)*(1/k + omega_e_free/k/K_vco), (kappa - eta)*(1/k - omega_e_free/k/K_vco), 'b.', 'MarkerSize', 30);

rec_max = sqrt((tau_1 + tau_2)*K_vco)*(1 + omega_e_free/K_vco);
rec1 = [y1_ht, y0_ht];
rec2 = [y1_ht, rec_max];
rec3 = [rec_max, rec_max];
rec4 = [rec_max, y0_ht];

plot([rec1(1) rec2(1)],[rec1(2) rec2(2)],'color','b')
plot([rec2(1) rec3(1)],[rec2(2) rec3(2)],'color','b')
plot([rec3(1) rec4(1)],[rec3(2) rec4(2)],'color','b')
plot([rec4(1) rec1(1)],[rec4(2) rec1(2)],'color','b')

grid on;
set(gca,'FontSize', 15)
legend('y_0 = y_0(y_1)', 'y_2 = y_2(y_1)', 'Location','best');

axis([2 12 0.6 10])





