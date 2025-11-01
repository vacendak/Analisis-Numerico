%% Interpolación de Lagrange con propagación del error
% El programa solicita:
%   - x_datos: valores en el eje X (nodos de interpolación)
%   - y_datos: valores tabulados asociados a cada x
%   - y_error: errores (incertidumbres experimentales) asociados a cada y
%   - x_eval: punto donde se desea interpolar
%
% Calcula:
%   - P(x_eval)  valor interpolado
%   - sigma_P    error propagado en ese punto
%
% Aplica redondeo:
%   * La incertidumbre se expresa con 1 cifra significativa
%   * Si empieza por 1 o 2 se permiten 2 cifras significativas
%   * El valor interpolado se redondea a los mismos decimales

clear;

fprintf("Introducir los valores de x (entre corchetes y separados por espacios):\n");
x_datos = input("x_datos = ");        % Ejemplo: [0.84 0.86 0.88]

fprintf("\nIntroducir los valores de y (entre corchetes y separados por espacios):\n");
y_datos = input("y_datos = ");        % Ejemplo: [1450 1250 1050]

fprintf("\nIntroducir el error de cada y (entre corchetes y separados por espacios):\n");
y_error = input("y_error = ");        % Ejemplo: [10 10 10]

fprintf("\nIntroducir el punto donde interpolar x_eval (valor escalar):\n");
x_eval = input("x_eval = ");          % Ejemplo: 0.8705


%% Comprobación de tamaños y estructura
if ~isvector(x_datos) || ~isvector(y_datos) || ~isvector(y_error)
    error("Las entradas deben ser vectores.");
end

x_datos = x_datos(:).';   % Se fuerzan a vector fila
y_datos = y_datos(:).';
y_error = y_error(:).';

n = numel(x_datos);

if numel(y_datos) ~= n || numel(y_error) ~= n
    error("x_datos, y_datos y y_error deben tener la MISMA longitud.");
end


%% Cálculo de los coeficientes L_i(x) de Lagrange en x_eval
L = pesos_lagrange(x_eval, x_datos);


%% Valor interpolado y error propagado
P = sum(y_datos .* L);
sigma_P = sqrt( sum( (L .* y_error).^2 ) );


%% Redondeo según norma de laboratorio
[P_red, sigma_red, num_decimales] = redondear_incertidumbre(P, sigma_P);


%% Salida formateada
fprintf("\n=============== RESULTADOS ===============\n");
fprintf("Valor sin redondeo:      P(x_eval) = %.12g\n", P);
fprintf("Error sin redondeo:     sigma_P    = %.12g\n\n", sigma_P);

% Construye el formato con el número de decimales
fmt = sprintf('%%.%df', num_decimales);

% Crea el formato completo como un solo string
formato_valor = sprintf('Valor con redondeo:      P(x) = %s\n', fmt);
formato_error = sprintf('Error con redondeo:      sigma_P = %s\n\n', fmt);

fprintf(formato_valor, P_red);
fprintf(formato_error, sigma_red);

% Si quieres la salida conjunta final
formato_final = sprintf('P(x_eval) = %s +- %s (años)\n', fmt, fmt);
fprintf(formato_final, P_red, sigma_red);


%% Mostrar contribuciones individuales al error
fprintf("\nContribución de cada dato |L_i(x_eval)| * error_y_i:\n");
for i = 1:n
    fprintf("  Nodo %d: x = %g   L_i = %.6f   error_y = %g   contrib = %.6f\n", ...
        i, x_datos(i), L(i), y_error(i), abs(L(i) * y_error(i)));
end
fprintf("==========================================\n");


%% ========= FUNCIONES AUXILIARES =========

function L = pesos_lagrange(x_eval, xnod)
% Calcula el vector L_i(x_eval) para los nodos de interpolación xnod
%
% Cada coeficiente L_i es:
%      L_i(x) = Prod (x - x_j) / (x_i - x_j)   para j != i

    n = numel(xnod);
    L = ones(1, n);

    for i = 1:n
        for j = 1:n
            if j ~= i
                L(i) = L(i) * (x_eval - xnod(j)) / (xnod(i) - xnod(j));
            end
        end
    end
end


function [valor_red, sigma_red, nd] = redondear_incertidumbre(valor, sigma)
% Redondeo estándar de laboratorio:
% * 1 cifra significativa en la incertidumbre (sigma)
% * Si empieza por 1 o 2 -> permitir 2 cifras significativas
% * El valor interpolado se redondea a los mismos decimales

    sigma = abs(sigma);

    if sigma == 0
        sigma_red = 0;
        valor_red = valor;
        nd = 6;
        return;
    end

    expo = floor(log10(sigma));           % orden de magnitud
    primer = floor(sigma / 10^expo);      % primer dígito significativo

    if primer == 1 || primer == 2
        cifras = 2;
    else
        cifras = 1;
    end

    nd = max(0, -expo + (cifras - 1));    % número de decimales
    sigma_red = round(sigma, nd);
    valor_red = round(valor, nd);
end
