% Programa para interpolación de Lagrange y Newton (diferencias progresivas)
% con estimación de error - Datos de Carbono-14
clear; clc;

% Datos de la tabla
N14C = [0.78; 0.80; 0.82; 0.84; 0.86; 0.88; 0.90; 0.92];
Antiguedad = [2050; 1850; 1650; 1450; 1250; 1050; 870; 690];

% Valor a interpolar (puedes cambiarlo)
x_interp = 0.85; % Ejemplo: interpolar en N14C = 0.85

fprintf(’========================================\n’);
fprintf(‘INTERPOLACIÓN CON ESTIMACIÓN DE ERROR\n’);
fprintf(’========================================\n\n’);

%% INTERPOLACIÓN 1: i = 2,3,4,5 (índices MATLAB: 3,4,5,6)
fprintf(‘INTERPOLACIÓN 1: Datos i=2,3,4,5\n’);
fprintf(’––––––––––––––––––––\n’);
indices1 = 3:6; % En MATLAB los índices empiezan en 1
indice_extra1 = 7; % i=6 para el error

x1 = N14C(indices1);
y1 = Antiguedad(indices1);
h1 = x1(2) - x1(1); % Paso (debe ser constante)

% Verificar que el paso es constante
if max(abs(diff(x1) - h1)) > 1e-10
warning(‘El paso no es constante en interpolación 1’);
end

fprintf(‘Paso h = %.4f\n\n’, h1);

% Lagrange con 4 puntos
y_lag1 = interpolacion_lagrange(x1, y1, x_interp);
fprintf(‘Lagrange (4 puntos): %.4f años\n’, y_lag1);

% Lagrange con 5 puntos (para calcular error)
x1_5pts = N14C([indices1, indice_extra1]);
y1_5pts = Antiguedad([indices1, indice_extra1]);
y_lag1_5pts = interpolacion_lagrange(x1_5pts, y1_5pts, x_interp);
error_lag1 = y_lag1_5pts - y_lag1;
fprintf(‘Error estimado Lagrange: %.4f años\n\n’, abs(error_lag1));

% Newton con diferencias progresivas
[tabla_diff1, coef_newton1] = diferencias_progresivas(y1);
s1 = (x_interp - x1(1)) / h1; % Variable normalizada
y_newton1 = evaluar_newton_progresivo(coef_newton1, s1);
fprintf(‘Newton progresivo (4 puntos): %.4f años\n’, y_newton1);
fprintf(‘Variable s = (x - x0)/h = %.6f\n’, s1);

% Calcular diferencias progresivas con 5 puntos para el error
[tabla_diff1_5pts, ~] = diferencias_progresivas(y1_5pts);
delta4_f0 = tabla_diff1_5pts(1, 5); % Δ⁴f₀
termino_s1 = prod(s1 - (0:3)) / factorial(4); % s(s-1)(s-2)(s-3)/4!
error_newton1 = delta4_f0 * termino_s1;
fprintf(‘Δ⁴f₀ = %.4f\n’, delta4_f0);
fprintf(‘s(s-1)(s-2)(s-3)/4! = %.6f\n’, termino_s1);
fprintf(‘Error estimado Newton: %.4f años\n\n’, abs(error_newton1));

%% INTERPOLACIÓN 2: i = 3,4,5,6 (índices MATLAB: 4,5,6,7)
fprintf(’\n========================================\n’);
fprintf(‘INTERPOLACIÓN 2: Datos i=3,4,5,6\n’);
fprintf(’––––––––––––––––––––\n’);
indices2 = 4:7;
indice_extra2 = 8; % i=7 para el error

x2 = N14C(indices2);
y2 = Antiguedad(indices2);
h2 = x2(2) - x2(1);

% Verificar que el paso es constante
if max(abs(diff(x2) - h2)) > 1e-10
warning(‘El paso no es constante en interpolación 2’);
end

fprintf(‘Paso h = %.4f\n\n’, h2);

% Lagrange con 4 puntos
y_lag2 = interpolacion_lagrange(x2, y2, x_interp);
fprintf(‘Lagrange (4 puntos): %.4f años\n’, y_lag2);

% Lagrange con 5 puntos (para calcular error)
x2_5pts = N14C([indices2, indice_extra2]);
y2_5pts = Antiguedad([indices2, indice_extra2]);
y2_lag2_5pts = interpolacion_lagrange(x2_5pts, y2_5pts, x_interp);
error_lag2 = y2_lag2_5pts - y_lag2;
fprintf(‘Error estimado Lagrange: %.4f años\n\n’, abs(error_lag2));

% Newton con diferencias progresivas
[tabla_diff2, coef_newton2] = diferencias_progresivas(y2);
s2 = (x_interp - x2(1)) / h2;
y_newton2 = evaluar_newton_progresivo(coef_newton2, s2);
fprintf(‘Newton progresivo (4 puntos): %.4f años\n’, y_newton2);
fprintf(‘Variable s = (x - x0)/h = %.6f\n’, s2);

% Calcular diferencias progresivas con 5 puntos para el error
[tabla_diff2_5pts, ~] = diferencias_progresivas(y2_5pts);
delta4_f0_2 = tabla_diff2_5pts(1, 5);
termino_s2 = prod(s2 - (0:3)) / factorial(4);
error_newton2 = delta4_f0_2 * termino_s2;
fprintf(‘Δ⁴f₀ = %.4f\n’, delta4_f0_2);
fprintf(‘s(s-1)(s-2)(s-3)/4! = %.6f\n’, termino_s2);
fprintf(‘Error estimado Newton: %.4f años\n\n’, abs(error_newton2));

%% Mostrar tablas de diferencias progresivas
fprintf(’\n========================================\n’);
fprintf(‘TABLA DE DIFERENCIAS PROGRESIVAS - Interpolación 1\n’);
fprintf(’––––––––––––––––––––\n’);
mostrar_tabla_progresivas(x1_5pts, y1_5pts, tabla_diff1_5pts);

fprintf(’\n========================================\n’);
fprintf(‘TABLA DE DIFERENCIAS PROGRESIVAS - Interpolación 2\n’);
fprintf(’––––––––––––––––––––\n’);
mostrar_tabla_progresivas(x2_5pts, y2_5pts, tabla_diff2_5pts);

%% FUNCIONES

function y_interp = interpolacion_lagrange(x, y, x_val)
n = length(x);
y_interp = 0;
for i = 1:n
L = 1;
for j = 1:n
if i ~= j
L = L * (x_val - x(j)) / (x(i) - x(j));
end
end
y_interp = y_interp + y(i) * L;
end
end

function [tabla, coef] = diferencias_progresivas(y)
n = length(y);
tabla = zeros(n, n);
tabla(:, 1) = y;

```
% Calcular diferencias progresivas
for j = 2:n
    for i = 1:(n-j+1)
        tabla(i, j) = tabla(i+1, j-1) - tabla(i, j-1);
    end
end

% Coeficientes son las diferencias en la primera fila
coef = tabla(1, :);
```

end

function y = evaluar_newton_progresivo(coef, s)
n = length(coef);
y = coef(1);

```
for k = 2:n
    % Calcular s(s-1)(s-2)...(s-k+2) / (k-1)!
    termino = 1;
    for j = 0:(k-2)
        termino = termino * (s - j);
    end
    termino = termino / factorial(k-1);
    
    y = y + coef(k) * termino;
end
```

end

function mostrar_tabla_progresivas(x, y, tabla)
n = length(x);

```
fprintf('   i      x_i        f_i        Δf_i      Δ²f_i     Δ³f_i     Δ⁴f_i\n');
fprintf('-----------------------------------------------------------------------------\n');
for i = 1:n
    fprintf('  %2d   %.4f   ', i-1, x(i));
    for j = 1:n
        if i+j-1 <= n
            fprintf(' %10.2f ', tabla(i,j));
        end
    end
    fprintf('\n');
end
```

end