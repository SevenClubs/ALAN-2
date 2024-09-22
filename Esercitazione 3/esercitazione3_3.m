% GRUPPO composto da:
% Lorenzo Aliotta, 5655762
% Riccardo Dal Seno, 5605031
% Teresa de Jesus Fernandes, 4190022


% ESERCITAZIONE 3 (autovalori e autovettori)
% Esercizio 3

% Gestione pulizia workspace
clc
clear
close all

% Salvataggio dell'output su file
diary('3.3_output.txt');
diary on;

% Definizione della matrice A
A = [1 -1 2; -2 0 5; 6 -3 6];

% Parametri comuni
tolleranza = 1e-6; % tolleranza per la convergenza
max_iterazioni = 1000; % numero massimo di iterazioni

% Vettore iniziale v1 = [1; 1; 1]
v1 = [1; 1; 1];

% Vettore per memorizzare la norma delle differenze
norm_diff_v1 = zeros(max_iterazioni, 1);

% Iterazione per il metodo delle potenze con v1
for k = 1:max_iterazioni
    w = A * v1;
    v1_new = w / norm(w); % normalizzazione del vettore
    norm_diff_v1(k) = norm(v1_new - v1); % memorizza la norma della differenza
    if norm_diff_v1(k) < tolleranza
        norm_diff_v1 = norm_diff_v1(1:k); % riduce il vettore alle iterazioni effettive
        break;
    end
    v1 = v1_new;
end

lambda1 = (v1' * A * v1) / (v1' * v1); % stima dell'autovalore

% Vettore iniziale v2 = [3; 10; 4]
v2 = [3; 10; 4];

% Vettore per memorizzare la norma delle differenze
norm_diff_v2 = zeros(max_iterazioni, 1);

% Iterazione per il metodo delle potenze con v2
for k = 1:max_iterazioni
    w = A * v2;
    v2_new = w / norm(w); % normalizzazione del vettore
    norm_diff_v2(k) = norm(v2_new - v2); % memorizza la norma della differenza
    if norm_diff_v2(k) < tolleranza
        norm_diff_v2 = norm_diff_v2(1:k); % riduce il vettore alle iterazioni effettive
        break;
    end
    v2 = v2_new;
end

lambda2 = (v2' * A * v2) / (v2' * v2); % stima dell'autovalore

% Stampa dei risultati del metodo delle potenze
disp('Metodo delle potenze con vettore iniziale [1; 1; 1]:');
disp(['Autovalore approssimato: ', num2str(lambda1)]);
disp(['Numero di iterazioni: ', num2str(length(norm_diff_v1))]);
disp('Autovettore approssimato:');
disp(v1);

disp('Metodo delle potenze con vettore iniziale [3; 10; 4]:');
disp(['Autovalore approssimato: ', num2str(lambda2)]);
disp(['Numero di iterazioni: ', num2str(length(norm_diff_v2))]);
disp('Autovettore approssimato:');
disp(v2);

% Analisi dettagliata del comportamento del metodo delle potenze con il vettore dato v2
disp('*** Analisi dettagliata del vettore [3; 10; 4] ***');

% Esperimenti con 50 e 200 iterazioni
iter_50 = 50;
iter_200 = 200;

% Metodo delle potenze con vettore v2 e 50 iterazioni
v2_50 = v2;
for k = 1:iter_50
    w = A * v2_50;
    v2_50 = w / norm(w);
end
lambda2_50 = (v2_50' * A * v2_50) / (v2_50' * v2_50);

% Metodo delle potenze con vettore v2 e 200 iterazioni
v2_200 = v2;
for k = 1:iter_200
    w = A * v2_200;
    v2_200 = w / norm(w);
end
lambda2_200 = (v2_200' * A * v2_200) / (v2_200' * v2_200);

% Stampa dei risultati per 50 e 200 iterazioni
disp(['Autovalore approssimato con 50 iterazioni: ', num2str(lambda2_50)]);
disp('Autovettore approssimato:');
disp(v2_50);

disp(['Autovalore approssimato con 200 iterazioni: ', num2str(lambda2_200)]);
disp('Autovettore approssimato:');
disp(v2_200);

% Grafico di confronto tra i due vettori v1 e v2 per il metodo delle potenze
fig1 = figure;
hold on;
plot(1:length(norm_diff_v1), norm_diff_v1, '-o', 'DisplayName', 'Metodo delle Potenze (v1 = [1; 1; 1])');
plot(1:length(norm_diff_v2), norm_diff_v2, '-x', 'DisplayName', 'Metodo delle Potenze (v2 = [3; 10; 4])');
xlabel('Numero di Iterazioni');
ylabel('Norma della differenza tra vettori consecutivi');
title('Convergenza del Metodo delle Potenze');
legend show;
grid on;
hold off;

saveas(fig1, '3.3_metodoPotenze.png');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PUNTO B (Metodo delle potenze inverse)

% PUNTO B1: Metodo delle potenze inverse con vettore iniziale [1; 1; 1] e p = 7
v_inv1 = [1; 1; 1]; % vettore iniziale
p1 = 7; % Stima iniziale dell'autovalore

% Vettore per memorizzare la norma delle differenze
norm_diff_inv1 = zeros(max_iterazioni, 1);

% Inizializzazione
lpr1 = p1;
w1 = v_inv1 / norm(v_inv1); % Normalizza il vettore iniziale
[L1, U1, P1] = lu(A - p1 * eye(size(A))); % Decomposizione LU della matrice (A - pI)

% Ciclo per il metodo delle potenze inverse
for m1 = 1:max_iterazioni
    v_inv1 = U1 \ (L1 \ (P1 * w1)); % Risoluzione del sistema lineare
    mu_max1 = w1' * v_inv1; % Calcola mu_max
    
    if abs(mu_max1) < tolleranza % Controllo per evitare divisioni per zero o numeri molto piccoli
        disp('Errore: mu_max1 è troppo piccolo, rischio di instabilità numerica.');
        lp1 = NaN;
        break;
    end
    
    lp1 = p1 + 1/mu_max1; % Aggiorna l'autovalore approssimato
    w1 = v_inv1 / norm(v_inv1); % Normalizza il nuovo vettore

    norm_diff_inv1(m1) = norm(v_inv1 - w1); % Memorizza la norma della differenza

    % Controllo della convergenza
    if abs(lp1 - lpr1) < tolleranza * abs(lp1)
        norm_diff_inv1 = norm_diff_inv1(1:m1); % riduce il vettore alle iterazioni effettive
        break;
    else
        lpr1 = lp1;
    end
end

% Visualizzazione dei risultati per il primo caso
if isnan(lp1)
    disp('Metodo delle potenze inverse (v = [1; 1; 1], p = 7): Calcolo fallito.');
else
    disp(['Metodo delle potenze inverse (v = [1; 1; 1], p = 7): Autovalore approssimato: ', num2str(lp1)]);
    disp(['Numero di iterazioni: ', num2str(m1)]);
    disp('Autovettore approssimato:');
    disp(w1);
end

% PUNTO B2: Metodo delle potenze inverse con vettore iniziale [3; 10; 4] e p = 2
v_inv2 = [3; 10; 4]; % vettore iniziale
p2 = 2; % Stima iniziale dell'autovalore

% Vettore per memorizzare la norma delle differenze
norm_diff_inv2 = zeros(max_iterazioni, 1);

% Inizializzazione
lpr2 = p2;
w2 = v_inv2 / norm(v_inv2); % Normalizza il vettore iniziale
[L2, U2, P2] = lu(A - p2 * eye(size(A))); % Decomposizione LU della matrice (A - pI)

% Ciclo per il metodo delle potenze inverse
for m2 = 1:max_iterazioni
    v_inv2 = U2 \ (L2 \ (P2 * w2)); % Risoluzione del sistema lineare
    mu_max2 = w2' * v_inv2; % Calcola mu_max
    
    if abs(mu_max2) < tolleranza % Controllo per evitare divisioni per zero o numeri molto piccoli
        disp('Errore: mu_max2 è troppo piccolo, rischio di instabilità numerica.');
        lp2 = NaN;
        break;
    end
    
    lp2 = p2 + 1/mu_max2; % Aggiorna l'autovalore approssimato
    w2 = v_inv2 / norm(v_inv2); % Normalizza il nuovo vettore

    norm_diff_inv2(m2) = norm(v_inv2 - w2); % Memorizza la norma della differenza

    % Controllo della convergenza
    if abs(lp2 - lpr2) < tolleranza * abs(lp2)
        norm_diff_inv2 = norm_diff_inv2(1:m2); % riduce il vettore alle iterazioni effettive
        break;
    else
        lpr2 = lp2;
    end
end

% Visualizzazione dei risultati per il secondo caso
if isnan(lp2)
    disp('Metodo delle potenze inverse (v = [3; 10; 4], p = 2): Calcolo fallito.');
else
    disp(['Metodo delle potenze inverse (v = [3; 10; 4], p = 2): Autovalore approssimato: ', num2str(lp2)]);
    disp(['Numero di iterazioni: ', num2str(m2)]);
    disp('Autovettore approssimato:');
    disp(w2);
end


% Grafico di confronto tra i due metodi delle potenze inverse
fig2 = figure;
hold on;
plot(1:length(norm_diff_inv1), norm_diff_inv1, '-o', 'DisplayName', 'Metodo delle Potenze Inverse (v = [1; 1; 1], p = 7)');
plot(1:length(norm_diff_inv2), norm_diff_inv2, '-x', 'DisplayName', 'Metodo delle Potenze Inverse (v = [3; 10; 4], p = 2)');
xlabel('Numero di Iterazioni');
ylabel('Norma della differenza tra vettori consecutivi');
title('Convergenza del Metodo delle Potenze Inverse');
legend show;
grid on;
hold off;

saveas(fig2, '3.3_metodoPotenzeInverse.png');

% Calcolo della velocità di convergenza per il metodo delle potenze con v1 e v2
convergenza_v1 = zeros(length(norm_diff_v1) - 1, 1);
for i = 2:length(norm_diff_v1)
    convergenza_v1(i-1) = norm_diff_v1(i) / norm_diff_v1(i-1);
end

convergenza_v2 = zeros(length(norm_diff_v2) - 1, 1);
for i = 2:length(norm_diff_v2)
    convergenza_v2(i-1) = norm_diff_v2(i) / norm_diff_v2(i-1);
end

% Calcolo della velocità di convergenza per il metodo delle potenze inverse con v1 e v2
convergenza_inv1 = zeros(length(norm_diff_inv1) - 1, 1);
for i = 2:length(norm_diff_inv1)
    convergenza_inv1(i-1) = norm_diff_inv1(i) / norm_diff_inv1(i-1);
end

convergenza_inv2 = zeros(length(norm_diff_inv2) - 1, 1);
for i = 2:length(norm_diff_inv2)
    convergenza_inv2(i-1) = norm_diff_inv2(i) / norm_diff_inv2(i-1);
end

% Visualizzazione dei risultati totali
disp('Velocità di convergenza per il metodo delle potenze (v1):');
disp(convergenza_v1);

disp('Velocità di convergenza per il metodo delle potenze (v2):');
disp(convergenza_v2);

disp('Velocità di convergenza per il metodo delle potenze inverse per p = 7 su (v1):');
disp(convergenza_inv1);

disp('Velocità di convergenza per il metodo delle potenze inverse per p = 2 su (v2):');
disp(convergenza_inv2);

diary off;