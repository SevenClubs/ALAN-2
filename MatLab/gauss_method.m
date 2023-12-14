% implementare una funzione per triangolarizzare una matrice con il metodo 
% di eliminazione di Gauss
% Parametri:
% input:  A, la matrice da triangolarizzare
% output: T, la matrice triangolarizzata

function T = gauss_method(A)
    
    % ottenere le dimensioni di A tramite il comando size
    [row, col] = size(A);
    
    % inizializzazione della matrice di output con una matrice di zeri
    
    T = zeros(row,col);
    disp(T);

    % gli indici in Matlab partono da 1...
    for j = 1:col-1
        %controllo se nella colonna ci sono solo zeri. se la colonna e
        %nulla, mi sposto sulla colonna a dx

        %cse ho superato il controllo, vado al passo 0
        if (A(j,j)=0)
            
        % elemento della diagonale di T
        
        %perno = ...;
        
        for i = j+1:row
                % calcola il moltiplicatore 
                %mult = ...;
                %T(i,:) = T(i,:) - mult*T(j, :);
        end
    end
%end

