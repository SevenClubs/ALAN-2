%% main script

% creazione matrici

%1. matrice esercitazione 1
A = [3 1 -1 0; 
      0 7 -3 0;
      0 -3 9 -2;
      0 0 4 -10];
%disp(A);

%2. matrice tridiagonale

n = 10*(2+1)+2;
a = 2;
b = -1;

A = a*diag(ones(n,1)) + b*diag(ones(n-1,1),1) +b*diag(ones(n-1,1),-1);
%disp(A);

% 3. matrice 6x6

m = 6;

A = [1:m;
      1:m;
      1:m;
      1:m;
      1:m;
      1:m];
     
A = A/m;
A = transpose(A);
%disp(A);

% gli indici in Matlab partono da 1...
for i = 1:m
%     % eleva all'i-esima potenza la i-esima colonna di A
%     % le operazioni elemento per elemento si effettuano aggiungendo un punto
%     % prima dell'operatore
     A(:, i)=A(:,i).^(i-1);
end
disp(A);
% 
T = gauss_method(A);
% 
% % stampa la matrice
% 
disp(T);
