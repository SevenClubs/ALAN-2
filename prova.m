clc
clear
close all

n = 5;
A = diag(ones(1, n-1), 1) + eye(n);
disp(A);
E = zeros(n);
E(n,1) = 2^(-n);
B = A + E;
VA = eig(A);
VB = eig(B);

matrixdistance = norm(A-B)/norm(A);
autovdistance = norm(VA-VB)/norm(VA);
disp(matrixdistance);
disp(autovdistance);

Ata = A'*A;
Btb = B'*B;
VAta = eig(Ata);
VBtb = eig(Btb);
matrixdistance = norm(Ata-Btb)/norm(Ata);
autovdistance = norm(VAta-VBtb)/norm(VAta);
