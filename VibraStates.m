function [Fr,Dr] = VibraStates(xyz,q,m,k)
dl = (abs(xyz(1,1) - xyz(1,2)) + abs(xyz(2,1) - xyz(2,2)) + abs(xyz(3,1) - xyz(3,2))) / 10000;
dx = [dl; 0; 0];
dy = [0; dl; 0];
dz = [0; 0; dl];
xyz_1 = xyz;
xyz_2 = xyz;
N = size(xyz, 2);
P = zeros(3*N);

cond = 1;
F_e = [0; 0; 0];
F_k = [0; 0; 0];
F = [0; 0; 0];
for i=1:N
    F_k = ElasticForce(i, xyz, k);
    F_e = ElectroStaticForce(i, xyz, q);
    F = F_e + F_k;
    if (abs((F'*F) / (F_k'*F_k)) > 0.001)
        cond = 0;
    end
end

if (cond == 0)
    disp('Нет равновесия');
else
    disp('Есть равновесие');

for i=1:3*N
    for j = 1:3*N
        xyz_1 = xyz;
        xyz_2 = xyz;
        xyz_2(i) = xyz(i) + dl;
        p1 = Potential((i-1 - mod(i-1, 3)) / 3 + 1, xyz_1, q, m, k);
        p2 = Potential((i-1 - mod(i-1, 3)) / 3 + 1, xyz_2, q, m, k);
        %(i-1 - mod(i-1, 3)) / 3 + 1
        dp_1 = (p2 - p1) / dl;
        xyz_1(j) = xyz_1(j) + dl;
        xyz_2(j) = xyz_2(j) + dl;
        p1 = Potential((i-1 - mod(i-1, 3)) / 3 + 1, xyz_1, q, m, k);
        p2 = Potential((i-1 - mod(i-1, 3)) / 3 + 1, xyz_2, q, m, k);
        dp_2 = (p2 - p1) / dl;
        ddp = (dp_2 - dp_1) / dl;
        P(i, j) = ddp;
    end
end
P = P./2;
Dr = zeros(3*N, 3*N); % 3*N-6
Fr = zeros(1, 3*N); % 3*N-6
M = zeros(3*N);
for i=1:N
    M(3*i - 2, 3*i - 2) = m(i);
    M(3*i - 1, 3*i - 1) = m(i);
    M(3*i, 3*i) = m(i);
end
x = sym('x');
X = solve(det(M * x + P) == 0);
X = double(X);
for i=1:length(X)
    if (abs(X(i) / X(3*N)) < 10 ^ (-7)) 
        X(i) = 0;
    end
end
j = 1;
Y = zeros(3*N-6, 1);
for i=1:length(X) % 3*N
    if (X(i) > 0) 
        Y(j) = X(i);
        j = j + 1;
    end
end
Fr = (sqrt(Y)).';
for i=1:length(Fr) % 3*N-6
    if (Fr(i) > 0)
        [A, L] = eig(M * Fr(i) ^ 2 + P);
        min = abs(L(1,1));
        for j = 2:length(L)
            if (abs(L(j, j)) < min)
                min = abs(L(j, j));
            end
        end
        Dr(:,i) = A(:,j);
    end
end
end

