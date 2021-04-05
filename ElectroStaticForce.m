function F_e = ElectroStaticForce(i, xyz, q)
F_e = [0;0;0];
%k = 1/(4*pi*8.85*10^(-12)*10^(-30)) * (1.6 * 10^(-19)) ^ 2 /((1.66 * 10 ^(-27)) * ((10^(-10)) ^ 2));
k = 1/(4*pi*8.85*10^(-12)*10^(-30)) * (1.6 * 10^(-19)) ^ 2;
for j = 1:size(xyz, 2)
    if(i ~= j)
        F_e = F_e + k*((xyz(:,i) - xyz(:,j)) * q(i) * q(j)) / ((xyz(1, i) - xyz(1, j)) ^ 2 + (xyz(2,i) - xyz(2,j))^2 + (xyz(3,i) - xyz(3,j)) ^ 2) ^ (3/2);
    end
end
end