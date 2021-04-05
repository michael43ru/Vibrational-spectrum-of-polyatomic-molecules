function P = Potential(i, xyz, q, m, k)
P = 0;
k_e = 1/(4*pi*8.85*10^(-12)) * (1.6 * 10^(-19)) ^ 2 /((1.66 * 10 ^(-27)) * ((10^(-10)) ^ 2));
for j = 1:size(xyz, 2)
    if(i ~= j)
        P = P - k_e * q(i) * q(j) / ((xyz(1, i) - xyz(1, j)) ^ 2 + (xyz(2,i) - xyz(2,j))^2 + (xyz(3,i) - xyz(3,j)) ^ 2);
    end
end

for n = 1:size(k, 2)
    if (k(1,n) == i || k(2,n) == i)
        if (k(1,n) == i)
            j = k(2,n);
        end
        if (k(2,n) == i)
            j = k(1, n);
        end
        r = ((xyz(1,i) - xyz(1,j)) ^ 2 + (xyz(2,i) + xyz(2,j)) ^ 2 + (xyz(3,i) + xyz(3,j)) ^ 2) ^ 0.5;
        P = P + k(3,n) * (r - k(4,n)) ^ 2 / 2;
    end
end
end