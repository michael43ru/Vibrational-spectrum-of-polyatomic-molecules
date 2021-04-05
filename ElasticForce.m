function F_k = ElasticForce(i, xyz, k)
F_k = [0;0;0];
for m = 1:size(k, 2)
    if (k(1,m) == i || k(2,m) == i)
        if (k(1,m) == i)
            j = k(2,m);
        end
        if (k(2,m) == i)
            j = k(1, m);
        end
        r = xyz(:,i) - xyz(:, j);
        F_k = F_k + k(3) * (-r)/norm(r) * (norm(r) - k(4));
    end
end
end