function    F = discreteF(F0)


[SF,I] = sort(-F0);
SF = -SF;

[m,n] = size(F0);
F = zeros(m,n);

for i = 1:n
   F(I(1,i),i) = 1;
end

end