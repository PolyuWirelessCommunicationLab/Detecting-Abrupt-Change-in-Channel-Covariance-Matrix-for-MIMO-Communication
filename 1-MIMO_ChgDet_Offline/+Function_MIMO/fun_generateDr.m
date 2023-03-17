function D = fun_generateDr(Num1,Num2)
D1 = zeros(Num1);
for p = 1:Num1
    for q = 1:Num1
            D1(p,q) = p-q;
    end
end

D2 = ones(Num2);

D = kron(D1,D2);