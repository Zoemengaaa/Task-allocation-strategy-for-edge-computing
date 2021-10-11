                  function Dag = matrixstructure(x_n)

z_n = (x_n - 1)*x_n;

Dag = zeros((x_n - 1)*x_n, x_n+z_n);

for i =1:x_n
    Dag((x_n-1)*i-(x_n-2):(x_n-1)*i, (x_n+1)*i-x_n) = 1;
    
end

for i = 1:x_n-1
    Dag(x_n*i-x_n+1:x_n*i,x_n*i-x_n+1+i:x_n*i+i) = eye(x_n);
    
end
end