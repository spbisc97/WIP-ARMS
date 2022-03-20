function z = pend(t, y)
    k=0.1;
    z =  [y(2); 4.9050*sin(y(1))-k*y(2) ];%
end