function twist = vee(g)

g = formula(g);
R = g(1:3,1:3);
v = g(1:3,4);


twist = [v;R(3,2);R(1,3);R(2,1)];
end