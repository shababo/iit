J = zeros(8,8);
J(1,2) = 1;
J(1,5) = 1;
J(2,3) = 1;
J(3,5) = 1;
J(4,3) = 1;
J(4,7) = 1;
J(4,8) = 1;
J(5,6) = 1;
J(6,8) = 1;
J(7,2) = 1;
J(7,6) = 1;
J(8,7) = 1;
J_BT = J/2;