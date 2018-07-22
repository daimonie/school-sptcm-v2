function H2=getHeisenberg(J)
    % Constructing the Hamiltonian in basis {up,down}
    Sz = diag([1/2,-1/2]);
    Sp = [0 1;0 0];
    Sm = Sp';

    % Two site Hamiltonian:
    H2 = J*( kron(Sz,Sz) + 0.5*(kron(Sm,Sp) + kron(Sp,Sm)));
end