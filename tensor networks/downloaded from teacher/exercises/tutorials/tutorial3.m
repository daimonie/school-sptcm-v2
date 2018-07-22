%% Tutorial 3: Exact diagonalization
% Let's do a simple exact diagonalization of a Heisenberg model with open
% boundary conditions:
%%
% 
% $$H = J \sum_{\langle i,j \rangle} S_i S_j = H=J \sum_{\langle i,j\rangle} S^z_i S^z_j + \frac{1}{2} (S^+_i S^-_j + S^-_i  S^+_j )$$
% 

%% Constructing the Hamiltonian in basis {up,down}
Sz = diag([1/2,-1/2]);
Sp = [0 1;0 0];
Sm = Sp';
Id = eye(2);

% Two site Hamiltonian:
H2 = kron(Sz,Sz) + 0.5*(kron(Sm,Sp) + kron(Sp,Sm))

%% Diagonalize and get lowest energy state
[V,D] = eigs(H2,1,'SA');
disp(['2-site ground state energy:',num2str(D)]);
disp('Eigenvector (singlet):')
disp(V);


%% And now construct L-site Hamiltonian recursively
L = 10; %example
% Use sparse representation for Hamiltonian:
H2 = sparse(H2);
HL = H2;
Id2 = speye(2); %sparse identity
for k=3:L
    % construct identity operator on k-2 sites
    NId = speye(2^(k-2));
    HL = kron(HL,Id2)+kron(NId,H2);
end

% Diagonalize and get lowest energy state
[V,D] = eigs(HL,1,'SA');
disp(['10-site ground state energy:',num2str(D)]);