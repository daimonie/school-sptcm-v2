function HL=getHL(H2,L)
    % get an L-site Hamiltonian with 2-site Hamiltonian as input
    H2 = sparse(H2);
    HL = H2;
    Id2 = speye(2); %sparse identity
    for k=3:L
        % construct identity operator on k-2 sites
        NId = speye(2^(k-2));
        HL = kron(HL,Id2)+kron(NId,H2);
    end
end
