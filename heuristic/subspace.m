function Un = subspace(Y, L)
    Cov = (Y * Y') / L;
    [eigenvectors,eigenvalues] = eig(Cov);
    [~, idx] = sort(diag(eigenvalues),'descend');
    eigenvectors = eigenvectors(:, idx);
    Un = eigenvectors(:, 2:end);
end