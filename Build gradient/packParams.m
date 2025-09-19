function qVec = packParams(ux, uy, T, p)
    % ux, uy, T, p — матрицы одинакового размера [Nx×Ny]
    qVec = [ux(:); uy(:); T(:); p(:)];
end
