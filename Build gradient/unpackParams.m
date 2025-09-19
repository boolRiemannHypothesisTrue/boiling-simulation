function [ux, uy, T, p] = unpackParams(qVec, sz)
    % sz = [Nx, Ny] — размер сетки
    Nx = sz(1);
    Ny = sz(2);
    nGrid = Nx * Ny;

    ux = reshape(qVec(1:nGrid), Nx, Ny);
    uy = reshape(qVec(nGrid+1:2*nGrid), Nx, Ny);
    T  = reshape(qVec(2*nGrid+1:3*nGrid), Nx, Ny);
    p  = reshape(qVec(3*nGrid+1:4*nGrid), Nx, Ny);
end
