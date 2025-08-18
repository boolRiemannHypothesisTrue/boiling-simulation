function fields = unpackFieldsOnly(x, M, N)
    total_points = M * N;

    ux_vec = x(1:total_points);
    uy_vec = x(total_points+1:2*total_points);
    T_vec  = x(2*total_points+1:3*total_points);
    P_vec  = x(3*total_points+1:4*total_points);

    % Восстановление полей из векторов
    fields.u_x = reshape(ux_vec, M, N);
    fields.u_y = reshape(uy_vec, M, N);
    fields.T   = reshape(T_vec, M, N);
    fields.P   = reshape(P_vec, M, N);
end
