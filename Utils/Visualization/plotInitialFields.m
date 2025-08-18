function plotInitialFields(X, Y, fields)
    figure;

    % u_x
    subplot(2,2,1);
    pcolor(X, Y, fields.u_x);
    shading interp;
    colorbar;
    colormap jet;

    title('Начальное поле u_x (м/c)');
    xlabel('x (м)');
    ylabel('y (м)');
    axis equal tight;

    % u_y
    subplot(2,2,2);
    pcolor(X, Y, fields.u_y);
    shading interp;
    colorbar;
    colormap jet;

    title('Начальное поле u_y (м/с)');
    xlabel('x (м)');
    ylabel('y (м)');
    axis equal tight;

    % T
    subplot(2,2,3);
    pcolor(X, Y, fields.T);
    shading interp;
    colorbar;
    colormap jet;

    title('Начальное поле температуры T (К)');
    xlabel('x (м)');
    ylabel('y (м)');
    axis equal tight;

    % p
    subplot(2,2,4);
    pcolor(X, Y, fields.P);
    shading interp;
    colorbar;
    colormap jet;
 
    title('Начальное поле давления p (Па)');
    xlabel('x (м)');
    ylabel('y (м)');
    axis equal tight;
end
