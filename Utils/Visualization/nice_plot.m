function nice_plot(a,X,Y) 
    figure
    
    pcolor(X, Y, a);
    shading interp;
    colorbar;
    colormap jet;
    
    xlabel('x (м)');
    ylabel('y (м)');
    axis equal tight
    
    
    norm(a)

end