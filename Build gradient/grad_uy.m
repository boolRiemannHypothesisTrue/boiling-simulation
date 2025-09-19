function val = grad_uy(Qu,QT,Qx,Qy,Qbuy1,Qbuy2,fields,params,phase,dx,dy)

        d1 = Qu_dQu_duy(params,phase,Qu,dx);

        d2 = Qx_dQx_duy(Qx,fields,params,phase,dx,dy);

        d3 = Qy_dQy_duy(Qy,fields,params,phase,dx,dy);

        d4 = QT_dQT_duy(QT, fields,params,phase,dx,dy);

        d5 = Qb_uy1_dQ_duy(Qbuy1);

        d6 = Qb_uy2_dQ_duy(Qbuy2);

        val = d1 + d2 + d3 + d4 +d5 + d6;
        
end
