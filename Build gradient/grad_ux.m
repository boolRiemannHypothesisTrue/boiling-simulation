function val = grad_ux(Qu,QT,Qx,Qy,Qbux1,Qbux2,fields,params,phase,dx,dy)

        d1 = Qu_dQu_dux(params,phase,Qu,dx);

        d2 = Qx_dQx_dux(Qx,fields,params,phase,dx,dy);

        d3 = Qy_dQy_dux(Qy,fields,params,phase,dx,dy);

        d4 = QT_dQT_dux(QT, fields,params,phase,dx,dy);

        d5 = Qb_ux1_dQ_dux(Qbux1);

        d6 = Qb_ux2_dQ_dux(Qbux2);

        val = d1 + d2 + d3 + d4 +d5 + d6;
        
end
