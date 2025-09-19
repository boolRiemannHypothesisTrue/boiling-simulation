function val = grad_T(QT,QbT1,QbT2,fields,params,phase,dx,dy)

        d1 = QT_dQT_dT(QT, fields,params,phase,dx,dy);

        d2 = Qb_T1_dQ_dT(QbT1);

        d3 = Qb_T1_dQ_dT(QbT2);

        val = d1 + d2 + d3;

end
