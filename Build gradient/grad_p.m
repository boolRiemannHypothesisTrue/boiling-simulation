function val = grad_p(Qx,Qy,Qbp,dx,dy)

        d1 = Qx_dQx_dp_PLUS_Qy_dQy_dp(Qx, Qy, dx, dy);

        d2 = Qb_p_dQ_dp(Qbp);



        val = d1 + d2;

end
