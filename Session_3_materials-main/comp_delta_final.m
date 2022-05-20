function delta_t_final = comp_delta_final(X_idx, m, l, Fi, m_coord, Ftotal)
    delt_t_all = zeros(1, length(X_idx));
        for ii=1:length(X_idx);
            i = m.var.edge_all(ii, 1);
            j = m.var.edge_all(ii, 2);
            l_pl = 0.5*(Fi.rg(Fi.in(X_idx(ii))+1)+Fi.ig(Fi.in(X_idx(ii))));
            l_mi = 0.5*(Fi.rg(Fi.in(X_idx(ii))-2)+Fi.rg(Fi.in(X_idx(ii))-1));

            u = (norm(l_pl)^2)-(l(ii)^2);
            dot_prod = dot((m_coord(j, :)-m_coord(i, :)), (Ftotal(j, :)-Ftotal(i, :)));
            d = 2*m.pm.mu*dot_prod;
            delt_t_pl = u/d;

            u = (norm(l_mi)^2)-(l(ii)^2);
            dot_prod = dot((m_coord(i, :)-m_coord(j, :)), (Ftotal(i, :)-Ftotal(j, :)));
            d = 2*m.pm.mu*dot_prod;
            delt_t_mi = u/d;

            delt_t = 0;
            if delt_t_pl > 0
                delt_t = delt_t_pl;
            else
                delt_t = delt_t_mi;
            end
            delt_t_all(ii) = delt_t;
        end
        delta_t_final = min(delt_t_all);
end