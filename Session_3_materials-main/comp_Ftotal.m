function Ftotal=comp_Ftotal(m, m_coord, f_of_l, l)
    Ftotal = zeros(length(m_coord), 3);
        for vortex=1:length(m_coord);
            ind = (m.var.edge_all==vortex);
            ind = sum(ind, 2);

            sublist = m.var.edge_all(logical(ind), :);
            sublist_f_l = f_of_l(logical(ind));
            sublist_l = l(logical(ind));
            Fsublist = zeros(length(sublist), 3);
            for sub_r=1:length(sublist);
                f_l = sublist_f_l(sub_r);
                i = sublist(sub_r, 1);
                j = sublist(sub_r, 2);
                ll = sublist(sub_r);
                dir_ij = (m_coord(j, :) - m_coord(i, :))/ll;

                if vortex==i;
                    F = f_l*(-1*dir_ij); %for vortex = i
                else
                    F = f_l*dir_ij; %for vortex = j
                end

                Fsublist(sub_r, :) = F;
            end

            Ftotal(vortex, :) = sum(Fsublist, 1);
        end
end