function l=edge_length(m_coord, edge_all)
    l = zeros(1, length(edge_all));

    for ii=1:length(edge_all);
        r_i = edge_all(ii, 1);
        r_j = edge_all(ii, 2);
        l(ii) = norm(m_coord(r_j,:)-m_coord(r_i, :));
    end
end