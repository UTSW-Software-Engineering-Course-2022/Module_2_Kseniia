function [X_idx, f_of_l] = X_idx_and_f_of_l(Fi, l)
    %read f(l) from Fi data structure
        dr = Fi.rn(2)-Fi.rn(1);
        r_min = Fi.rn(1);
        f_of_l = zeros(1, length(l));
        X_idx = zeros(1, length(l));
        n_shift = floor(r_min/dr+0.5);
        for ii=1:length(l);
            X_idx(ii) = floor(l(1, ii)/dr+0.5)-n_shift;
            f_of_l(ii) = Fi.fn(X_idx(ii));
        end
end