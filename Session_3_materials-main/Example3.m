%% requirements:
% 1. gcc for compiling c codes. Load before starting Matlab:
%    module load gcc/6.3.0
% 2. Use Matlab 2020a
%% Setup the directory where the membrane object is located and add the directory to Matlab's function pool 
%dir_mod = '/home2/s171152/codes/matlab/mine/git/memCompCourse/memcompcourse';
dir_mod = '/archive/course/SWE22/train15/Downloads/Session_3_materials-main2/';
addpath(dir_mod);
%--------------------------------------------------------------------------
% create 'unit' u using the unit module, and 'membrane' m using the membrane module 
u=ComUnit('erg',ComUnit.nm_to_cm(1000),300,ComUnit.kBT_to_erg(10,300)); 
m=ModMembrane(2,'unit',u);
[Fi] = Finternal(m, 'plot_or_not', true);
m.pm.Vdh.V0 = 0.02;
m_coord = m.var.coord;
n_iter = 20;
kc = 100;
ks = 0.001; %100 makes surface irregular, 0.001 makes it flat? use 0.0001
kv = 100;
m.pm.k_c = kc;
r=mean(sqrt(sum(m.var.coord(:,1).^2+m.var.coord(:,2).^2+m.var.coord(:,3).^2,2)));
V0 = 0.6*4/3*r^3;
Vtotal = kv*((sum(Volume(m))-V0)^2)/V0;
A0 = 4*pi*r^2;
Atotal = ks*((sum(Area(m))-A0)^2)/A0;
Htotal = sum(Helfrich(m));

stds = zeros(n_iter, 1);
for iter=1:n_iter
    %calculate l for each edge
    l = zeros(1, length(m.var.edge_all));
    
    for ii=1:length(m.var.edge_all);
        r_i = m.var.edge_all(ii, 1);
        r_j = m.var.edge_all(ii, 2);
        l(ii) = norm(m_coord(r_j,:)-m_coord(r_i, :));
    end

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

    %compute F for each vortex of each edge and store in Fs
    
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
    
    
    %compute Fb
    delta_x = 0.000001;
    delta_y = 0.000001;
    delta_z = 0.000001;
          
    Fb= zeros(162, 3);
    for i=1:162;
        %compute x axis
        m.var.coord(i,1) = m.var.coord(i,1)+delta_x;
        H_new = sum(Helfrich(m));
        Fb_i = -1*(H_new-Htotal)/delta_x;
        Fb(i,1) = Fb_i;
        m.var.coord = m_coord;
        
        %compute y axis
        m.var.coord(i,2) = m.var.coord(i,2)+delta_y;
        H_new = sum(Helfrich(m));
        Fb_i = -1*(H_new-Htotal)/delta_y;
        Fb(i, 2) = Fb_i;
        m.var.coord = m_coord;
        
        %compute z axis
        m.var.coord(i,3) = m.var.coord(i,3)+delta_z;
        H_new = sum(Helfrich(m));
        Fb_i = -1*(H_new-Htotal)/delta_z;
        Fb(i,3) = Fb_i;
        m.var.coord = m_coord;
    end
    %quiver3(m.var.coord(:,1), m.var.coord(:,1), m.var.coord(:,3), Fb(:,1), Fb(:,2), Fb(:,3))
    
    
    %Fa = ks*((A-A0)^2)/A0
    Fa = zeros(162, 3);
    for i=1:162;
        %compute x axis
        m.var.coord(i,1) = m.var.coord(i,1)+delta_x;
        A=sum(Area(m));
        A_i = -1*(ks*((A-A0)^2)/A0-Atotal)/delta_x;
        Fa(i,1) = A_i;
        m.var.coord = m_coord;
        
        %compute y axis
        m.var.coord(i,2) = m.var.coord(i,2)+delta_y;
        A=sum(Area(m));
        A_i = -1*(ks*((A-A0)^2)/A0-Atotal)/delta_y;
        Fa(i,2) = A_i;
        m.var.coord = m_coord;
        
        %compute z axis
        m.var.coord(i,3) = m.var.coord(i,3)+delta_z;
        A=sum(Area(m));
        A_i = -1*(ks*((A-A0)^2)/A0-Atotal)/delta_z;
        Fa(i,3) = A_i;
        m.var.coord = m_coord;
    end
    %quiver3(m.var.coord(:,1), m.var.coord(:,1), m.var.coord(:,3), Fa(:,1), Fa(:,2), Fa(:,3))
    
    %compute Fv
    %Fv = kv*((V-V0)^2)/V0
    Fv = zeros(162,3);
    for i=1:162;
        m.var.coord(i,1) = m.var.coord(i,1)+delta_x;
        V=sum(Volume(m));
        V_i = -1*(kv*((V-V0)^2)/V0-Vtotal)/delta_x;
        Fv(i,1) = V_i;
        m.var.coord = m_coord;
        
        m.var.coord(i,2) = m.var.coord(i,2)+delta_y;
        V=sum(Volume(m));
        V_i = -1*(kv*((V-V0)^2)/V0-Vtotal)/delta_y;
        Fv(i,2) = V_i;
        m.var.coord = m_coord;
        
        m.var.coord(i,3) = m.var.coord(i,3)+delta_z;
        V=sum(Volume(m));
        V_i = -1*(kv*((V-V0)^2)/V0-Vtotal)/delta_z;
        Fv(i,3) = V_i;
        m.var.coord = m_coord;
    end
    %quiver3(m.var.coord(:,1), m.var.coord(:,1), m.var.coord(:,3), Fv(:,1), Fv(:,2), Fv(:,3))
    
    %add random force
    k = 0;
    F_random = k*randn(length(m.var.coord), 3);
    
    Ftotal = Ftotal+F_random+Fv+Fa+Fb;
    %Ftotal = Ftotal+F_random+Fb;
    
    %compute deltat_final
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
    delta_t_final = min(delt_t_all)
    
    
    m_coord = m_coord+m.pm.mu*Ftotal*delta_t_final;
    m.var.coord=m_coord;
    std(l)
    stds(iter) = std(l);
    
    rLim=(min(Fi.rg, max(Fi.rg)));
    
    [m] = RemeshCtrl(m, Fi,rLim);
    
end

m.var.coord = m_coord;



%%

quiver3(m.var.coord(:,1), m.var.coord(:,1), m.var.coord(:,3), Fb(:,1), Fb(:,2), Fb(:,3))

%%
m.var.coord(i,1):Fa(i,1)
%%
plot(1:n_iter, stds)
%%
%Plot the membranel 'm'. Note that Matlab autonatically recognize m is an 'object' and apply m's own plot function  
fig=figure;
subplot(1,2,1);
plot(m,'f',fig);
subplot(1,2,2);
col=rand(m.var.n_coord,3);
plot(m,'f',fig,'col',col,'col_min',0,'col_max',1,'colBar',true);
%--------------------------------------------------------------------------

