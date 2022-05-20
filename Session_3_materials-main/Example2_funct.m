%% requirements:
% 1. gcc for compiling c codes. Load before starting Matlab:
%    module load gcc/6.3.0
% 2. Use Matlab 2020a
%% Setup the directory where the membrane object is located and add the directory to Matlab's function pool 
%dir_mod = '/home2/s171152/codes/matlab/mine/git/memCompCourse/memcompcourse';
dir_mod = '/archive/course/SWE22/train15/Downloads/Session_3_materials-main/';
addpath(dir_mod);
%--------------------------------------------------------------------------
% create 'unit' u using the unit module, and 'membrane' m using the membrane module 
u=ComUnit('erg',ComUnit.nm_to_cm(1000),300,ComUnit.kBT_to_erg(10,300)); 
m=ModMembrane(2,'unit',u);

m.pm.Vdh.V0=0.02
[Fi] = Finternal(m, 'plot_or_not', true)
m_coord = m.var.coord;
n_iter = 2
stds = zeros(n_iter, 1);
for iter=1:n_iter
    %calculate length of each edge
    l = edge_length(m_coord, m.var.edge_all);

    X_idx, f_of_l = X_idx_and_f_of_l(Fi, l);

    %compute F total
    Ftotal=comp_Ftotal(m, m_coord, f_of_l, l);
    
    %compute delta t final
    delta_t_final = comp_delta_final(X_idx, m, l, Fi, m_coord, Ftotal)
    
    %compute random force
    k = 0;
    F_random = k*randn(length(m.var.coord), 3);
    
    Ftotal = Ftotal+F_random;
    m_coord = m_coord+m.pm.mu*Ftotal*delta_t_final;
    m.var.coord = m_coord;
    stds(iter) = std(l);
    
end


%%
plot(1:n_iter, stds)
%% Plot the membranel 'm'.
fig=figure;
subplot(1,2,1);
plot(m,'f',fig);
subplot(1,2,2);
col=rand(m.var.n_coord,3);
plot(m,'f',fig,'col',col,'col_min',0,'col_max',1,'colBar',true);

