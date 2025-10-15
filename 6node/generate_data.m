close all;
clear all;
clc;

% === 1) RUTA de tu GAMS ===
gamsHome = 'C:\GAMS\50';  % <-- CAMBIA esto

[basedir,~,~] = fileparts(mfilename('fullpath'));
basedir = fullfile(basedir, 'export_csv');   % subcarpeta de exportación
if ~exist(basedir,'dir'); mkdir(basedir); end

% === Tus datos ===
[n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    load_factor,op_link_cost,congestion_coef_stations,...
    congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,eta,...
    a_max,candidasourcertes] = parameters_6node_network();

M = 1e4;
nreg = 10;
eps = 1e-3;
vals_regs = linspace(0.005,0.995,nreg-1);
[lin_coef,bord,b] = get_linearization(n,nreg,alt_time,alt_price,-0.2,-0.2,vals_regs);

candidates = zeros(n);
for i=1:n
    candidates(i,candidasourcertes{i}) = 1;
end

%%

% === Helpers para escribir CSV ===
write_matrix_csv = @(A, fn) writetable( ...
  array2table(A, ...
    'VariableNames', cellstr(compose('j%d', 1:size(A,2))), ...
    'RowNames',     cellstr(compose('i%d', (1:size(A,1))')) ), ...
  fullfile(basedir, fn), 'WriteRowNames', true );


write_vector_csv = @(v, fn, prefix) ( ...
    writetable( table( cellstr(prefix+string((1:numel(v)).')), v(:), ...
               'VariableNames', {'idx','value'} ), ...
               fullfile(basedir, fn)) ...
);

write_scalar_csv_append = @(name, val, fn) ( ...
    writetable( table( string(name), double(val), ...
               'VariableNames', {'name','value'} ), ...
               fullfile(basedir, fn), 'WriteMode','append') ...
);


write_gams_param_iii('./export_txt/lin_coef.txt', lin_coef);
write_gams_param_iii('./export_txt/b.txt', b);
write_gams_param_iii('./export_txt/bord.txt', bord);

% === 2D matrices ===
write_gams_param_ii('./export_txt/demand.txt', demand);
write_gams_param_ii('./export_txt/travel_time.txt', travel_time);
write_gams_param_ii('./export_txt/alt_time.txt', alt_time);
write_gams_param_ii('./export_txt/alt_price.txt', alt_price);
write_gams_param_ii('./export_txt/link_cost.txt', link_cost);
write_gams_param_ii('./export_txt/link_capacity_slope.txt', link_capacity_slope);
write_gams_param_ii('./export_txt/prices.txt', prices);
write_gams_param_ii('./export_txt/op_link_cost.txt', op_link_cost);
write_gams_param_ii('./export_txt/candidates.txt', candidates);
write_gams_param_ii('./export_txt/congestion_coefs_links.txt', congestion_coef_links);



% === 1D vectores ===
write_gams_param1d_full('./export_txt/station_cost.txt', station_cost);
write_gams_param1d_full('./export_txt/station_capacity_slope.txt', station_capacity_slope);
write_gams_param1d_full('./export_txt/congestion_coefs_stations.txt', congestion_coef_stations);

alfa = 0.5;
betas = 1:16;
budgets = 1e4.*linspace(1.5,4,16);
lam = 5;
%% run MIPREG


for bb=1:length(betas)
    eps = 1e-3;
    alfa = 0.5;
    lam = 5;
    beta = betas(bb);
    bud = budgets(bb);
    disp(bud);
    tic;
    [sprim,s,deltas,aprim,a,deltaa,f,fext,fij,comp_time] = compute_sim_MIP_entr(lam,beta,alfa,n,bud);
    [obj_val,pax_obj,op_obj] = get_obj_val(alfa, op_link_cost,...
            travel_time,prices,alt_time,alt_price,aprim,deltaa, ...
            sprim,deltas,fij,f,fext,demand);
       budget = get_budget(s,sprim,a,aprim,n,...
            station_cost,station_capacity_slope,link_cost,link_capacity_slope,lam);
       filename = sprintf('./6node_rebutal_MIPREG/beta=%d_lam=%d.mat',beta,lam);
        save(filename,'s','sprim','deltas', ...
        'a','aprim','deltaa','f','fext','fij','comp_time','budget', ...
        'pax_obj','op_obj','obj_val');
    
    disp(['bb = ',num2str(bb),', budget = ',num2str(budget),', obj_val = ',num2str(obj_val),', pax_obj = ',num2str(pax_obj), ...
        ', op_obj = ',num2str(op_obj),', nlinks = ', num2str(sum(sum(a > eps)))]);

end
%%

for bb=1:length(betas)
    beta = betas(bb);
    bud = budgets(bb);
    [sprim,s,deltas,aprim,a,deltaa,f,fext,fij] = compute_sim_MIP(lam,beta,alfa,n,bud);
end

%% Comparison known network

dif = zeros(length(betas),1);
opt_gaps = dif;

for bb=1:length(betas)
    disp(bb);
    beta = betas(bb);
    filename = sprintf('./6node_rebutal_MIPREG/beta=%d_lam=%d.mat',beta,lam);
    load(filename);    
    obj_val_MIPREG = obj_val;
    filename = sprintf('./6node_rebutal_MIP_10min/beta=%d_lam=%d.mat',beta,lam);
    load(filename);
    obj_val_MIP = obj_val;
    opt_gaps(bb) = 100*mipgap;

    dif(bb) = 100*(obj_val_MIPREG-obj_val_MIP)/obj_val_MIP;


end


%% random networks
betas = 1:8;
budgets_aval = 1e4.*linspace(1.5,4,8);

lams = [5];
n_runs = 10;
budgets = zeros(length(betas),length(lams),n_runs);


for ll = 1:length(lams)
   lam = lams(ll);
   for bb=1:length(betas)
       beta = betas(bb);
       budget_aval = budgets_aval(bb);
       for rr=1:n_runs
           rng(rr);
            [n,link_cost,station_cost,link_capacity_slope,...
            station_capacity_slope,demand,prices,...
            load_factor,op_link_cost,travel_time,alt_time,alt_price,a_nom,tau,eta,...
            a_max,candidates] = parameters_6node_network_rand();
           tic;
           [sprim,s,deltas,aprim,a,deltaa,f,fext,fij,comp_time] = compute_sim_MIP_entr(lam,beta,alfa,n,budget_aval);
           [obj_val,pax_obj,op_obj] = get_obj_val(alfa, op_link_cost,...
                travel_time,prices,alt_time,alt_price,aprim,deltaa, ...
                sprim,deltas,fij,f,fext,demand);

           budget = get_budget(s,sprim,a,aprim,n,...
                station_cost,station_capacity_slope,link_cost,link_capacity_slope,lam);
           budgets(bb,ll) = budget;
           filename = sprintf('./6node_rebutal_MIPREG/beta=%d_lam=%d_run=%d.mat',beta,lam,rr);
            save(filename,'s','sprim','deltas', ...
            'a','aprim','deltaa','f','fext','fij','comp_time','budget', ...
            'pax_obj','op_obj','obj_val');
            % now, run MIP
            [sprim,s,deltas,aprim,a,deltaa,f,fext,fij] = compute_sim_MIP_rand(lam,beta,alfa,n,budget,rr, ...
    op_link_cost, travel_time,prices,alt_time,alt_price,demand,station_cost,link_cost,...
    station_capacity_slope,link_capacity_slope);

       end
      
   end

end
%%
n_runs = 10;betas = 1:8;lams = [5];
dif = zeros(length(betas),n_runs);
gaps = dif;
times_MIP = dif;
budgets_aval = 1e4.*linspace(1.5,4,8);
times_MIPREG = dif;
links_MIP = dif;
links_MIPREG = dif;
gaps_10 = gaps;
links_MIP_10 = dif;
dif_10 = dif;
times_MIP_10 = dif;
lam = 5;

for bb=1:length(betas)
    beta = betas(bb);
    for rr=1:n_runs
        filename = sprintf('./6node_rebutal_MIPREG/beta=%d_lam=%d_run=%d.mat',beta,lam,rr);
        load(filename);
        times_MIPREG(bb,rr) = comp_time;
        links_MIPREG(bb,rr) = sum(sum(a>1e-2));
        obj_val_MIPREG = obj_val;
        filename = sprintf('./6node_rebutal_MIP_30min/beta=%d_lam=%d_run=%d.mat',beta,lam,rr);
        load(filename);
        obj_val_MIP = obj_val;
        gaps(bb,rr) = 100*mipgap;
        links_MIP(bb,rr) = sum(sum(a>1e-2));
        times_MIP(bb,rr) = comp_time;
        dif(bb,rr) = 100*(obj_val_MIPREG-obj_val_MIP)/obj_val_MIP;

        filename = sprintf('./6node_rebutal_MIP_10min/beta=%d_lam=%d_run=%d.mat',beta,lam,rr);
        load(filename);
        obj_val_MIP_10 = obj_val;
        times_MIP_10(bb,rr) = comp_time;
        gaps_10(bb,rr) = 100*mipgap;
        links_MIP_10(bb,rr) = sum(sum(a>1e-2));
        times_MIP_10(bb,rr) = comp_time;
        dif_10(bb,rr) = 100*(obj_val_MIPREG-obj_val_MIP_10)/obj_val_MIP_10;
    end

end



dif_10_vals = mean(dif_10,2);
dif_vals = mean(dif,2);
gaps_10_vals = mean(gaps_10,2);
gaps_vals = mean(gaps,2);

gaps_10_stds = std(gaps_10,0,2);
gaps_stds = std(gaps,0,2);

dif_10_stds = std(dif_10,0,2);
dif_stds = std(dif,0,2);


times_MIP_vals = mean(times_MIP,2);
times_MIP_10_vals = mean(times_MIP_10,2);

times_MIP_std = std(times_MIP,0,2);
times_MIP_10_std = std(times_MIP_10,0,2);

times_MIPREG_vals = mean(times_MIPREG,2);
times_MIPREG_std = std(times_MIPREG,0,2);


for bb=1:length(betas)
    beta = betas(bb);
    avail_bud = budgets_aval(bb);
    dif_10_val = dif_10_vals(bb);
    dif_val = dif_vals(bb);
    dif_10_std = dif_10_stds(bb);
    dif_std = dif_stds(bb);
    disp([sprintf('%.1e', avail_bud),'&',sprintf('%.2f $\\pm$ %.1f',dif_10_val,dif_10_std),'&',...
        sprintf('%.2f $\\pm$ %.1f',dif_val,dif_std),...
        '&',...
        sprintf('%.1f $\\pm$ %.1f',gaps_10_vals(bb),gaps_10_stds(bb)), ...
        '&',sprintf('%.1f $\\pm$ %.1f',gaps_vals(bb),gaps_stds(bb)),'&', ...%sprintf('%.2f',-opt_gap_MIP_8h),'&',...
        sprintf('%.1e $\\pm$ %.1e',times_MIP_10_vals(bb),times_MIP_10_std(bb)),'&',...sprintf('%.2e',t_comp_MIP_8h),...'&',
        sprintf('%.1e $\\pm$ %.1e',times_MIP_vals(bb),times_MIP_std(bb)),'&',...
        sprintf('%.1e $\\pm$ %.1e',times_MIPREG_vals(bb),times_MIPREG_std(bb)),'\\ \hline']);

end







%%


close all;
figure('Position', [100, 100, 1000, 800]);

azul_col = [0 0.4470 0.7410]; %icvx
rojo_col = [0.8500 0.3250 0.0980]; %mip30
naranja_col = [0.9290 0.6940 0.1250]; %MIP10
verde_col = [0.4660 0.6740 0.1880]; %mipreg

subplot(221);
% A y B: nObs x nX  (cada columna = un valor de X)
[nObs, nX] = size(dif');

posA = (1:nX) - 0.15;   % desplazamiento izquierda para A
posB = (1:nX) + 0.15;   % desplazamiento derecha para B

cla; hold on
boxplot(dif', 'positions', posA, 'Widths', 0.25, 'Colors',[0 0 0], 'Symbol','');
boxplot(dif_10', 'positions', posB, 'Widths', 0.25, 'Colors',[0.5 0.5 0.5], 'Symbol','');

xlim([0.5, nX+0.5]); xticks(1:nX)
eur =['[',char(8364),']'];
xl = xlabel(['Budget'],'interpreter','latex');
yl = ylabel('Diff [\%]','Interpreter','latex');
set(gca, 'FontSize', 9);
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'XTick', 1:8, 'XTickLabel', budgets_aval);

% leyenda “dummy”
plot(nan,nan,'Color',[0 0 0],'LineWidth',2); plot(nan,nan,'Color',[0.5 0.5 0.5],'LineWidth',2);
legend({'MIPREG w.r.t. MIP 30 min','MIPREG w.r.t. MIP 10 min'},'Location','best','Interpreter','latex','FontSize',9); hold off
xlim([0.7 8.3]);


subplot(222);
x = 1:size(times_MIP,1);

% Estadísticos
mu_MIP = mean(times_MIP,2);               
sigma_MIP = std(times_MIP,0,2);   


mu_MIP_10 = mean(times_MIP_10,2);
sigma_MIP_10 = std(times_MIP_10,0,2);

mu_MIPREG = mean(times_MIPREG,2);               
sigma_MIPREG = std(times_MIPREG,0,2);         

% Curvas de +1 y -1 desviación
upper_MIP = mu_MIP + sigma_MIP;
lower_MIP = mu_MIP - sigma_MIP;

upper_MIP_10 = mu_MIP_10 + sigma_MIP_10;
lower_MIP_10 = mu_MIP_10 - sigma_MIP_10;

upper_MIPREG = mu_MIPREG + sigma_MIPREG;
lower_MIPREG = mu_MIPREG - sigma_MIPREG;

% Área sombreada (patch o fill)
fill([x fliplr(x)], [upper_MIP' fliplr(lower_MIP')], ...
     rojo_col, 'EdgeColor','none', 'FaceAlpha',0.1);

hold on;
% Media
h1 = plot(x, mu_MIP, 's--', 'LineWidth',2,'Color',rojo_col);
hold on;

fill([x fliplr(x)], [upper_MIP_10' fliplr(lower_MIP_10')], ...
     naranja_col, 'EdgeColor','none', 'FaceAlpha',0.1);
hold on;
% Media
h2 = plot(x, mu_MIP_10, '*-.', 'LineWidth',2,'Color',naranja_col);
hold on;

fill([x fliplr(x)], [upper_MIPREG' fliplr(lower_MIPREG')], ...
     verde_col, 'EdgeColor','none', 'FaceAlpha',0.1);

hold on;
% Media
h3 = plot(x, mu_MIPREG, 'x-', 'LineWidth',2,'Color',verde_col);



grid on;
xl = xlabel(['Budget'],'interpreter','latex');
yl = ylabel('$t_{comp} [s]$','Interpreter','latex');
set(gca, 'FontSize', 9);
set(gca, 'XTick', 1:8, 'XTickLabel', budgets_aval);
set(gca, 'TickLabelInterpreter', 'latex');
legend([h1 h2 h3], {'MIP 30 min','MIP 10 min','MIPREG'}, 'Interpreter','latex','Location','best','FontSize',9 );
xlim([1 8]);


subplot(223);

x = 1:size(links_MIP,1);

% Estadísticos
mu_MIP = mean(links_MIP,2);               
sigma_MIP = std(links_MIP,0,2);   

mu_MIP_10 = mean(links_MIP_10,2);               
sigma_MIP_10 = std(links_MIP_10,0,2); 

mu_MIPREG = mean(links_MIPREG,2);               
sigma_MIPREG = std(links_MIPREG,0,2);         

% Curvas de +1 y -1 desviación
upper_MIP = mu_MIP + sigma_MIP;
lower_MIP = mu_MIP - sigma_MIP;

upper_MIP_10 = mu_MIP_10 + sigma_MIP_10;
lower_MIP_10 = mu_MIP_10 - sigma_MIP_10;

upper_MIPREG = mu_MIPREG + sigma_MIPREG;
lower_MIPREG = mu_MIPREG - sigma_MIPREG;

% Área sombreada (patch o fill)
fill([x fliplr(x)], [upper_MIP' fliplr(lower_MIP')], ...
     rojo_col, 'EdgeColor','none', 'FaceAlpha',0.1);

hold on;
% Media
h1 = plot(x, mu_MIP, 's--', 'LineWidth',2,'Color',rojo_col);

hold on;
% Área sombreada (patch o fill)
fill([x fliplr(x)], [upper_MIP_10' fliplr(lower_MIP_10')], ...
     rojo_col, 'EdgeColor','none', 'FaceAlpha',0.1);

hold on;
% Media
h2 = plot(x, mu_MIP_10, '*-.', 'LineWidth',2,'Color',naranja_col);
hold on;

fill([x fliplr(x)], [upper_MIPREG' fliplr(lower_MIPREG')], ...
     verde_col, 'EdgeColor','none', 'FaceAlpha',0.1);

hold on;
% Media
h3 = plot(x, mu_MIPREG, 'x-', 'LineWidth',2,'Color',verde_col);
grid on;
xl = xlabel(['Budget'],'interpreter','latex');
yl = ylabel('$N_{arcs}$','Interpreter','latex');
set(gca, 'FontSize', 9);
set(gca, 'XTick', 1:8, 'XTickLabel', budgets_aval);
set(gca, 'TickLabelInterpreter', 'latex');
legend([h1 h2 h3], {'MIP 30 min','MIP 10 min','MIPREG'}, 'Interpreter','latex','Location','best','FontSize',9 );
xlim([1 8]);


subplot(224);
x = 1:size(gaps,1);

% Estadísticos
mu_MIP = mean(gaps,2);               
sigma_MIP = std(gaps,0,2);  

mu_MIP_10 = mean(gaps_10,2);               
sigma_MIP_10 = std(gaps_10,0,2); 

% Curvas de +1 y -1 desviación
upper_MIP = mu_MIP + sigma_MIP;
lower_MIP = mu_MIP - sigma_MIP;

upper_MIP_10 = mu_MIP_10 + sigma_MIP_10;
lower_MIP_10 = mu_MIP_10 - sigma_MIP_10;

% Área sombreada (patch o fill)
fill([x fliplr(x)], [upper_MIP' fliplr(lower_MIP')], ...
     rojo_col, 'EdgeColor','none', 'FaceAlpha',0.1);

hold on;
% Media
h1 = plot(x, mu_MIP, 's--', 'LineWidth',2,'Color',rojo_col);

hold on;

fill([x fliplr(x)], [upper_MIP_10' fliplr(lower_MIP_10')], ...
     naranja_col, 'EdgeColor','none', 'FaceAlpha',0.1);

hold on;
% Media
h2 = plot(x, mu_MIP_10, '*-.', 'LineWidth',2,'Color',naranja_col);

grid on;
xl = xlabel(['Budget'],'interpreter','latex');
yl = ylabel('Optimality gap [\%]','Interpreter','latex');
set(gca, 'FontSize', 9);
set(gca, 'XTick', 1:8, 'XTickLabel', budgets_aval);
set(gca, 'TickLabelInterpreter', 'latex');
legend([h1,h2], {'MIP 30 min','MIP 10 min'}, 'Interpreter','latex','Location','best','FontSize',9 )
xlim([1 8]);


saveas(gcf, './6node_random_results.png');

%% Functions


function budget = get_budget(s,s_prim,a,a_prim,n,...
    station_cost,station_capacity_slope,link_cost,link_capacity_slope,lam)
    budget = 0;
    for i=1:n
        if s_prim(i) > 1e-2
            budget = budget + lam*station_cost(i)+ ...
                station_capacity_slope(i)*s_prim(i);
        end
        for j=1:n
            if a_prim(i,j) > 1e-2
                budget = budget + lam*link_cost(i,j)+ ...
                    link_capacity_slope(i,j) * a_prim(i,j);
            end
        end
    end
end

function [pax_obj] = get_entr_val(travel_time,prices,alt_time,alt_price,a_prim,delta_a,...
    s_prim,delta_s,fij,f,fext,demand,dm_pax,dm_op,n)
    
    pax_obj = 0;
    for o=1:n
        for d=1:n
            pax_obj = pax_obj + 1e-6*(demand(o,d).*sum(sum((travel_time+prices).*fij(:,:,o,d)))); 
        end
    end
    pax_obj = pax_obj + 1e-6*(sum(sum(demand.*(alt_time+alt_price).*fext)));
    pax_obj = pax_obj + 1e-6*(sum(sum(demand.*(-entr(f) - f))));
    pax_obj = pax_obj + 1e-6*(sum(sum(demand.*(-entr(fext) - fext))));
    pax_obj = 1e6.*pax_obj;

end

function [obj_val,pax_obj,op_obj] = get_obj_val(alfa, op_link_cost,...
    travel_time,prices,alt_time,alt_price,a_prim,delta_a, ...
    s_prim,delta_s,fij,f,fext,demand)
    n = 6;
    pax_obj = 0;
    op_obj = 0;
    eps = 1e-3;
    dm_pax = 1;
    dm_op = 0.008;


    op_obj = op_obj + (sum(sum(op_link_cost.*a_prim))); %operational costs
    for o=1:n
        for d=1:n
            pax_obj = pax_obj + demand(o,d).*fext(o,d);
        end
    end
    obj_val = (alfa/(dm_pax))*pax_obj + ((1-alfa)/(dm_op))*op_obj;
end



function [sprim,s,deltas,aprim,a,deltaa,f,fext,fij,comp_time] = compute_sim_MIP_entr(lam,beta,alfa,n,budget)

    tic; 
    fid = fopen("./export_txt/lam.txt",'w');
    if fid < 0, error('No puedo abrir %s', filename); end
        fprintf(fid, '%d', lam);
    fclose(fid);

    fid = fopen("./export_txt/alfa.txt",'w');
    if fid < 0, error('No puedo abrir %s', filename); end
        fprintf(fid, '%d', alfa);
    fclose(fid);

    fid = fopen("./export_txt/beta.txt",'w');
    if fid < 0, error('No puedo abrir %s', filename); end
        fprintf(fid, '%d', beta);
    fclose(fid);


    a_prev = 1e4*ones(n);
    s_prev = 1e4*ones(1,n);
    write_gams_param_ii('./export_txt/a_prev.txt', a_prev);
    write_gams_param1d_full('./export_txt/s_prev.txt', s_prev);

    dm_pax = 0.01;
    dm_op = 0.008;

    fid = fopen("./export_txt/dm_pax.txt",'w');
    if fid < 0, error('No puedo abrir %s', filename); end
        fprintf(fid, '%d', dm_pax);
    fclose(fid);

    fid = fopen("./export_txt/dm_op.txt",'w');
    if fid < 0, error('No puedo abrir %s', filename); end
        fprintf(fid, '%d', dm_op);
    fclose(fid);

    fid = fopen("./export_txt/budget.txt",'w');
    if fid < 0, error('No puedo abrir %s', filename); end
        fprintf(fid, '%d', budget);
    fclose(fid);

    % fid = fopen("./export_txt/niters.txt",'w');
    % if fid < 0, error('No puedo abrir %s', filename); end
    %     fprintf(fid, '%d', niters);
    % fclose(fid);


    cmdpath = "cd C:\GAMS\50";

    system(cmdpath);
    
    gmsFile  = 'C:\Users\freal\MATLAB\Projects\untitled\code\6node\mipreg_mosek.gms';

    gamsExe = 'C:\GAMS\50\gams.exe'; 



    cmd = sprintf('%s %s ', ...
              gamsExe, gmsFile);
    [status,out] = system(cmd);
    disp(out);

    results_file_ctime = readtable('./output_all.xlsx','Sheet','solver_time');
    comp_time = table2array(results_file_ctime);


    results_file_sprim = readtable('./output_all.xlsx','Sheet','sprim_level');
    sprim = table2array(results_file_sprim(1,:));

    results_file_s = readtable('./output_all.xlsx','Sheet','s_level');
    s = table2array(results_file_s(1,:));

    results_file_deltas = readtable('./output_all.xlsx','Sheet','deltas_level');
    deltas = table2array(results_file_deltas(1,:));

    results_file_aprim = readtable('./output_all.xlsx','Sheet','aprim_level');
    aprim = table2array(results_file_aprim(1:n,2:(n+1)));

    results_file_a = readtable('./output_all.xlsx','Sheet','a_level');
    a = table2array(results_file_a(1:n,2:(n+1)));

    results_file_deltaa = readtable('./output_all.xlsx','Sheet','deltaa_level');
    deltaa = table2array(results_file_deltaa(1:n,2:(n+1)));

    results_file_f = readtable('./output_all.xlsx','Sheet','f_level');
    f = table2array(results_file_f(1:n,2:(n+1)));

    results_file_fext = readtable('./output_all.xlsx','Sheet','fext_level');
    fext = table2array(results_file_fext(1:n,2:(n+1)));

    T = readtable('fij_long.csv');      % columnas: i, j, o, d, value (strings/números)
    [iU,~,iIdx] = unique(T.i,'stable');
    [jU,~,jIdx] = unique(T.j,'stable');
    [oU,~,oIdx] = unique(T.o,'stable');
    [dU,~,dIdx] = unique(T.d,'stable');

    
    fij = accumarray([iIdx,jIdx,oIdx,dIdx], T.value, ...
                   [numel(iU), numel(jU), numel(oU), numel(dU)], @sum, 0);

    a_bin = zeros(n);
    a_bin(aprim > 1e-2) = 1;

    s_bin = zeros(n,1);
    s_bin(sprim > 1e-2) = 1;

    write_gams_param_ii('./export_txt/a_bin_mipreg.txt', a_bin);
    write_gams_param_ii('./export_txt/a_prim_mipreg.txt', max(0,aprim));
    write_gams_param1d_full('./export_txt/s_bin_mipreg.txt', s_bin);
    write_gams_param1d_full('./export_txt/s_prim_mipreg.txt', max(0,sprim));

    gmsFile  = 'C:\Users\freal\MATLAB\Projects\untitled\code\6node\flow_assignment.gms';

    gamsExe = 'C:\GAMS\50\gams.exe'; 



    cmd = sprintf('%s %s ', ...
              gamsExe, gmsFile);
    [status,out] = system(cmd);
    disp(out);


    results_file_f = readtable('./output_all.xlsx','Sheet','f_level');
    f = table2array(results_file_f(1:n,2:(n+1)));

    results_file_fext = readtable('./output_all.xlsx','Sheet','fext_level');
    fext = table2array(results_file_fext(1:n,2:(n+1)));

    T = readtable('fij_long.csv');      % columnas: i, j, o, d, value (strings/números)
    [iU,~,iIdx] = unique(T.i,'stable');
    [jU,~,jIdx] = unique(T.j,'stable');
    [oU,~,oIdx] = unique(T.o,'stable');
    [dU,~,dIdx] = unique(T.d,'stable');

    
    fij = accumarray([iIdx,jIdx,oIdx,dIdx], T.value, ...
                   [numel(iU), numel(jU), numel(oU), numel(dU)], @sum, 0);



end


function [sprim,s,deltas,aprim,a,deltaa,f,fext,fij] = compute_sim_MIP(lam,beta,alfa,n,budget)

[n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    load_factor,op_link_cost,congestion_coef_stations,...
    congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,eta,...
    a_max,candidasourcertes] = parameters_6node_network();

    tic;

    dm_pax = 0.01;
    dm_op = 0.008;

    fid = fopen("./export_txt/lam.txt",'w');
    if fid < 0, error('No puedo abrir %s', fid); end
        fprintf(fid, '%d', lam);
    fclose(fid);

    fid = fopen("./export_txt/alfa.txt",'w');
    if fid < 0, error('No puedo abrir %s', fid); end
        fprintf(fid, '%d', alfa);
    fclose(fid);

    fid = fopen("./export_txt/beta.txt",'w');
    if fid < 0, error('No puedo abrir %s', fid); end
        fprintf(fid, '%d', beta);
    fclose(fid);

    fid = fopen("./export_txt/budget.txt",'w');
    if fid < 0, error('No puedo abrir %s', fid); end
        fprintf(fid, '%d', budget);
    fclose(fid);

    disp(budget);

    a_prev = 1e4*ones(n);
    s_prev = 1e4*ones(1,n);
    write_gams_param_ii('./export_txt/a_prev.txt', a_prev);
    write_gams_param1d_full('./export_txt/s_prev.txt', s_prev);


    fid = fopen("./export_txt/dm_pax.txt",'w');
    if fid < 0, error('No puedo abrir %s', fid); end
        fprintf(fid, '%d', dm_pax);
    fclose(fid);

    fid = fopen("./export_txt/dm_op.txt",'w');
    if fid < 0, error('No puedo abrir %s', fid); end
        fprintf(fid, '%d', dm_op);
    fclose(fid);


    

    cmdpath = "cd C:\GAMS\50";

    system(cmdpath);
    
    gmsFile  = 'C:\Users\freal\MATLAB\Projects\untitled\code\6node\logit_res_mosek.gms';

    gamsExe = 'C:\GAMS\50\gams.exe'; 



    cmd = sprintf('%s %s ', ...
              gamsExe, gmsFile);
    [status,out] = system(cmd);
    disp(out);

    results_file_sprim = readtable('./output_all.xlsx','Sheet','sprim_level');
    sprim = table2array(results_file_sprim(1,:));

    results_file_s = readtable('./output_all.xlsx','Sheet','s_level');
    s = table2array(results_file_s(1,:));

    results_file_deltas = readtable('./output_all.xlsx','Sheet','deltas_level');
    deltas = table2array(results_file_deltas(1,:));

    results_file_aprim = readtable('./output_all.xlsx','Sheet','aprim_level');
    aprim = table2array(results_file_aprim(1:n,2:(n+1)));

    results_file_a = readtable('./output_all.xlsx','Sheet','a_level');
    a = table2array(results_file_a(1:n,2:(n+1)));

    results_file_deltaa = readtable('./output_all.xlsx','Sheet','deltaa_level');
    deltaa = table2array(results_file_deltaa(1:n,2:(n+1)));

    results_file_f = readtable('./output_all.xlsx','Sheet','f_level');
    f = table2array(results_file_f(1:n,2:(n+1)));

    results_file_fext = readtable('./output_all.xlsx','Sheet','fext_level');
    fext = table2array(results_file_fext(1:n,2:(n+1)));

    results_file_mipgap = readtable('./output_all.xlsx','Sheet','mip_opt_gap');
    mipgap = table2array(results_file_mipgap);

    T = readtable('fij_long.csv');      % columnas: i, j, o, d, value (strings/números)
    [iU,~,iIdx] = unique(T.i,'stable');
    [jU,~,jIdx] = unique(T.j,'stable');
    [oU,~,oIdx] = unique(T.o,'stable');
    [dU,~,dIdx] = unique(T.d,'stable');

    
    fij = accumarray([iIdx,jIdx,oIdx,dIdx], T.value, ...
                   [numel(iU), numel(jU), numel(oU), numel(dU)], @sum, 0);

    results_file_ctime = readtable('./output_all.xlsx','Sheet','solver_time');
    comp_time = table2array(results_file_ctime);
    [obj_val,pax_obj,op_obj] = get_obj_val(alfa, op_link_cost,...
            travel_time,prices,alt_time,alt_price,aprim,deltaa, ...
            sprim,deltas,fij,f,fext,demand);
       budget = get_budget(s,sprim,a,aprim,n,...
            station_cost,station_capacity_slope,link_cost,link_capacity_slope,lam);
   filename = sprintf('./6node_rebutal_MIP_30min/beta=%d_lam=%d.mat',beta,lam);
    save(filename,'s','sprim','deltas', ...
    'a','aprim','deltaa','f','fext','fij','comp_time','budget', ...
    'pax_obj','op_obj','obj_val','mipgap');

end

function [sprim,s,deltas,aprim,a,deltaa,f,fext,fij] = compute_sim_MIP_rand(lam,beta,alfa,n,budget,rr, ...
    op_link_cost, travel_time,prices,alt_time,alt_price,demand,station_cost,link_cost,...
    station_capacity_slope,link_capacity_slope)


    tic;

    fid = fopen("./export_txt/lam.txt",'w');
    if fid < 0, error('No puedo abrir %s', fid); end
        fprintf(fid, '%d', lam);
    fclose(fid);

    fid = fopen("./export_txt/alfa.txt",'w');
    if fid < 0, error('No puedo abrir %s', fid); end
        fprintf(fid, '%d', alfa);
    fclose(fid);

    fid = fopen("./export_txt/beta.txt",'w');
    if fid < 0, error('No puedo abrir %s', fid); end
        fprintf(fid, '%d', beta);
    fclose(fid);

    fid = fopen("./export_txt/budget.txt",'w');
    if fid < 0, error('No puedo abrir %s', fid); end
        fprintf(fid, '%d', budget);
    fclose(fid);

    disp(budget);

    a_prev = 1e4*ones(n);
    s_prev = 1e4*ones(1,n);
    write_gams_param_ii('./export_txt/a_prev.txt', a_prev);
    write_gams_param1d_full('./export_txt/s_prev.txt', s_prev);

    dm_pax = 1e4;
    dm_op = 1e2;

    fid = fopen("./export_txt/dm_pax.txt",'w');
    if fid < 0, error('No puedo abrir %s', fid); end
        fprintf(fid, '%d', dm_pax);
    fclose(fid);

    fid = fopen("./export_txt/dm_op.txt",'w');
    if fid < 0, error('No puedo abrir %s', fid); end
        fprintf(fid, '%d', dm_op);
    fclose(fid);


    cmdpath = "cd C:\GAMS\50";

    system(cmdpath);
    
    gmsFile  = 'C:\Users\freal\MATLAB\Projects\untitled\code\6node\logit_res_mosek.gms';

    gamsExe = 'C:\GAMS\50\gams.exe'; 



    cmd = sprintf('%s %s ', ...
              gamsExe, gmsFile);
    [status,out] = system(cmd);
    disp(out);

    results_file_sprim = readtable('./output_all.xlsx','Sheet','sprim_level');
    sprim = table2array(results_file_sprim(1,:));

    results_file_s = readtable('./output_all.xlsx','Sheet','s_level');
    s = table2array(results_file_s(1,:));

    results_file_deltas = readtable('./output_all.xlsx','Sheet','deltas_level');
    deltas = table2array(results_file_deltas(1,:));

    results_file_aprim = readtable('./output_all.xlsx','Sheet','aprim_level');
    aprim = table2array(results_file_aprim(1:n,2:(n+1)));

    results_file_a = readtable('./output_all.xlsx','Sheet','a_level');
    a = table2array(results_file_a(1:n,2:(n+1)));

    results_file_deltaa = readtable('./output_all.xlsx','Sheet','deltaa_level');
    deltaa = table2array(results_file_deltaa(1:n,2:(n+1)));

    results_file_f = readtable('./output_all.xlsx','Sheet','f_level');
    f = table2array(results_file_f(1:n,2:(n+1)));

    results_file_fext = readtable('./output_all.xlsx','Sheet','fext_level');
    fext = table2array(results_file_fext(1:n,2:(n+1)));

     results_file_mipgap = readtable('./output_all.xlsx','Sheet','mip_opt_gap');
    mipgap = table2array(results_file_mipgap);

    T = readtable('fij_long.csv');      % columnas: i, j, o, d, value (strings/números)
    [iU,~,iIdx] = unique(T.i,'stable');
    [jU,~,jIdx] = unique(T.j,'stable');
    [oU,~,oIdx] = unique(T.o,'stable');
    [dU,~,dIdx] = unique(T.d,'stable');

    
    fij = accumarray([iIdx,jIdx,oIdx,dIdx], T.value, ...
                   [numel(iU), numel(jU), numel(oU), numel(dU)], @sum, 0);

   results_file_ctime = readtable('./output_all.xlsx','Sheet','solver_time');
    comp_time = table2array(results_file_ctime);
   [obj_val,pax_obj,op_obj] = get_obj_val(alfa, op_link_cost,...
        travel_time,prices,alt_time,alt_price,aprim,deltaa, ...
        sprim,deltas,fij,f,fext,demand);
   budget = get_budget(s,sprim,a,aprim,n,...
        station_cost,station_capacity_slope,link_cost,link_capacity_slope,lam);
   filename = sprintf('./6node_rebutal_MIP_30min/beta=%d_lam=%d_run=%d.mat',beta,lam,rr);
    save(filename,'s','sprim','deltas', ...
    'a','aprim','deltaa','f','fext','fij','comp_time','budget', ...
    'pax_obj','op_obj','obj_val','mipgap');

end



function [n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    load_factor,op_link_cost,congestion_coef_stations,...
    congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,eta,...
    a_max,candidates] = parameters_6node_network()

    n = 6;
    
    %candidates to construct a link for each neighbor
    candidates = {[2,3],[1,3,4],[1,2,4,5],[2,3,5,6],[3,4,6],[4,5]};
    
    %cost of using the alternative network for each o-d pair
    alt_cost = [0,1.6,0.8,2,1.6,2.5; 
                2,0,0.9,1.2,1.5,2.5; 
                1.5,1.4,0,1.3,0.9,2; 
                1.9,2,1.9,0,1.8,2; 
                3,1.5,2,2,0,1.5; 
                2.1,2.7,2.2,1,1.5,0];
    
    %fixed cost for constructing links
    link_cost = (1e6/(25*365.25)).*[0,1.7,2.7,0,0,0; 
                 1.7,0,2.1,3,0,0;
                 2.7,2.1,0,2.6,1.7,0; 
                 0,3,2.6,0,2.8,2.4; 
                 0,0,1.7,2.8,0,1.9; 
                 0,0,0,2.4,1.9,0];
    link_cost (link_cost ==0) = 1e4;

    
    %fixed cost for constructing stations
    station_cost = (1e6/(25*365.25)).*[2, 3, 2.2, 3, 2.5, 1.3];
    
    link_capacity_slope = 0.04.*link_cost; 
    station_capacity_slope = 0.04.*station_cost;
    
    %demand between od pairs
    demand = 1e3.*[0,9,26,19,13,12;
              11,0,14,26,7,18;
              30,19,0,30,24,8;
              21,9,11,0,22,16;
              14,14,8,9,0,20;
              26,1,22,24,13,0];
    
    distance = 10000 * ones(n, n); % Distances between arcs
    
    for i = 1:n
        distance(i, i) = 0;
    end
    
    distance(1, 2) = 0.75;
    distance(1, 3) = 0.7;

    
    distance(2, 3) = 0.6;
    distance(2, 4) = 1.1;
    
    distance(3, 4) = 1.1;
    distance(3, 5) = 0.5;

    
    distance(4,5) = 0.8;
    distance(4,6) = 0.7;

    
    distance(5,6) = 0.5;

    
    for i = 1:n
        for j = i+1:n
            distance(j, i) = distance(i, j); % Distances are symmetric
        end
    end
    
    %Load factor on stations
    load_factor = 0.25 .* ones(1, n);
    
    % Op Link Cost
    op_link_cost = 4.*distance;
    
    % Congestion Coefficients
    congestion_coef_stations = 0.1 .* ones(1, n);
    congestion_coef_links = 0.1 .* ones(n);
    
    % Prices
    prices = 0.1.*(distance).^(0.7);
    %prices = zeros(n);
    
    % Travel Time
    travel_time = 60 .* distance ./ 30; % Time in minutes

    
    % Alt Time
    alt_time = 60 .* alt_cost ./ 30; % Time in minutes

    alt_price = 0.1.*(alt_cost).^(0.7); %price

    
    
    a_nom = 588;             
    
    tau = 0.57;
    eta = 0.25;
    a_max = 1e9;
    eps = 1e-3;
end

function [n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    load_factor,op_link_cost,travel_time,alt_time,alt_price,a_nom,tau,eta,...
    a_max,candidates] = parameters_6node_network_rand()

    n = 6;
    
    %candidates to construct a link for each neighbor
    candidates_cell = {[2,3],[1,3,4],[1,2,4,5],[2,3,5,6],[3,4,6],[4,5]};
   %cost of using the alternative network for each o-d pair
    alt_cost = [0,1.6,0.8,2,1.6,2.5; 
                2,0,0.9,1.2,1.5,2.5; 
                1.5,1.4,0,1.3,0.9,2; 
                1.9,2,1.9,0,1.8,2; 
                3,1.5,2,2,0,1.5; 
                2.1,2.7,2.2,1,1.5,0];

    alt_cost = max(0,alt_cost + randn(n).*(1-eye(n)));
    
    %fixed cost for constructing links
    link_cost = (1e6/(25*365.25)).*[0,1.7,2.7,0,0,0; 
                 1.7,0,2.1,3,0,0;
                 2.7,2.1,0,2.6,1.7,0; 
                 0,3,2.6,0,2.8,2.4; 
                 0,0,1.7,2.8,0,1.9; 
                 0,0,0,2.4,1.9,0];
    link_cost (link_cost ==0) = 1e4;
    link_cost = max(0,link_cost + randn(n).*(1-eye(n)));

    
    %fixed cost for constructing stations
    station_cost = [2, 3, 2.2, 3, 2.5, 1.3];

    station_cost = max(0,station_cost + randn(1,n));

    station_cost = station_cost.*(1e6/(25*365.25));
    
    link_capacity_slope = 0.04.*link_cost; 
    station_capacity_slope = 0.04.*station_cost;
    
    %demand between od pairs
    demand = 1e3.*[0,9,26,19,13,12;
              11,0,14,26,7,18;
              30,19,0,30,24,8;
              21,9,11,0,22,16;
              14,14,8,9,0,20;
              26,1,22,24,13,0];

    demand = max(0,demand + 3e3.*randn(n).*(1-eye(n)));
    
    distance = 10000 * ones(n, n); % Distances between arcs
    
    for i = 1:n
        distance(i, i) = 0;
    end
    
    distance(1, 2) = 0.75;
    distance(1, 3) = 0.7;

    
    distance(2, 3) = 0.6;
    distance(2, 4) = 1.1;
    
    distance(3, 4) = 1.1;
    distance(3, 5) = 0.5;

    
    distance(4,5) = 0.8;
    distance(4,6) = 0.7;

    
    distance(5,6) = 0.5;

    
    for i = 1:n
        for j = i+1:n
            distance(j, i) = distance(i, j); % Distances are symmetric
        end
    end
    
    %Load factor on stations
    load_factor = 0.25 .* ones(1, n);
    
    % Op Link Cost
    op_link_cost = 4.*distance;
    

    
    % Prices
    prices = 0.1.*(distance).^(0.7);
    %prices = zeros(n);
    
    % Travel Time
    travel_time = 60 .* distance ./ 30; % Time in minutes

    
    % Alt Time
    alt_time = 60 .* alt_cost ./ 30; % Time in minutes

    alt_price = 0.1.*(alt_cost).^(0.7); %price

    
    
    a_nom = 588;             
    
    tau = 0.57;
    eta = 0.25;
    a_max = 1e9;
    eps = 1e-3;

    M = 1e4;
    nreg = 10;
    eps = 1e-3;
    vals_regs = linspace(0.005,0.995,nreg-1);
    [lin_coef,bord,b] = get_linearization(n,nreg,alt_time,alt_price,-0.2,-0.2,vals_regs);
    
    candidates = zeros(n);
    for i=1:n
        candidates(i,candidates_cell{i}) = 1;
    end
    write_gams_param_iii('./export_txt/lin_coef.txt', lin_coef);
    write_gams_param_iii('./export_txt/b.txt', b);
    write_gams_param_iii('./export_txt/bord.txt', bord);
    
    % === 2D matrices ===
    write_gams_param_ii('./export_txt/demand.txt', demand);
    write_gams_param_ii('./export_txt/travel_time.txt', travel_time);
    write_gams_param_ii('./export_txt/alt_time.txt', alt_time);
    write_gams_param_ii('./export_txt/alt_price.txt', alt_price);
    write_gams_param_ii('./export_txt/link_cost.txt', link_cost);
    write_gams_param_ii('./export_txt/link_capacity_slope.txt', link_capacity_slope);
    write_gams_param_ii('./export_txt/prices.txt', prices);
    write_gams_param_ii('./export_txt/op_link_cost.txt', op_link_cost);
    write_gams_param_ii('./export_txt/candidates.txt', candidates);
    
    
    
    % === 1D vectores ===
    write_gams_param1d_full('./export_txt/station_cost.txt', station_cost);
    write_gams_param1d_full('./export_txt/station_capacity_slope.txt', station_capacity_slope);
end


function [lin_coef,bord,b] = get_linearization(n,nreg,alt_time,alt_price,omega_t,omega_p,vals_regs)
    dmax = zeros(nreg,n,n);
    dmin = dmax;
    lin_coef = dmax;
    bord = zeros(nreg,n,n);

    
    for o=1:n
        for d=1:n
            u = omega_t*alt_time(o,d) + omega_p*alt_price(o,d);
            for r=1:(nreg-1)
                dmax(r,o,d) = min(0,u + log(vals_regs(r)/(1-vals_regs(r)))  );
            end
            dmax(nreg,o,d) = 0;
            dmin(1,o,d) = -3e1;
            for r=2:nreg
                dmin(r,o,d) = dmax(r-1,o,d);
            end

            for r=2:(nreg-1)
                if (dmax(r,o,d) == dmin(r,o,d))
                    lin_coef(r,o,d) = 0;
                    bord(r,o,d) = vals_regs(r);
                else
                    lin_coef(r,o,d) = (vals_regs(r)-vals_regs(r-1))/(dmax(r,o,d)-dmin(r,o,d));
                    bord(r,o,d) = vals_regs(r-1);
                end
            end
            lin_coef(1,o,d) = (vals_regs(1))/(dmax(1,o,d)-dmin(1,o,d));
            bord(1,o,d) = 0;
            if dmin(nreg,o,d)==0
                lin_coef(nreg,o,d) = 0;
            else
                lin_coef(nreg,o,d) = (1-vals_regs(nreg-1))/(0-dmin(nreg,o,d));
            end
            bord(nreg,o,d) = vals_regs(nreg-1);
        end
    end
    b = dmin;


end

function val = logit(x,omega_t,omega_p,time,price)
    val = exp(x)./( exp(x) + exp(omega_t*time + omega_p*price) );
end

function write_gams_param_iii(filename, M)
    % filename: ruta del .txt (p.ej. 'demand.txt')
    % M: matriz NxN (puede ser sparse)
    % zero_tol: umbral para considerar cero (p.ej. 0 o 1e-12)

    [n1, n2, n3] = size(M);


    fid = fopen(filename,'w');
    if fid < 0, error('No puedo abrir %s', filename); end

    % Si M es dispersa, recorre solo no-ceros

     for s=1:n1
            for r = 1:n2
                for c = 1:n3 
                    val = M(s,r,c);
                        fprintf(fid, 'seg%d.i%d.i%d %.12g\n', s, r, c, val);
                end
            end
     end
    fclose(fid);
end


function write_gams_param_ii(filename, M)
    % filename: ruta del .txt (p.ej. 'demand.txt')
    % M: matriz NxN (puede ser sparse)
    % zero_tol: umbral para considerar cero (p.ej. 0 o 1e-12)

    if nargin < 3; end
    [n1, n2] = size(M);
    if n1 ~= n2
        error('M debe ser cuadrada para dominio (i,i).');
    end

    fid = fopen(filename,'w');
    if fid < 0, error('No puedo abrir %s', filename); end

    % Si M es dispersa, recorre solo no-ceros
    if issparse(M)
        [r,c,v] = find(M);
        for k = 1:numel(v)
                fprintf(fid, 'i%d.i%d %.12g\n', r(k), c(k), v(k));
        end
    else
        for r = 1:n1
            for c = 1:n2
                val = M(r,c);
                    fprintf(fid, 'i%d.i%d %.12g\n', r, c, val);
            end
        end
    end
    fclose(fid);
end

function write_gams_param1d_full(filename, v)
% Escribe un parámetro 1D en formato GAMS:
%   Parameter <paramName> /
%   i1 <valor>
%   i2 <valor>
%   ...
%   /;
%
% filename : ruta del .txt (p.ej. 'station_cost.txt')
% v        : vector (Nx1 o 1xN)
% paramName: nombre del parámetro en GAMS (def: 'station_cost')
% prefix   : prefijo de la etiqueta (def: 'i' -> i1, i2, ...)
% zero_tol : umbral para omitir ~0 (def: 0 => no escribe los ceros exactos)
v = v(:);
n = numel(v);

fid = fopen(filename,'w');
if fid < 0, error('No puedo abrir %s', filename); end

for k = 1:n
    val = v(k);
    fprintf(fid, 'i%d %.12g\n', k, val);
end
fclose(fid);
end

