close all;
clear all;
clc;

% === 1) RUTA de tu GAMS ===
gamsHome = 'C:\GAMS\50';  % <-- CAMBIA esto

[basedir,~,~] = fileparts(mfilename('fullpath'));
basedir = fullfile(basedir, 'export_csv');   % subcarpeta de exportación
if ~exist(basedir,'dir'); mkdir(basedir); end




%%

% === Tus datos ===
[n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    load_factor,op_link_cost,congestion_coef_stations,...
    congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,eta,...
    a_max,candidasourcertes] = parameters_9node_network();

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



% === 1D vectores ===
write_gams_param1d_full('./export_txt/station_cost.txt', station_cost);
write_gams_param1d_full('./export_txt/station_capacity_slope.txt', station_capacity_slope);

alfa = 0.5;

%%

betas = [0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.5]; %estos son los de verdad
lams = [5,6];
runs = 1:10;

niters = 10;
budgets = zeros(length(betas),length(lams),length(runs));
for ll=1:length(lams)
    lam = lams(ll);
    for bb=1:length(betas)
        for rr=1:length(runs)
            eps = 1e-3;
            alfa = 0.5;
            beta = betas(bb);
            rng(rr);
           [n,link_cost,station_cost,link_capacity_slope,...
                station_capacity_slope,demand,prices,...
                load_factor,op_link_cost,congestion_coef_stations,...
                congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,eta,...
                a_max,candidates] = parameters_9node_network_rand();

           [sprim,s,deltas,aprim,a,deltaa,f,fext,fij,comp_time] = compute_sim(niters,lam,beta,alfa,n);
          [obj_val,pax_obj,op_obj] = get_obj_val(alfa, op_link_cost,...
        travel_time,prices,alt_time,alt_price,aprim,deltaa, ...
        sprim,deltas,fij,f,fext,demand,1,1);
           budget = get_budget(s,sprim,a,aprim,n,...
                station_cost,station_capacity_slope,link_cost,link_capacity_slope,lam);
           budgets(bb,ll,rr) = budget;
           filename = sprintf('./9node_rebutal_ICVX/beta=%d_lam=%d_run=%d.mat',beta,lam,rr);
            save(filename,'s','sprim','deltas', ...
            'a','aprim','deltaa','f','fext','fij','comp_time','budget', ...
            'pax_obj','op_obj','obj_val');
            
            budget = budgets(bb,ll,rr);
    
            [sprim,s,deltas,aprim,a,deltaa,f,fext,fij,comp_time] = compute_sim_MIP_entr(lam,beta,alfa,n,budget);
            [obj_val,pax_obj,op_obj] = get_obj_val(alfa, op_link_cost,...
                travel_time,prices,alt_time,alt_price,aprim,deltaa, ...
                sprim,deltas,fij,f,fext,demand,1,1);
               budget = get_budget(s,sprim,a,aprim,n,...
                    station_cost,station_capacity_slope,link_cost,link_capacity_slope,lam);
               filename = sprintf('./9node_rebutal_MIPREG/beta=%d_lam=%d_run=%d.mat',beta,lam,rr);
                save(filename,'s','sprim','deltas', ...
                'a','aprim','deltaa','f','fext','fij','comp_time','budget', ...
                'pax_obj','op_obj','obj_val');
            
            disp(['bb = ',num2str(bb),', budget = ',num2str(budget),', obj_val = ',num2str(obj_val),', pax_obj = ',num2str(pax_obj), ...
                ', op_obj = ',num2str(op_obj),', nlinks = ', num2str(sum(sum(a > eps)))]);
        end
        
    end
end

%%
betas = [0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.5]; %estos son los de verdad
runs = 1:10;
lams = [5,6];
budgets = zeros(length(betas),length(lams),length(runs));
dif = zeros(length(betas),length(lams),length(runs));
opt_gaps = dif;
times_MIPREG = dif;
times_ICVX = dif;
links_MIPREG = dif;
links_ICVX = dif;

SD_MIPREG = dif;
SD_ICVX = dif;



for ll=1:length(lams)
    lam = lams(ll);
    for bb=1:length(betas)
        for rr=1:length(runs)
            rng(rr);
           [n,link_cost,station_cost,link_capacity_slope,...
                station_capacity_slope,demand,prices,...
                load_factor,op_link_cost,congestion_coef_stations,...
                congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,eta,...
                a_max,candidates] = parameters_9node_network_rand();
            beta = betas(bb);
            filename = sprintf('./9node_rebutal_MIPREG/beta=%d_lam=%d_run=%d.mat',beta,lam,rr);
            load(filename);    
            budgets(bb,ll,rr) = budget;
            links_MIPREG(bb,ll,rr) = sum(sum(a>1e-2));
            obj_val_MIPREG = obj_val;
            times_MIPREG(bb,ll,rr) = comp_time;

            SD_MIPREG(bb,ll,rr) = 100.*sum(sum(f.*demand))./(sum(sum(demand)));
            
            filename = sprintf('./9node_rebutal_ICVX/beta=%d_lam=%d_run=%d.mat',beta,lam,rr);
            load(filename);
            times_ICVX(bb,ll,rr) = comp_time;
            links_ICVX(bb,ll,rr) = sum(sum(a>1e-2));
            obj_val_ICVX = obj_val;
            dif(bb,ll,rr) = 100*(obj_val_ICVX-obj_val_MIPREG)/obj_val_MIPREG;
            SD_ICVX(bb,ll,rr) = 100.*sum(sum(f.*demand))./(sum(sum(demand)));
        end
    end
end
%%
for ll=1:length(lams)
    lam = lams(ll);
    for bb=1:length(betas)
        beta = betas(bb);
        budget = mean(budgets(bb,ll,:));
        budget_std = std(budgets(bb,ll,:));
        dif_val = mean(dif(bb,ll,:));
        dif_std = std(dif(bb,ll,:));
        nlinks_MIPREG = mean(links_MIPREG(bb,ll,:));
        nlinks_ICVX = mean(links_ICVX(bb,ll,:));
        nlinks_MIPREG_std = std(links_MIPREG(bb,ll,:));
        nlinks_ICVX_std = std(links_ICVX(bb,ll,:));
        times_MIPREG_val = mean(times_MIPREG(bb,ll,:));
        times_MIPREG_std = std(times_MIPREG(bb,ll,:));
        times_ICVX_val = mean(times_ICVX(bb,ll,:));
        times_ICVX_std = std(times_ICVX(bb,ll,:));

        disp([ sprintf('%d',lam),'&', sprintf('%.2f',beta),'&',...
            sprintf('%.2e $\\pm$ %.0e', budget,budget_std),'&',sprintf('%.2f $\\pm$ %.1f',dif_val,dif_std),'&',...
            sprintf('%.2f $\\pm$ %.2f', times_ICVX_val,times_ICVX_std),'&',sprintf('%.2e $\\pm$ %.1e', times_MIPREG_val,times_MIPREG_std),'&',...
            sprintf('%.2f $\\pm$ %.1f', nlinks_ICVX,nlinks_ICVX_std),'&',sprintf('%.2f $\\pm$ %.1f', nlinks_MIPREG,nlinks_MIPREG_std),'&',...
            sprintf('%.2f $\\pm$ %.1f', mean(SD_ICVX(bb,ll,:)),std(SD_ICVX(bb,ll,:))),'&',sprintf('%.2f $\\pm$ %.1f', mean(SD_MIPREG(bb,ll,:)),std(SD_MIPREG(bb,ll,:))),...
            '\\ \hline']);
    
    end
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
[nObs, nX] = size(squeeze(dif(:,1,:))');

posA = (1:nX) - 0.15;   % desplazamiento izquierda para A
posB = (1:nX) + 0.15;   % desplazamiento derecha para B

cla; hold on
boxplot(squeeze(dif(:,1,:))', 'positions', posA, 'Widths', 0.25, 'Colors',[0 0 0], 'Symbol',''); 
boxplot(squeeze(dif(:,2,:))', 'positions', posB, 'Widths', 0.25, 'Colors',[0.5 0.5 0.5], 'Symbol','');

xlim([0.5, nX+0.5]); xticks(1:nX)
eur =['[',char(8364),']'];
xl = xlabel(['$\beta$'],'interpreter','latex');
yl = ylabel('Diff [\%]','Interpreter','latex');
set(gca, 'FontSize', 9);
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'XTick', 1:9, 'XTickLabel', betas);

% leyenda “dummy”
plot(nan,nan,'Color',[0 0 0],'LineWidth',2); plot(nan,nan,'Color',[0.5 0.5 0.5],'LineWidth',2);
legend({'ICVX w.r.t. MIPREG, $\lambda = 5$','ICVX w.r.t. MIPREG, $\lambda = 6$'},'Location','best','Interpreter','latex','FontSize',9); hold off
xlim([0.7 9.3]);

subplot(222);

% times_MIPREG, times_ICVX: [nBeta x nLambda(=2) x nRuns]
[nB, nL, ~] = size(times_MIPREG);
assert(nL==2, 'Se esperaban exactamente 2 valores de lambda');

x = 1:nB;

% Medias sobre runs -> matrices [nB x 2]
mu_MIP  = mean(times_MIPREG, 3, 'omitnan');    % [nB x 2]
sd_MIP  = std(times_MIPREG, 0, 3, 'omitnan');  % [nB x 2]

mu_ICVX = mean(times_ICVX, 3, 'omitnan');
sd_ICVX = std(times_ICVX, 0, 3, 'omitnan');

% --- Límites ±1σ ---
upper_MIP = mu_MIP + sd_MIP;
lower_MIP = mu_MIP - sd_MIP;

upper_ICVX = mu_ICVX + sd_ICVX;
lower_ICVX = mu_ICVX - sd_ICVX;

% Etiquetas de lambda (opcional, si tienes vector 'lambdas' o 'lambda_vals')
if exist('lambdas','var') && numel(lambdas)>=2
    lbl1 = sprintf('\\lambda=%g', lambdas(1));
    lbl2 = sprintf('\\lambda=%g', lambdas(2));
elseif exist('lambda_vals','var') && numel(lambda_vals)>=2
    lbl1 = sprintf('\\lambda=%g', lambda_vals(1));
    lbl2 = sprintf('\\lambda=%g', lambda_vals(2));
else
    lbl1 = '$\lambda = 5$'; lbl2 = '$\lambda=6$';
end

hold on;
% MIPREG en verde: dos lambdas
fill([x fliplr(x)], [upper_MIP(:,1)' fliplr(lower_MIP(:,1)')], ...
    verde_col, 'EdgeColor','none', 'FaceAlpha',0.2); hold on;
h1 = plot(x, mu_MIP(:,1), 'x-', 'LineWidth',2, 'Color', verde_col, 'DisplayName', ['MIPREG, ' lbl1]);
hold on;



fill([x fliplr(x)], [upper_MIP(:,2)' fliplr(lower_MIP(:,2)')], ...
    verde_col, 'EdgeColor','none', 'FaceAlpha',0.2); hold on;
h2 = plot(x, mu_MIP(:,2), 'x:', 'LineWidth',2, 'Color', verde_col, 'DisplayName', ['MIPREG, ' lbl2]);

% ICVX en azul: dos lambdas
hold on;
fill([x fliplr(x)], [upper_ICVX(:,1)' fliplr(lower_ICVX(:,1)')], ...
    azul_col, 'EdgeColor','none', 'FaceAlpha',0.2); hold on;
h3 = plot(x, mu_ICVX(:,1), 'o-',  'LineWidth',2, 'Color', azul_col,  'DisplayName', ['ICVX, ' lbl1]);

hold on;
fill([x fliplr(x)], [upper_ICVX(:,2)' fliplr(lower_ICVX(:,2)')], ...
    azul_col, 'EdgeColor','none', 'FaceAlpha',0.2); hold on;
h4 = plot(x, mu_ICVX(:,2), 'Marker','o','LineStyle',':',  'LineWidth',2, 'Color', azul_col,  'DisplayName', ['ICVX, ' lbl2]);

grid on;
xlabel(['$\beta$'],'Interpreter','latex');
ylabel('$t_{comp} \,[s]$','Interpreter','latex');
set(gca, 'FontSize', 9, 'TickLabelInterpreter', 'latex');

% Eje X con tus betas si existen
if exist('betas','var')
    set(gca, 'XTick', 1:numel(betas), 'XTickLabel', betas);
else
    set(gca, 'XTick', 1:nB);
end

legend([h1 h2 h3 h4], 'Interpreter','latex', 'Location','best', 'FontSize',9);
xlim([1 nB]);
ylim([-2 300]);
hold off;


subplot(223);

% times_MIPREG, times_ICVX: [nBeta x nLambda(=2) x nRuns]
[nB, nL, ~] = size(links_MIPREG);
assert(nL==2, 'Se esperaban exactamente 2 valores de lambda');

x = 1:nB;

% Medias sobre runs -> matrices [nB x 2]
mu_MIP  = mean(links_MIPREG, 3, 'omitnan');    % [nB x 2]
sd_MIP  = std(links_MIPREG, 0, 3, 'omitnan');  % [nB x 2]

mu_ICVX = mean(links_ICVX, 3, 'omitnan');
sd_ICVX = std(links_ICVX, 0, 3, 'omitnan');

% --- Límites ±1σ ---
upper_MIP = mu_MIP + sd_MIP;
lower_MIP = mu_MIP - sd_MIP;

upper_ICVX = mu_ICVX + sd_ICVX;
lower_ICVX = mu_ICVX - sd_ICVX;

% Etiquetas de lambda (opcional, si tienes vector 'lambdas' o 'lambda_vals')
if exist('lambdas','var') && numel(lambdas)>=2
    lbl1 = sprintf('\\lambda=%g', lambdas(1));
    lbl2 = sprintf('\\lambda=%g', lambdas(2));
elseif exist('lambda_vals','var') && numel(lambda_vals)>=2
    lbl1 = sprintf('\\lambda=%g', lambda_vals(1));
    lbl2 = sprintf('\\lambda=%g', lambda_vals(2));
else
    lbl1 = '$\lambda = 5$'; lbl2 = '$\lambda=6$';
end

hold on;
% MIPREG en verde: dos lambdas
fill([x fliplr(x)], [upper_MIP(:,1)' fliplr(lower_MIP(:,1)')], ...
    verde_col, 'EdgeColor','none', 'FaceAlpha',0.2); hold on;
h1 = plot(x, mu_MIP(:,1), 'x-', 'LineWidth',2, 'Color', verde_col, 'DisplayName', ['MIPREG, ' lbl1]);
hold on;



fill([x fliplr(x)], [upper_MIP(:,2)' fliplr(lower_MIP(:,2)')], ...
    verde_col, 'EdgeColor','none', 'FaceAlpha',0.2); hold on;
h2 = plot(x, mu_MIP(:,2), 'x:', 'LineWidth',2, 'Color', verde_col, 'DisplayName', ['MIPREG, ' lbl2]);

% ICVX en azul: dos lambdas
hold on;
fill([x fliplr(x)], [upper_ICVX(:,1)' fliplr(lower_ICVX(:,1)')], ...
    azul_col, 'EdgeColor','none', 'FaceAlpha',0.2); hold on;
h3 = plot(x, mu_ICVX(:,1), 'o-',  'LineWidth',2, 'Color', azul_col,  'DisplayName', ['ICVX, ' lbl1]);

hold on;
fill([x fliplr(x)], [upper_ICVX(:,2)' fliplr(lower_ICVX(:,2)')], ...
    azul_col, 'EdgeColor','none', 'FaceAlpha',0.2); hold on;
h4 = plot(x, mu_ICVX(:,2), 'Marker','o','LineStyle',':',  'LineWidth',2, 'Color', azul_col,  'DisplayName', ['ICVX, ' lbl2]);

grid on;
xlabel(['$\beta$'],'Interpreter','latex');
ylabel('$N_{arcs}$','Interpreter','latex');
set(gca, 'FontSize', 9, 'TickLabelInterpreter', 'latex');

% Eje X con tus betas si existen
if exist('betas','var')
    set(gca, 'XTick', 1:numel(betas), 'XTickLabel', betas);
else
    set(gca, 'XTick', 1:nB);
end

legend([h1 h2 h3 h4], 'Interpreter','latex', 'Location','best', 'FontSize',9);
xlim([1 nB]);
ylim([0 30]);
hold off;



subplot(224);
% times_MIPREG, times_ICVX: [nBeta x nLambda(=2) x nRuns]
[nB, nL, ~] = size(SD_MIPREG);
assert(nL==2, 'Se esperaban exactamente 2 valores de lambda');

x = 1:nB;

% Medias sobre runs -> matrices [nB x 2]
mu_MIP  = mean(SD_MIPREG, 3, 'omitnan');    % [nB x 2]
sd_MIP  = std(SD_MIPREG, 0, 3, 'omitnan');  % [nB x 2]

mu_ICVX = mean(SD_ICVX, 3, 'omitnan');
sd_ICVX = std(SD_ICVX, 0, 3, 'omitnan');

% --- Límites ±1σ ---
upper_MIP = mu_MIP + sd_MIP;
lower_MIP = mu_MIP - sd_MIP;

upper_ICVX = mu_ICVX + sd_ICVX;
lower_ICVX = mu_ICVX - sd_ICVX;

% Etiquetas de lambda (opcional, si tienes vector 'lambdas' o 'lambda_vals')
if exist('lambdas','var') && numel(lambdas)>=2
    lbl1 = sprintf('\\lambda=%g', lambdas(1));
    lbl2 = sprintf('\\lambda=%g', lambdas(2));
elseif exist('lambda_vals','var') && numel(lambda_vals)>=2
    lbl1 = sprintf('\\lambda=%g', lambda_vals(1));
    lbl2 = sprintf('\\lambda=%g', lambda_vals(2));
else
    lbl1 = '$\lambda = 5$'; lbl2 = '$\lambda=6$';
end

hold on;
% MIPREG en verde: dos lambdas
fill([x fliplr(x)], [upper_MIP(:,1)' fliplr(lower_MIP(:,1)')], ...
    verde_col, 'EdgeColor','none', 'FaceAlpha',0.2); hold on;
h1 = plot(x, mu_MIP(:,1), 'x-', 'LineWidth',2, 'Color', verde_col, 'DisplayName', ['MIPREG, ' lbl1]);
hold on;



fill([x fliplr(x)], [upper_MIP(:,2)' fliplr(lower_MIP(:,2)')], ...
    verde_col, 'EdgeColor','none', 'FaceAlpha',0.2); hold on;
h2 = plot(x, mu_MIP(:,2), 'x:', 'LineWidth',2, 'Color', verde_col, 'DisplayName', ['MIPREG, ' lbl2]);

% ICVX en azul: dos lambdas
hold on;
fill([x fliplr(x)], [upper_ICVX(:,1)' fliplr(lower_ICVX(:,1)')], ...
    azul_col, 'EdgeColor','none', 'FaceAlpha',0.2); hold on;
h3 = plot(x, mu_ICVX(:,1), 'o-',  'LineWidth',2, 'Color', azul_col,  'DisplayName', ['ICVX, ' lbl1]);

hold on;
fill([x fliplr(x)], [upper_ICVX(:,2)' fliplr(lower_ICVX(:,2)')], ...
    azul_col, 'EdgeColor','none', 'FaceAlpha',0.2); hold on;
h4 = plot(x, mu_ICVX(:,2), 'Marker','o','LineStyle',':',  'LineWidth',2, 'Color', azul_col,  'DisplayName', ['ICVX, ' lbl2]);

grid on;
xlabel(['$\beta$'],'Interpreter','latex');
ylabel('SD [\%]','Interpreter','latex');
set(gca, 'FontSize', 9, 'TickLabelInterpreter', 'latex');

% Eje X con tus betas si existen
if exist('betas','var')
    set(gca, 'XTick', 1:numel(betas), 'XTickLabel', betas);
else
    set(gca, 'XTick', 1:nB);
end

legend([h1 h2 h3 h4], 'Interpreter','latex', 'Location','best', 'FontSize',9);
xlim([1 nB]);
ylim([0 42]);
hold off;

saveas(gcf, './9node_random_results.png');


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
    s_prim,delta_s,fij,f,fext,demand,dm_pax,dm_op)
    n = 9;
    pax_obj = 0;
    op_obj = 0;
    eps = 1e-3;
    dm_pax = 1.2;
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


    %a_prev = 1e4*ones(n);
    %s_prev = 1e4*ones(1,n);
    %write_gams_param_ii('./export_txt/a_prev.txt', a_prev);
    %write_gams_param1d_full('./export_txt/s_prev.txt', s_prev);

    dm_pax = 1.2;
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
    
    gmsFile  = 'C:\Users\freal\MATLAB\Projects\untitled\code\9node\mipreg_mosek.gms';

    gamsExe = 'C:\GAMS\50\gams.exe'; 



    cmd = sprintf('%s %s ', ...
              gamsExe, gmsFile);
    [status,out] = system(cmd);
    disp(out);


    results_file_ctime = readtable('./output_all.xlsx','Sheet','solver_time');
    comp_time = table2array(results_file_ctime);
    %comp_time = toc;


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


end


function [sprim,s,deltas,aprim,a,deltaa,f,fext,fij] = compute_sim_MIP(lam,beta,alfa,n,budget)

    [n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    load_factor,op_link_cost,travel_time,alt_time,alt_price,a_nom,tau,eta,...
    a_max,candidasourcertes] = parameters_3node_network();

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
    
    gmsFile  = 'C:\Users\freal\MATLAB\Projects\untitled\code\3node\logit_res_mosek.gms';

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

   comp_time = toc;
    [obj_val,pax_obj,op_obj] = get_obj_val(alfa, op_link_cost,...
            travel_time,prices,alt_time,alt_price,aprim,deltaa, ...
            sprim,deltas,fij,f,fext,demand);
       budget = get_budget(s,sprim,a,aprim,n,...
            station_cost,station_capacity_slope,link_cost,link_capacity_slope,lam);
   filename = sprintf('./3node_rebutal_MIP_10min/beta=%d_lam=%d.mat',beta,lam);
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
    
    gmsFile  = 'C:\Users\freal\MATLAB\Projects\untitled\code\3node\logit_res_mosek.gms';

    gamsExe = 'C:\GAMS\50\gams.exe'; 



    cmd = sprintf('%s %s ', ...
              gamsExe, gmsFile);
    [status,out] = system(cmd);
    disp(out);
    %comp_time = toc;

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

    results_file_ctime = readtable('./output_all.xlsx','Sheet','solver_time');
    comp_time = table2array(results_file_ctime);

   [obj_val,pax_obj,op_obj] = get_obj_val(alfa, op_link_cost,...
        travel_time,prices,alt_time,alt_price,aprim,deltaa, ...
        sprim,deltas,fij,f,fext,demand,1,1);
   budget = get_budget(s,sprim,a,aprim,n,...
        station_cost,station_capacity_slope,link_cost,link_capacity_slope,lam);
   filename = sprintf('./3node_rebutal_MIP_10min/beta=%d_lam=%d_run=%d.mat',beta,lam,rr);
    save(filename,'s','sprim','deltas', ...
    'a','aprim','deltaa','f','fext','fij','comp_time','budget', ...
    'pax_obj','op_obj','obj_val');

end


function [n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    load_factor,op_link_cost,travel_time,alt_time,alt_price,a_nom,tau,eta,...
    a_max,candidates] = parameters_3node_network()

    n = 3;
    
    %candidates to construct a link for each neighbor
    candidates = {[2,3],[1,3],[1,2]};
    
    %cost of using the alternative network for each o-d pair
    alt_cost = [0,1.5,2;...
                1.5,0,3;...
                2,3,0];
    
    %fixed cost for constructing links
    link_cost = [0,2,5;...
                2,0,1;...
                5,1,0];
    link_cost (link_cost ==0) = 1e4;

    
    %fixed cost for constructing stations
    station_cost = [2,3,2];
    
    link_capacity_slope = link_cost; 
    station_capacity_slope = station_cost;
    
    %demand between od pairs
    demand = 1e3.*[0,2,1;...
                   2,0,4;...
                   1,4,0];
    
    distance = [0,1,1;...
                1,0,2;...
                1,2,0];
    
    
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
end

function [n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    load_factor,op_link_cost,travel_time,alt_time,alt_price,a_nom,tau,eta,...
    a_max,candidates] = parameters_3node_network_rand()

    n = 3;
    
    %candidates to construct a link for each neighbor
    candidates_cell = {[2,3],[1,3],[1,2]};
    
    %cost of using the alternative network for each o-d pair
    alt_cost = [0,1.5,2;...
                1.5,0,3;...
                2,3,0];

    alt_cost = max(0.2,alt_cost + randn(n).*(1-eye(n)));
    
    %fixed cost for constructing links
    link_cost = [0,2,5;...
                2,0,1;...
                5,1,0];
    link_cost (link_cost ==0) = 1e4;

    link_cost = max(0,link_cost + randn(n).*(1-eye(n)));

    
    %fixed cost for constructing stations
    station_cost = [2,3,2];

    station_cost = max(0,station_cost + randn(1,n));
    
    link_capacity_slope = link_cost; 
    station_capacity_slope = station_cost;
    
    %demand between od pairs
    demand = 1e3.*[0,2,1;...
                   2,0,4;...
                   1,4,0];

    demand = max(0,demand + 1e3.*randn(n).*(1-eye(n)));
    
    distance = [0,1,1;...
                1,0,2;...
                1,2,0];
    
    
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
    nreg = 200;
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



function [n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    load_factor,op_link_cost,congestion_coef_stations,...
    congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,eta,...
    a_max,candidates] = parameters_9node_network()

    n = 9;
    
    %candidates to construct a link for each neighbor
    candidates = {[2,3,9],[1,3,4],[9,1,2,4,5],[2,3,5,6,8],[3,4,6,7],[4,5,7,8],[5,6],[4,6],[1,3]};
    
    %cost of using the alternative network for each o-d pair
    alt_cost = [0,1.6,0.8,2,1.6,2.5,3,2.5,0.8; 
                2,0,0.9,1.2,1.5,2.5,2.7,2.4,1.8; 
                1.5,1.4,0,1.3,0.9,2,1.6,2.3,0.9; 
                1.9,2,1.9,0,1.8,2,1.9,1.2,2; 
                3,1.5,2,2,0,1.5,1.1,1.8,1.7; 
                2.1,2.7,2.2,1,1.5,0,0.9,0.9,2.9; 
                2.8,2.3,1.5,1.8,0.9,0.8,0,1.3,2.1; 
                2.8,2.2,2,1.1,1.5,0.8,1.9,0,0.3; 
                1,1.5,1.1,2.7,1.9,1.8,2.4,3,0];
    
    %fixed cost for constructing links
    link_cost = (1e6/(25*365.25)).*[0,1.7,2.7,0,0,0,0,0,2.9; 
                 1.7,0,2.1,3,0,0,0,0,0; 
                 2.7,2.1,0,2.6,1.7,0,0,0,2.5; 
                 0,3,2.6,0,2.8,2.4,0,3.2,0; 
                 0,0,1.7,2.8,0,1.9,3,0,0; 
                 0,0,0,2.4,1.9,0,2.7,2.8,0; 
                 0,0,0,0,3,2.7,0,0,0; 
                 0,0,0,3.2,0,2.8,0,0,0; 
                 2.9,0,2.5,0,0,0,0,0,0];
    link_cost (link_cost ==0) = 1e4;

    
    %fixed cost for constructing stations
    station_cost = (1e6/(25*365.25)).*[2, 3, 2.2, 3, 2.5, 1.3, 2.8, 2.2, 3.1];
    
    link_capacity_slope = 0.04.*link_cost; 
    station_capacity_slope = 0.04.*station_cost;
    
    %demand between od pairs
    demand = 1e3.*[0,9,26,19,13,12,13,8,11;
              11,0,14,26,7,18,3,6,12;
              30,19,0,30,24,8,15,12,5;
              21,9,11,0,22,16,25,21,23;
              14,14,8,9,0,20,16,22,21;
              26,1,22,24,13,0,16,14,12;
              8,6,9,23,6,13,0,11,11;
              9,2,14,20,18,16,11,0,4;
              8,7,11,22,27,17,8,12,0];
    
    distance = 10000 * ones(n, n); % Distances between arcs
    
    for i = 1:n
        distance(i, i) = 0;
    end
    
    distance(1, 2) = 0.75;
    distance(1, 3) = 0.7;
    distance(1, 9) = 0.9;
    
    distance(2, 3) = 0.6;
    distance(2, 4) = 1.1;
    
    distance(3, 4) = 1.1;
    distance(3, 5) = 0.5;
    distance(3, 9) = 0.7;
    
    distance(4,5) = 0.8;
    distance(4,6) = 0.7;
    distance(4,8) = 0.8;
    
    distance(5,6) = 0.5;
    distance(5,7) = 0.7;
    
    distance(6,7) = 0.5;
    distance(6,8) = 0.4;
    
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
    load_factor,op_link_cost,congestion_coef_stations,...
    congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,eta,...
    a_max,candidates] = parameters_9node_network_rand()

    n = 9;
    
    %candidates to construct a link for each neighbor
    candidasourcertes = {[2,3,9],[1,3,4],[9,1,2,4,5],[2,3,5,6,8],[3,4,6,7],[4,5,7,8],[5,6],[4,6],[1,3]};
    
    %cost of using the alternative network for each o-d pair
    alt_cost = [0,1.6,0.8,2,1.6,2.5,3,2.5,0.8; 
                2,0,0.9,1.2,1.5,2.5,2.7,2.4,1.8; 
                1.5,1.4,0,1.3,0.9,2,1.6,2.3,0.9; 
                1.9,2,1.9,0,1.8,2,1.9,1.2,2; 
                3,1.5,2,2,0,1.5,1.1,1.8,1.7; 
                2.1,2.7,2.2,1,1.5,0,0.9,0.9,2.9; 
                2.8,2.3,1.5,1.8,0.9,0.8,0,1.3,2.1; 
                2.8,2.2,2,1.1,1.5,0.8,1.9,0,0.3; 
                1,1.5,1.1,2.7,1.9,1.8,2.4,3,0];

    alt_cost = max(0.2,alt_cost + 0.5.*randn(n));
    
    %fixed cost for constructing links
    link_cost = (1e6/(25*365.25)).*[0,1.7,2.7,0,0,0,0,0,2.9; 
                 1.7,0,2.1,3,0,0,0,0,0; 
                 2.7,2.1,0,2.6,1.7,0,0,0,2.5; 
                 0,3,2.6,0,2.8,2.4,0,3.2,0; 
                 0,0,1.7,2.8,0,1.9,3,0,0; 
                 0,0,0,2.4,1.9,0,2.7,2.8,0; 
                 0,0,0,0,3,2.7,0,0,0; 
                 0,0,0,3.2,0,2.8,0,0,0; 
                 2.9,0,2.5,0,0,0,0,0,0];
    link_cost (link_cost ==0) = 1e4;

    link_cost = max(0,link_cost + (1e6/(25*365.25)).*randn(n));

    
    %fixed cost for constructing stations
    station_cost = (1e6/(25*365.25)).*[2, 3, 2.2, 3, 2.5, 1.3, 2.8, 2.2, 3.1];

    station_cost = max(0,station_cost + (1e6/(25*365.25)).*randn(1,n));
    
    link_capacity_slope = 0.04.*link_cost; 
    station_capacity_slope = 0.04.*station_cost;
    
    %demand between od pairs
    demand = 1e3.*[0,9,26,19,13,12,13,8,11;
              11,0,14,26,7,18,3,6,12;
              30,19,0,30,24,8,15,12,5;
              21,9,11,0,22,16,25,21,23;
              14,14,8,9,0,20,16,22,21;
              26,1,22,24,13,0,16,14,12;
              8,6,9,23,6,13,0,11,11;
              9,2,14,20,18,16,11,0,4;
              8,7,11,22,27,17,8,12,0];

    demand = max(0,demand + 1e3.*3.*eye(n).*randn(n));
    
    distance = 10000 * ones(n, n); % Distances between arcs
    
    for i = 1:n
        distance(i, i) = 0;
    end
    
    distance(1, 2) = 0.75;
    distance(1, 3) = 0.7;
    distance(1, 9) = 0.9;
    
    distance(2, 3) = 0.6;
    distance(2, 4) = 1.1;
    
    distance(3, 4) = 1.1;
    distance(3, 5) = 0.5;
    distance(3, 9) = 0.7;
    
    distance(4,5) = 0.8;
    distance(4,6) = 0.7;
    distance(4,8) = 0.8;
    
    distance(5,6) = 0.5;
    distance(5,7) = 0.7;
    
    distance(6,7) = 0.5;
    distance(6,8) = 0.4;
    
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

    candidates = zeros(n);
    for i=1:n
        candidates(i,candidasourcertes{i}) = 1;
    end

    
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

function [sprim,s,deltas,aprim,a,deltaa,f,fext,fij,comp_time] = compute_sim(niters,lam,beta,alfa,n)

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


    fid = fopen("./export_txt/niters.txt",'w');
    if fid < 0, error('No puedo abrir %s', filename); end
        fprintf(fid, '%d', niters);
    fclose(fid);


    cmdpath = "cd C:\GAMS\50";

    system(cmdpath);
    
    gmsFile  = 'C:\Users\freal\MATLAB\Projects\untitled\code\9node\new1.gms';

    gamsExe = 'C:\GAMS\50\gams.exe'; 

    comp_time = 0;


    for iter=1:niters


        fid = fopen("./export_txt/current_iter.txt",'w');
        if fid < 0, error('No puedo abrir %s', filename); end
            fprintf(fid, '%d', iter);
        fclose(fid);
    
        cmd = sprintf('%s %s ', ...
                  gamsExe, gmsFile);
        [status,out] = system(cmd);
        disp(out);
        
        results_file_ctime = readtable('./output_all.xlsx','Sheet','solver_time');
        comp_time = comp_time + table2array(results_file_ctime);

        results_file_a = readtable('./output_all.xlsx','Sheet','aprim_level');
        a_prev = table2array(results_file_a(1:n,2:(n+1)));
        a_prev = max(1e-2,a_prev);
    
        results_file_s = readtable('./output_all.xlsx','Sheet','sprim_level');
        s_prev = table2array(results_file_s(1,:));
        s_prev = max(0,s_prev);
        if iter < niters
            write_gams_param_ii('./export_txt/a_prev.txt', a_prev);
            write_gams_param1d_full('./export_txt/s_prev.txt', s_prev);
        end
    
    end
    
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

end
