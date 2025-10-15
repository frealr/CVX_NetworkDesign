close all;
clear all;
clc;

% === 1) RUTA de tu GAMS ===
gamsHome = 'C:\GAMS\50';  % <-- CAMBIA esto

[basedir,~,~] = fileparts(mfilename('fullpath'));
basedir = fullfile(basedir, 'export_csv');   % subcarpeta de exportación
if ~exist(basedir,'dir'); mkdir(basedir); end

% === Tus datos ===
[n,link_cost,station_cost,link_capacity_slope, ...
  station_capacity_slope,demand,prices, ...
  op_link_cost,congestion_coef_stations, ...
  congestion_coef_links,travel_time,alt_time,alt_price, ...
  a_nom,tau,sigma, ...
  a_max,candidasourcertes,pond_coefs] = parameters_sevilla_network();

candidates = zeros(n);
for i=1:n
    candidates(i,candidasourcertes{i}) = 1;
end

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

% === 2D matrices ===
write_gams_param_ii('./export_txt/demand.txt', demand);
write_gams_param_ii('./export_txt/travel_time.txt', travel_time);
write_gams_param_ii('./export_txt/alt_time.txt', alt_time);
write_gams_param_ii('./export_txt/alt_price.txt', alt_price);
write_gams_param_ii('./export_txt/link_cost.txt', link_cost);
write_gams_param_ii('./export_txt/link_capacity_slope.txt', link_capacity_slope);
write_gams_param_ii('./export_txt/prices.txt', prices);
write_gams_param_ii('./export_txt/op_link_cost.txt', op_link_cost);
write_gams_param_ii('./export_txt/congestion_coefs_links.txt', congestion_coef_links);
write_gams_param_ii('./export_txt/candidates.txt', candidates);



% === 1D vectores ===
write_gams_param1d_full('./export_txt/station_cost.txt', station_cost);
write_gams_param1d_full('./export_txt/station_capacity_slope.txt', station_capacity_slope);
write_gams_param1d_full('./export_txt/congestion_coefs_stations.txt', congestion_coef_stations);





niters = 10;
alfa = 0.6;
n=24;

%% lam 5

betas = [1e-2,2.5e-2,5e-2,7.5e-2,1e-1:1e-2:3e-1,0.13:1e-3:0.15]; %betas for lam = 5
betas = 1e-1:2.5e-3:1.3e-1;
lams = [5];

for ll = 1:length(lams)
   lam = lams(ll);
   for bb=1:length(betas)
       beta = betas(bb);
       tic;
       [sprim,s,deltas,aprim,a,deltaa,f,fext,fij] = compute_sim(niters,lam,beta,alfa,n);
       comp_time = toc;
       [obj_val,pax_obj,op_obj] = get_obj_val(alfa, op_link_cost,...
            congestion_coef_links, ...
            congestion_coef_stations,travel_time,prices,alt_time,alt_price,aprim,deltaa, ...
            sprim,deltas,fij,f,fext,demand);
       budget = get_budget(s,sprim,a,aprim,n,...
            station_cost,station_capacity_slope,link_cost,link_capacity_slope,lam)
       filename = sprintf('./sevilla_rebutal/beta=%d_lam=%d.mat',beta,lam);
        save(filename,'s','sprim','deltas', ...
        'a','aprim','deltaa','f','fext','fij','comp_time','budget', ...
        'pax_obj','op_obj','obj_val');
   end

end


%% lam 10

%betas for lam=10
betas = [0,5e-2:1e-2:1e-1,6e-2:1e-3:8.2e-2];
lams = [10];
% 
% for ll = 1:length(lams)
%    lam = lams(ll);
%    for bb=1:length(betas)
%        beta = betas(bb);
%        tic;
%        [sprim,s,deltas,aprim,a,deltaa,f,fext,fij] = compute_sim(niters,lam,beta,alfa,n);
%        comp_time = toc;
%        [obj_val,pax_obj,op_obj] = get_obj_val(alfa, op_link_cost,...
%             congestion_coef_links, ...
%             congestion_coef_stations,travel_time,prices,alt_time,alt_price,aprim,deltaa, ...
%             sprim,deltas,fij,f,fext,demand);
%        budget = get_budget(s,sprim,a,aprim,n,...
%             station_cost,station_capacity_slope,link_cost,link_capacity_slope,lam)
%        filename = sprintf('./sevilla_rebutal/beta=%d_lam=%d.mat',beta,lam);
%         save(filename,'s','sprim','deltas', ...
%         'a','aprim','deltaa','f','fext','fij','comp_time','budget', ...
%         'pax_obj','op_obj','obj_val');
%    end
% 
% end


%% lam 15

betas = [1e-2:1e-2:7e-2,4e-2:1e-3:6e-2]; %betas for lam=15
lams = [15];

% 
% for ll = 1:length(lams)
%    lam = lams(ll);
%    for bb=1:length(betas)
%        beta = betas(bb);
%        tic;
%        [sprim,s,deltas,aprim,a,deltaa,f,fext,fij] = compute_sim(niters,lam,beta,alfa,n);
%        comp_time = toc;
%        [obj_val,pax_obj,op_obj] = get_obj_val(alfa, op_link_cost,...
%             congestion_coef_links, ...
%             congestion_coef_stations,travel_time,prices,alt_time,alt_price,aprim,deltaa, ...
%             sprim,deltas,fij,f,fext,demand);
%        budget = get_budget(s,sprim,a,aprim,n,...
%             station_cost,station_capacity_slope,link_cost,link_capacity_slope,lam)
%        filename = sprintf('./sevilla_rebutal/beta=%d_lam=%d.mat',beta,lam);
%         save(filename,'s','sprim','deltas', ...
%         'a','aprim','deltaa','f','fext','fij','comp_time','budget', ...
%         'pax_obj','op_obj','obj_val');
%    end
% 
% end






%% Load obtained results
lams = [5];
betas = [1e-3,5e-3,1e-2,5e-2,1e-1,3e-1,5e-1,7e-1,1];
betas = 1e-2:1e-2:7e-2;
betas = 4e-2:1e-3:6e-2;
betas = 0.046;
betas = [5e-2:1e-2:1e-1,6e-2:1e-3:8.2e-2];
betas = [1e-2,2.5e-2,5e-2,7.5e-2,1e-1:1e-2:3e-1,0.13:1e-3:0.15]; 
betas = [1e-2:1e-2:7e-2,4e-2:1e-3:6e-2]; %lam 15
betas = [0.04,0.042,0.043,0.044,0.045,0.046,0.047,0.048]; %present with lam15
%betas = [1e-2,5e-2:1e-2:1e-1,6e-2:1e-3:8.2e-2]; %lam 10
betas = [0.05,0.06,0.061,0.063,0.064,0.067,0.068,0.069]; %present with lam10
%betas = [1e-2,2.5e-2,5e-2,7.5e-2,1e-1:1e-2:3e-1,0.13:1e-3:0.15, 1e-1:2.5e-3:1.3e-1]; %lam 5
betas = [0.075,0.105,0.11,0.1125,0.115, 0.12, 0.1225, 0.125]; %present with lam5
%betas = 1e-2:2.5e-4:1.3e-2;
close all;

% to plot
% lams = [5];
% betas = [0.075, 0.105, 0.11,0.115];

% lams = [10];
% betas = [0.05, 0.06,0.063,0.064];

% lams = [15];
% betas = [0.04, 0.042, 0.043,0.045];



close all;
figure('Position', [100, 100, 1000, 400]);

azul_col = [0 0.4470 0.7410]; %icvx
rojo_col = [0.8500 0.3250 0.0980]; %mip30
naranja_col = [0.9290 0.6940 0.1250]; %MIP10
verde_col = [0.4660 0.6740 0.1880]; %mipreg






% xlim([0.5, nX+0.5]); xticks(1:nX)
% eur =['[',char(8364),']'];
% xl = xlabel(['$\beta$'],'interpreter','latex');
% yl = ylabel('Diff [\%]','Interpreter','latex');
% set(gca, 'FontSize', 9);
% set(gca, 'TickLabelInterpreter', 'latex');
% set(gca, 'XTick', 1:9, 'XTickLabel', betas);
% 
% % leyenda “dummy”
% plot(nan,nan,'Color',[0 0 0],'LineWidth',2); plot(nan,nan,'Color',[0.5 0.5 0.5],'LineWidth',2);
% legend({'ICVX w.r.t. MIPREG, $\lambda = 5$','ICVX w.r.t. MIPREG, $\lambda = 6$'},'Location','best','Interpreter','latex','FontSize',9); hold off
% xlim([0.7 9.3]);

tipos = {'s-','o-','*-'};

color_line = {azul_col,rojo_col,naranja_col};


for ll=1:length(lams)
    lam = lams(ll);
    for bb=1:length(betas)
        beta_or = betas(bb);
        [n,link_cost,station_cost,link_capacity_slope,...
        station_capacity_slope,demand,prices,...
        op_link_cost,congestion_coef_stations,...
        congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,sigma,...
        a_max,candidates,pond_coefs_mat] = parameters_sevilla_network();
      
        filename = sprintf('./sevilla_rebutal/beta=%d_lam=%d.mat',beta_or,lam);
        load(filename);
        att_d(bb) = 100*sum(sum(f.*demand))/sum(sum(demand));
        bud(bb) = budget;
        nroutes = 0;
        dis_rut = 0;
        for o=1:n
            for d=1:n
                dis(o,d) = sum(sum(travel_time.*squeeze(fij(:,:,o,d))))*demand(o,d);
                if f(o,d) > 0.01
                    dis_rut = dis_rut + sum(sum(travel_time(squeeze(fij(:,:,o,d)) > 0.02)));
                    nroutes = nroutes + 1;
                end
                uu(o,d) = demand(o,d)*alt_time(o,d)*fext(o,d);
            end
        end
        long_mean(bb) = dis_rut/nroutes;
        d_med(bb) = sum(sum(dis))/(sum(sum(f.*demand)));
        u_med(bb) = sum(sum(uu))/(sum(sum(fext.*demand)));
        att_dem = round(100*sum(sum(f.*demand))/sum(sum(demand)),1);
        n_links = sum(sum(aprim));
        n_nodes = sum(sprim > 0.1);
        served_demand(bb) = round(100*sum(sum(f.*demand))./sum(sum(demand)),1);
        % disp(['beta = ',num2str(beta_or), ', lam = ', num2str(lam)]);
        % disp(['att_dem = ',num2str(served_demand),' %']);
        % disp(['Pres. = ',num2str(budget)]);
        % disp(['Arcos = ',num2str(sum(sum(aprim > 0.1)))]);
        % disp(['Nodos = ',num2str(sum(sprim > 0.1))]);
        % disp(['Cap. = ',num2str(sum(sum(aprim)))]);
        % disp(['distancia por pasajero red nueva = ',num2str(d_med(bb))]);
        % disp(['distancia por pasajero red actual = ',num2str(u_med(bb))]);
        % disp(['distancia por ruta = ',num2str(long_mean(bb))]);
        % disp(['tiempo computacional = ',num2str(comp_time)]);
        % disp('\n');
        % plot_network(aprim,lam,beta_or);
         % disp([sprintf('%.3f',beta_or),'&', ...
         %     num2str(served_demand), ...
         %     '&',sprintf('%.3e',budget),'&',num2str(sum(sum(aprim > 0.1))), ...
         %     '&',num2str(n_nodes), ...
         %     '&',num2str(n_links), ...
         %     '&',sprintf('%.2f',d_med(bb)),'&',sprintf('%.2f',u_med(bb)),'&', ...
         %     sprintf('%.2f',long_mean(bb)),'&',sprintf('%.2e',comp_time),'\\ \hline']);
       
    end
    subplot(1,2,1);
    plot(bud,served_demand, tipos{1}, 'LineWidth',2,'Color',color_line{1}); hold on;
    subplot(1,2,2);
    plot(betas,served_demand,tipos{1}, 'LineWidth',2,'Color',color_line{1}); hold on;



end

lams = [10];
betas = [0.05,0.06,0.061,0.063,0.064,0.067,0.068,0.069];

for ll=1:length(lams)
    lam = lams(ll);
    for bb=1:length(betas)
        beta_or = betas(bb);
        [n,link_cost,station_cost,link_capacity_slope,...
        station_capacity_slope,demand,prices,...
        op_link_cost,congestion_coef_stations,...
        congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,sigma,...
        a_max,candidates,pond_coefs_mat] = parameters_sevilla_network();
      
        filename = sprintf('./sevilla_rebutal/beta=%d_lam=%d.mat',beta_or,lam);
        load(filename);
        att_d(bb) = 100*sum(sum(f.*demand))/sum(sum(demand));
        bud(bb) = budget;
        nroutes = 0;
        dis_rut = 0;
        for o=1:n
            for d=1:n
                dis(o,d) = sum(sum(travel_time.*squeeze(fij(:,:,o,d))))*demand(o,d);
                if f(o,d) > 0.01
                    dis_rut = dis_rut + sum(sum(travel_time(squeeze(fij(:,:,o,d)) > 0.02)));
                    nroutes = nroutes + 1;
                end
                uu(o,d) = demand(o,d)*alt_time(o,d)*fext(o,d);
            end
        end
        long_mean(bb) = dis_rut/nroutes;
        d_med(bb) = sum(sum(dis))/(sum(sum(f.*demand)));
        u_med(bb) = sum(sum(uu))/(sum(sum(fext.*demand)));
        att_dem = round(100*sum(sum(f.*demand))/sum(sum(demand)),1);
        n_links = sum(sum(aprim));
        n_nodes = sum(sprim > 0.1);
        served_demand(bb) = round(100*sum(sum(f.*demand))./sum(sum(demand)),1);
        % disp(['beta = ',num2str(beta_or), ', lam = ', num2str(lam)]);
        % disp(['att_dem = ',num2str(served_demand),' %']);
        % disp(['Pres. = ',num2str(budget)]);
        % disp(['Arcos = ',num2str(sum(sum(aprim > 0.1)))]);
        % disp(['Nodos = ',num2str(sum(sprim > 0.1))]);
        % disp(['Cap. = ',num2str(sum(sum(aprim)))]);
        % disp(['distancia por pasajero red nueva = ',num2str(d_med(bb))]);
        % disp(['distancia por pasajero red actual = ',num2str(u_med(bb))]);
        % disp(['distancia por ruta = ',num2str(long_mean(bb))]);
        % disp(['tiempo computacional = ',num2str(comp_time)]);
        % disp('\n');
        % plot_network(aprim,lam,beta_or);
         % disp([sprintf('%.3f',beta_or),'&', ...
         %     num2str(served_demand), ...
         %     '&',sprintf('%.3e',budget),'&',num2str(sum(sum(aprim > 0.1))), ...
         %     '&',num2str(n_nodes), ...
         %     '&',num2str(n_links), ...
         %     '&',sprintf('%.2f',d_med(bb)),'&',sprintf('%.2f',u_med(bb)),'&', ...
         %     sprintf('%.2f',long_mean(bb)),'&',sprintf('%.2e',comp_time),'\\ \hline']);
       
    end
    subplot(1,2,1);
    plot(bud,served_demand, tipos{2}, 'LineWidth',2,'Color',color_line{2}); hold on;
    subplot(1,2,2);
    plot(betas,served_demand,tipos{2}, 'LineWidth',2,'Color',color_line{2}); hold on;


end

lams = [15];
betas = [0.04,0.042,0.043,0.044,0.046,0.047,0.048];
clear bud served_demand;

for ll=1:length(lams)
    lam = lams(ll);
    for bb=1:length(betas)
        beta_or = betas(bb);
        [n,link_cost,station_cost,link_capacity_slope,...
        station_capacity_slope,demand,prices,...
        op_link_cost,congestion_coef_stations,...
        congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,sigma,...
        a_max,candidates,pond_coefs_mat] = parameters_sevilla_network();
      
        filename = sprintf('./sevilla_rebutal/beta=%d_lam=%d.mat',beta_or,lam);
        load(filename);
        att_d(bb) = 100*sum(sum(f.*demand))/sum(sum(demand));
        bud(bb) = budget;
        nroutes = 0;
        dis_rut = 0;
        for o=1:n
            for d=1:n
                dis(o,d) = sum(sum(travel_time.*squeeze(fij(:,:,o,d))))*demand(o,d);
                if f(o,d) > 0.01
                    dis_rut = dis_rut + sum(sum(travel_time(squeeze(fij(:,:,o,d)) > 0.02)));
                    nroutes = nroutes + 1;
                end
                uu(o,d) = demand(o,d)*alt_time(o,d)*fext(o,d);
            end
        end
        long_mean(bb) = dis_rut/nroutes;
        d_med(bb) = sum(sum(dis))/(sum(sum(f.*demand)));
        u_med(bb) = sum(sum(uu))/(sum(sum(fext.*demand)));
        att_dem = round(100*sum(sum(f.*demand))/sum(sum(demand)),1);
        n_links = sum(sum(aprim));
        n_nodes = sum(sprim > 0.1);
        served_demand(bb) = round(100*sum(sum(f.*demand))./sum(sum(demand)),1);
        % disp(['beta = ',num2str(beta_or), ', lam = ', num2str(lam)]);
        % disp(['att_dem = ',num2str(served_demand),' %']);
        % disp(['Pres. = ',num2str(budget)]);
        % disp(['Arcos = ',num2str(sum(sum(aprim > 0.1)))]);
        % disp(['Nodos = ',num2str(sum(sprim > 0.1))]);
        % disp(['Cap. = ',num2str(sum(sum(aprim)))]);
        % disp(['distancia por pasajero red nueva = ',num2str(d_med(bb))]);
        % disp(['distancia por pasajero red actual = ',num2str(u_med(bb))]);
        % disp(['distancia por ruta = ',num2str(long_mean(bb))]);
        % disp(['tiempo computacional = ',num2str(comp_time)]);
        % disp('\n');
        % plot_network(aprim,lam,beta_or);
         % disp([sprintf('%.3f',beta_or),'&', ...
         %     num2str(served_demand), ...
         %     '&',sprintf('%.3e',budget),'&',num2str(sum(sum(aprim > 0.1))), ...
         %     '&',num2str(n_nodes), ...
         %     '&',num2str(n_links), ...
         %     '&',sprintf('%.2f',d_med(bb)),'&',sprintf('%.2f',u_med(bb)),'&', ...
         %     sprintf('%.2f',long_mean(bb)),'&',sprintf('%.2e',comp_time),'\\ \hline']);
       
    end
    subplot(1,2,1);
    plot(bud,served_demand, tipos{3}, 'LineWidth',2,'Color',color_line{3}); hold on;
    subplot(1,2,2);
    plot(betas,served_demand,tipos{3}, 'LineWidth',2,'Color',color_line{3}); hold on;


end

subplot(1,2,1);

xl = xlabel(['Budget [EUR]'],'interpreter','latex');
yl = ylabel('SD [\%]','Interpreter','latex');
set(gca, 'FontSize', 9);
set(gca, 'TickLabelInterpreter', 'latex');
grid on; legend({'$\lambda = 5$','$\lambda = 10$','$\lambda = 15$'},'Location','best','Interpreter','latex','FontSize',9);
subplot(1,2,2);
xl = xlabel(['$\beta$'],'interpreter','latex');
yl = ylabel('SD [\%]','Interpreter','latex');
set(gca, 'FontSize', 9);
set(gca, 'TickLabelInterpreter', 'latex');
grid on; legend({'$\lambda = 5$','$\lambda = 10$','$\lambda = 15$'},'Location','best','Interpreter','latex','FontSize',9);

saveas(gcf, './sevilla_dem_pres.png');


%%
lams = [5];
betas = [1e-3,5e-3,1e-2,5e-2,1e-1,3e-1,5e-1,7e-1,1];
betas = 1e-2:1e-2:7e-2;
betas = 4e-2:1e-3:6e-2;
betas = 0.046;
betas = [5e-2:1e-2:1e-1,6e-2:1e-3:8.2e-2];
betas = [1e-2,2.5e-2,5e-2,7.5e-2,1e-1:1e-2:3e-1,0.13:1e-3:0.15]; 
betas = [1e-2:1e-2:7e-2,4e-2:1e-3:6e-2]; %lam 15
betas = [0.04,0.042,0.043,0.044,0.045,0.046,0.047,0.048]; %present with lam15
%betas = [1e-2,5e-2:1e-2:1e-1,6e-2:1e-3:8.2e-2]; %lam 10
betas = [0.05,0.06,0.061,0.063,0.064,0.067,0.068,0.069]; %present with lam10
%betas = [1e-2,2.5e-2,5e-2,7.5e-2,1e-1:1e-2:3e-1,0.13:1e-3:0.15, 1e-1:2.5e-3:1.3e-1]; %lam 5
betas = [0.075,0.105,0.11,0.1125,0.115, 0.12, 0.1225, 0.125]; %present with lam5
%betas = 1e-2:2.5e-4:1.3e-2;
close all;

% to plot
% lams = [5];
% betas = [0.075, 0.105, 0.11,0.115];

% lams = [10];
% betas = [0.05, 0.06,0.063,0.064];

% lams = [15];
% betas = [0.04, 0.042, 0.043,0.045];



close all;
figure('Position', [100, 100, 1000, 800]);

azul_col = [0 0.4470 0.7410]; %icvx
rojo_col = [0.8500 0.3250 0.0980]; %mip30
naranja_col = [0.9290 0.6940 0.1250]; %MIP10
verde_col = [0.4660 0.6740 0.1880]; %mipreg






% xlim([0.5, nX+0.5]); xticks(1:nX)
% eur =['[',char(8364),']'];
% xl = xlabel(['$\beta$'],'interpreter','latex');
% yl = ylabel('Diff [\%]','Interpreter','latex');
% set(gca, 'FontSize', 9);
% set(gca, 'TickLabelInterpreter', 'latex');
% set(gca, 'XTick', 1:9, 'XTickLabel', betas);
% 
% % leyenda “dummy”
% plot(nan,nan,'Color',[0 0 0],'LineWidth',2); plot(nan,nan,'Color',[0.5 0.5 0.5],'LineWidth',2);
% legend({'ICVX w.r.t. MIPREG, $\lambda = 5$','ICVX w.r.t. MIPREG, $\lambda = 6$'},'Location','best','Interpreter','latex','FontSize',9); hold off
% xlim([0.7 9.3]);

tipos = {'s-','o-','*-'};

color_line = {azul_col,rojo_col,naranja_col};


for ll=1:length(lams)
    lam = lams(ll);
    for bb=1:length(betas)
        beta_or = betas(bb);
        [n,link_cost,station_cost,link_capacity_slope,...
        station_capacity_slope,demand,prices,...
        op_link_cost,congestion_coef_stations,...
        congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,sigma,...
        a_max,candidates,pond_coefs_mat] = parameters_sevilla_network();
      
        filename = sprintf('./sevilla_rebutal/beta=%d_lam=%d.mat',beta_or,lam);
        load(filename);
        att_d(bb) = 100*sum(sum(f.*demand))/sum(sum(demand));
        bud(bb) = budget;
        nroutes = 0;
        dis_rut = 0;
        for o=1:n
            for d=1:n
                dis(o,d) = sum(sum(travel_time.*squeeze(fij(:,:,o,d))))*demand(o,d);
                if f(o,d) > 0.01
                    dis_rut = dis_rut + sum(sum(travel_time(squeeze(fij(:,:,o,d)) > 0.02)));
                    nroutes = nroutes + 1;
                end
                uu(o,d) = demand(o,d)*alt_time(o,d)*fext(o,d);
            end
        end
        long_mean(bb) = dis_rut/nroutes;
        d_med(bb) = sum(sum(dis))/(sum(sum(f.*demand)));
        u_med(bb) = sum(sum(uu))/(sum(sum(fext.*demand)));
        att_dem = round(100*sum(sum(f.*demand))/sum(sum(demand)),1);
        n_links = sum(sum(aprim));
        n_nodes = sum(sprim > 0.1);
        served_demand(bb) = round(100*sum(sum(f.*demand))./sum(sum(demand)),1);
        % disp(['beta = ',num2str(beta_or), ', lam = ', num2str(lam)]);
        % disp(['att_dem = ',num2str(served_demand),' %']);
        % disp(['Pres. = ',num2str(budget)]);
        % disp(['Arcos = ',num2str(sum(sum(aprim > 0.1)))]);
        % disp(['Nodos = ',num2str(sum(sprim > 0.1))]);
        % disp(['Cap. = ',num2str(sum(sum(aprim)))]);
        % disp(['distancia por pasajero red nueva = ',num2str(d_med(bb))]);
        % disp(['distancia por pasajero red actual = ',num2str(u_med(bb))]);
        % disp(['distancia por ruta = ',num2str(long_mean(bb))]);
        % disp(['tiempo computacional = ',num2str(comp_time)]);
        % disp('\n');
        % plot_network(aprim,lam,beta_or);
         % disp([sprintf('%.3f',beta_or),'&', ...
         %     num2str(served_demand), ...
         %     '&',sprintf('%.3e',budget),'&',num2str(sum(sum(aprim > 0.1))), ...
         %     '&',num2str(n_nodes), ...
         %     '&',num2str(n_links), ...
         %     '&',sprintf('%.2f',d_med(bb)),'&',sprintf('%.2f',u_med(bb)),'&', ...
         %     sprintf('%.2f',long_mean(bb)),'&',sprintf('%.2e',comp_time),'\\ \hline']);
       
    end
    subplot(2,3,1);
    plot(betas,d_med, tipos{1}, 'LineWidth',2,'Color',color_line{1}); hold on;
    plot(betas,u_med, tipos{2}, 'LineWidth',2,'Color',color_line{2});
    xl = xlabel(['$\beta$'],'interpreter','latex');
    yl = ylabel('$\bar{t}_{PAX}$ [min]','Interpreter','latex');
    set(gca, 'FontSize', 9);
    set(gca, 'TickLabelInterpreter', 'latex');
    grid on; legend({'New net, $\lambda = 5$','Alternative net, $\lambda = 5$'},'Location','best','Interpreter','latex','FontSize',9);

    subplot(2,3,4);
    plot(betas,long_mean, tipos{1}, 'LineWidth',2,'Color',color_line{1}); hold on;
    xl = xlabel(['$\beta$'],'interpreter','latex');
    yl = ylabel('$\bar{t}_{route}$ [min]','Interpreter','latex');
    set(gca, 'FontSize', 9);
    set(gca, 'TickLabelInterpreter', 'latex');
    grid on; legend({'New net, $\lambda = 5$'},'Location','best','Interpreter','latex','FontSize',9);


end

lams = [10];
betas = [0.05,0.06,0.061,0.063,0.064,0.067,0.068,0.069];

for ll=1:length(lams)
    lam = lams(ll);
    for bb=1:length(betas)
        beta_or = betas(bb);
        [n,link_cost,station_cost,link_capacity_slope,...
        station_capacity_slope,demand,prices,...
        op_link_cost,congestion_coef_stations,...
        congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,sigma,...
        a_max,candidates,pond_coefs_mat] = parameters_sevilla_network();
      
        filename = sprintf('./sevilla_rebutal/beta=%d_lam=%d.mat',beta_or,lam);
        load(filename);
        att_d(bb) = 100*sum(sum(f.*demand))/sum(sum(demand));
        bud(bb) = budget;
        nroutes = 0;
        dis_rut = 0;
        for o=1:n
            for d=1:n
                dis(o,d) = sum(sum(travel_time.*squeeze(fij(:,:,o,d))))*demand(o,d);
                if f(o,d) > 0.01
                    dis_rut = dis_rut + sum(sum(travel_time(squeeze(fij(:,:,o,d)) > 0.02)));
                    nroutes = nroutes + 1;
                end
                uu(o,d) = demand(o,d)*alt_time(o,d)*fext(o,d);
            end
        end
        long_mean(bb) = dis_rut/nroutes;
        d_med(bb) = sum(sum(dis))/(sum(sum(f.*demand)));
        u_med(bb) = sum(sum(uu))/(sum(sum(fext.*demand)));
        att_dem = round(100*sum(sum(f.*demand))/sum(sum(demand)),1);
        n_links = sum(sum(aprim));
        n_nodes = sum(sprim > 0.1);
        served_demand(bb) = round(100*sum(sum(f.*demand))./sum(sum(demand)),1);
        % disp(['beta = ',num2str(beta_or), ', lam = ', num2str(lam)]);
        % disp(['att_dem = ',num2str(served_demand),' %']);
        % disp(['Pres. = ',num2str(budget)]);
        % disp(['Arcos = ',num2str(sum(sum(aprim > 0.1)))]);
        % disp(['Nodos = ',num2str(sum(sprim > 0.1))]);
        % disp(['Cap. = ',num2str(sum(sum(aprim)))]);
        % disp(['distancia por pasajero red nueva = ',num2str(d_med(bb))]);
        % disp(['distancia por pasajero red actual = ',num2str(u_med(bb))]);
        % disp(['distancia por ruta = ',num2str(long_mean(bb))]);
        % disp(['tiempo computacional = ',num2str(comp_time)]);
        % disp('\n');
        % plot_network(aprim,lam,beta_or);
         % disp([sprintf('%.3f',beta_or),'&', ...
         %     num2str(served_demand), ...
         %     '&',sprintf('%.3e',budget),'&',num2str(sum(sum(aprim > 0.1))), ...
         %     '&',num2str(n_nodes), ...
         %     '&',num2str(n_links), ...
         %     '&',sprintf('%.2f',d_med(bb)),'&',sprintf('%.2f',u_med(bb)),'&', ...
         %     sprintf('%.2f',long_mean(bb)),'&',sprintf('%.2e',comp_time),'\\ \hline']);
       
    end
    subplot(2,3,2);
    plot(betas,d_med, tipos{1}, 'LineWidth',2,'Color',color_line{1}); hold on;
    plot(betas,u_med, tipos{2}, 'LineWidth',2,'Color',color_line{2});
    xl = xlabel(['$\beta$'],'interpreter','latex');
    yl = ylabel('$\bar{t}_{PAX}$ [min]','Interpreter','latex');
    set(gca, 'FontSize', 9);
    set(gca, 'TickLabelInterpreter', 'latex');
    grid on; legend({'New net, $\lambda = 10$','Alternative net, $\lambda = 10$'},'Location','best','Interpreter','latex','FontSize',9);

    subplot(2,3,5);
    plot(betas,long_mean, tipos{1}, 'LineWidth',2,'Color',color_line{1}); hold on;
    xl = xlabel(['$\beta$'],'interpreter','latex');
    yl = ylabel('$\bar{t}_{route}$ [min]','Interpreter','latex');
    set(gca, 'FontSize', 9);
    set(gca, 'TickLabelInterpreter', 'latex');
    grid on; legend({'New net, $\lambda = 10$'},'Location','best','Interpreter','latex','FontSize',9);


end

lams = [15];
betas = [0.04,0.042,0.043,0.044,0.046,0.047,0.048];
clear bud served_demand d_med u_med long_mean;

for ll=1:length(lams)
    lam = lams(ll);
    for bb=1:length(betas)
        beta_or = betas(bb);
        [n,link_cost,station_cost,link_capacity_slope,...
        station_capacity_slope,demand,prices,...
        op_link_cost,congestion_coef_stations,...
        congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,sigma,...
        a_max,candidates,pond_coefs_mat] = parameters_sevilla_network();
      
        filename = sprintf('./sevilla_rebutal/beta=%d_lam=%d.mat',beta_or,lam);
        load(filename);
        att_d(bb) = 100*sum(sum(f.*demand))/sum(sum(demand));
        bud(bb) = budget;
        nroutes = 0;
        dis_rut = 0;
        for o=1:n
            for d=1:n
                dis(o,d) = sum(sum(travel_time.*squeeze(fij(:,:,o,d))))*demand(o,d);
                if f(o,d) > 0.01
                    dis_rut = dis_rut + sum(sum(travel_time(squeeze(fij(:,:,o,d)) > 0.02)));
                    nroutes = nroutes + 1;
                end
                uu(o,d) = demand(o,d)*alt_time(o,d)*fext(o,d);
            end
        end
        long_mean(bb) = dis_rut/nroutes;
        d_med(bb) = sum(sum(dis))/(sum(sum(f.*demand)));
        u_med(bb) = sum(sum(uu))/(sum(sum(fext.*demand)));
        att_dem = round(100*sum(sum(f.*demand))/sum(sum(demand)),1);
        n_links = sum(sum(aprim));
        n_nodes = sum(sprim > 0.1);
        served_demand(bb) = round(100*sum(sum(f.*demand))./sum(sum(demand)),1);
        % disp(['beta = ',num2str(beta_or), ', lam = ', num2str(lam)]);
        % disp(['att_dem = ',num2str(served_demand),' %']);
        % disp(['Pres. = ',num2str(budget)]);
        % disp(['Arcos = ',num2str(sum(sum(aprim > 0.1)))]);
        % disp(['Nodos = ',num2str(sum(sprim > 0.1))]);
        % disp(['Cap. = ',num2str(sum(sum(aprim)))]);
        % disp(['distancia por pasajero red nueva = ',num2str(d_med(bb))]);
        % disp(['distancia por pasajero red actual = ',num2str(u_med(bb))]);
        % disp(['distancia por ruta = ',num2str(long_mean(bb))]);
        % disp(['tiempo computacional = ',num2str(comp_time)]);
        % disp('\n');
        % plot_network(aprim,lam,beta_or);
         % disp([sprintf('%.3f',beta_or),'&', ...
         %     num2str(served_demand), ...
         %     '&',sprintf('%.3e',budget),'&',num2str(sum(sum(aprim > 0.1))), ...
         %     '&',num2str(n_nodes), ...
         %     '&',num2str(n_links), ...
         %     '&',sprintf('%.2f',d_med(bb)),'&',sprintf('%.2f',u_med(bb)),'&', ...
         %     sprintf('%.2f',long_mean(bb)),'&',sprintf('%.2e',comp_time),'\\ \hline']);
       
    end
    subplot(2,3,3);
    plot(betas,d_med, tipos{1}, 'LineWidth',2,'Color',color_line{1}); hold on;
    plot(betas,u_med, tipos{2}, 'LineWidth',2,'Color',color_line{2});
    xl = xlabel(['$\beta$'],'interpreter','latex');
    yl = ylabel('$\bar{t}_{PAX}$ [min]','Interpreter','latex');
    set(gca, 'FontSize', 9);
    set(gca, 'TickLabelInterpreter', 'latex');
    grid on; legend({'New net, $\lambda = 15$','Alternative net, $\lambda = 15$'},'Location','best','Interpreter','latex','FontSize',9);

    subplot(2,3,6);
    plot(betas,long_mean, tipos{1}, 'LineWidth',2,'Color',color_line{1}); hold on;
    xl = xlabel(['$\beta$'],'interpreter','latex');
    yl = ylabel('$\bar{t}_{route}$ [min]','Interpreter','latex');
    set(gca, 'FontSize', 9);
    set(gca, 'TickLabelInterpreter', 'latex');
    grid on; legend({'New net, $\lambda = 15$'},'Location','best','Interpreter','latex','FontSize',9);


end

saveas(gcf, './sevilla_mean_timeroute_split.png');
%%
close all;
 plot_network(zeros(24),0,0);
%% Functions

function budget = get_budget(s,s_prim,a,a_prim,n,...
    station_cost,station_capacity_slope,link_cost,link_capacity_slope,lam)
    budget = 0;
    for i=1:n
        if s_prim(i) > 1e-3
            budget = budget + lam*station_cost(i) + ...
                station_capacity_slope(i)*s_prim(i);
        end
        for j=1:n
            if a_prim(i,j) > 1e-3
                budget = budget + lam*link_cost(i,j) + ...
                    link_capacity_slope(i,j) * a_prim(i,j);
            end
        end
    end
end

function [obj_val,pax_obj,op_obj] = get_obj_val(alfa, op_link_cost,...
    congestion_coef_links, ...
    congestion_coef_stations,travel_time,prices,alt_time,alt_price,a_prim,delta_a, ...
    s_prim,delta_s,fij,f,fext,demand)
    n = 24;
    logit_coef = 0.2;
    pax_obj = 0;
    op_obj = 0;
    eps = 1e-3;
    op_obj = op_obj + 1e-3*(sum(sum(op_link_cost.*a_prim))); %operational costs
    for i=1:n
        if s_prim(i) > eps
            pax_obj = pax_obj - 1e-3*sum(log(congestion_coef_stations(i)*delta_s(i) + eps));
        end

        for j=1:n
            if a_prim(i,j) > eps
                pax_obj = pax_obj - 1e-3*sum(sum(log(congestion_coef_links(i,j)*delta_a(i,j) + eps)));
            end
        end
    end
    for o=1:n
        for d=1:n
            for i=1:n
                for j=1:n
                    pax_obj = pax_obj + 1e-3*(demand(o,d).*logit_coef*(travel_time(i,j)+prices(i,j)).*fij(i,j,o,d));
                end
            end
        end
    end
    pax_obj = pax_obj + 1e-3*(sum(sum(demand.*logit_coef*(alt_time+alt_price).*fext)));

    entro = max(f.*(log(f+eps)-1),0);
    entro_ext = max(fext.*(log(fext+eps)-1),0);
    pax_obj = pax_obj + 1e-3*(sum(sum(demand.*(entro))));
    pax_obj = pax_obj + 1e-3*(sum(sum(demand.*(entro_ext))));
    obj_val = (alfa/(1))*pax_obj + ((1-alfa)/(1))*op_obj;
end



function [sprim,s,deltas,aprim,a,deltaa,f,fext,fij] = compute_sim(niters,lam,beta,alfa,n)

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
    
    gmsFile  = 'C:\Users\freal\MATLAB\Projects\untitled\code\sevilla_net\new1.gms';

    for iter=1:niters


        fid = fopen("./export_txt/current_iter.txt",'w');
        if fid < 0, error('No puedo abrir %s', filename); end
            fprintf(fid, '%d', iter);
        fclose(fid);
    
        cmd = sprintf("gams %s", ...
                      gmsFile);
        
        [status,stdoutText] = system(cmd);
        
        fprintf('Exit code: %d\n', status);          % 0 = OK
        disp(stdoutText);                             % log por consola
        
        results_file_a = readtable('./output_all.xlsx','Sheet','aprim_level');
        a_prev = table2array(results_file_a(1:n,2:(n+1)));
    
        results_file_s = readtable('./output_all.xlsx','Sheet','sprim_level');
        s_prev = table2array(results_file_s(1,:));
    
        write_gams_param_ii('./export_txt/a_prev.txt', a_prev);
        write_gams_param1d_full('./export_txt/s_prev.txt', s_prev);
    
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

function [n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    op_link_cost,congestion_coef_stations,...
    congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,sigma,...
    a_max,candidates,pond_coefs_mat] = parameters_sevilla_network()

    n = 24;
    
    %candidates to construct a link for each neighbor
    candidates = {
    [4,23,2,18,3,9,8,10,24,5,17];
    [1,4,23,18,3,9,8,5,24];
    [1,2,18,16,14,21,19,9,8,24];
    [1,2,23,16,15,20,11,22,6];
    [1,2,10,9,21,14,17,12,13,24];
    [4,22,16,14,19,7,20];
    [12,13,17,19,21,14,11,22,6,20];
    [1,2,3,9,10,24];
    [1,2,3,16,21,20,19,17,5,10,8,12];
    [1,8,9,20,19,12,17,5,24,21];
    [4,22,20,7,19,14,16,15,23];
    [13,5,17,10,9,19,14,7,20];
    [24,5,17,19,7,12];
    [21,3,18,16,11,22,6,20,7,12,19,17,5,15];
    [18,23,4,22,11,20,14,16];
    [18,23,4,15,11,6,14,21,9,3,22];
    [5,24,1,10,9,21,14,19,7,12,13];
    [2,23,15,16,14,21,19,3,1];
    [13,17,10,9,3,18,21,14,11,6,20,7,12];
    [6,22,11,4,15,14,21,9,10,19,7,12];
    [3,18,23,16,14,20,7,19,17,5,10,9];
    [4,6,20,7,14,16,11,15,23];
    [1,4,22,11,15,16,21,18,24,2];
    [1,2,23,3,8,10,17,5,13];
    };
    

    population_file = readtable('./population.xlsx');
    population = table2array(population_file(1:24,2));
    coordinates = readtable('./coordenadas_Sevilla.xlsx');
    coor_x = table2array(coordinates(1:24,3));
    coor_y = table2array(coordinates(1:24,7));
    rng(1,"twister"); %seed
    distance = 1e6.*ones(n);
    for i=1:n
        distance(i,i) = 0;
        cand = candidates(i);
        cand = cand{1};
        cand = cand(cand > i);
        for j=i+1:n
            if sum(j == cand) > 0
                distance(i,j) = haversine(coor_y(i), coor_x(i), coor_y(j), coor_x(j));
                distance(j,i) = distance(i,j);
            end
            non_stop = rand < 0.4;

            alt_cost(i,j) = haversine(coor_y(i), coor_x(i), coor_y(j), coor_x(j));
            alt_cost(i,j) = alt_cost(i,j) + 0.2*non_stop*alt_cost(i,j);
            alt_cost(j,i) = alt_cost(i,j);
        end
    end
    demand = [0, 272, 272, 272, 272, 553, 272, 272, 272, 553, 553, 272, 553, 272, 553, 553, 553, 272, 937, 937, 937, 937, 937, 937;
           272, 0, 272, 272, 272, 553, 272, 272, 272, 553, 553, 272, 553, 272, 553, 553, 553, 272, 937, 937, 937, 937, 937, 937;
           327, 327, 0, 327, 327, 664, 327, 327, 327, 664, 664, 327, 664, 327, 664, 664, 664, 327, 1125, 1125, 1125, 1125, 1125, 1125;
           185, 185, 185, 0, 185, 376, 185, 185, 185, 376, 376, 185, 376, 185, 376, 376, 376, 185, 637, 637, 637, 637, 637, 637;
           272, 272, 272, 272, 0, 553, 272, 272, 272, 553, 553, 272, 553, 272, 553, 553, 553, 272, 937, 937, 937, 937, 937, 937;
           225, 225, 225, 225, 225, 0, 225, 225, 225, 188, 188, 225, 188, 225, 188, 188, 188, 225, 284, 284, 284, 284, 284, 284;
           283, 283, 283, 283, 283, 575, 0, 283, 283, 575, 575, 283, 575, 283, 575, 575, 575, 283, 975, 975, 975, 975, 975, 975;
           272, 272, 272, 272, 272, 553, 272, 0, 272, 553, 553, 272, 553, 272, 553, 553, 553, 272, 937, 937, 937, 937, 937, 937;
           272, 272, 272, 272, 272, 553, 272, 272, 0, 553, 553, 272, 553, 272, 553, 553, 553, 272, 937, 937, 937, 937, 937, 937;
           511, 511, 511, 511, 511, 428, 511, 511, 511, 0, 428, 511, 428, 511, 428, 428, 428, 511, 645, 645, 645, 645, 645, 645;
           225, 225, 225, 225, 225, 188, 225, 225, 225, 188, 0, 225, 188, 225, 188, 188, 188, 225, 284, 284, 284, 284, 284, 284;
           272, 272, 272, 272, 272, 553, 272, 272, 272, 553, 553, 0, 553, 272, 553, 553, 553, 272, 937, 937, 937, 937, 937, 937;
           306, 306, 306, 306, 306, 257, 306, 306, 306, 257, 306, 257, 0, 257, 306, 306, 306, 257, 387, 387, 387, 387, 387, 387;
           294, 294, 294, 294, 294, 597, 294, 294, 294, 597, 597, 294, 597, 0, 597, 597, 597, 297, 1012, 1012, 1012, 1012, 1012, 1012;
           409, 409, 409, 409, 409, 342, 409, 409, 409, 342, 342, 409, 342, 409, 0, 342, 342, 409, 516, 516, 516, 516, 516, 516;
           511, 511, 511, 511, 511, 428, 511, 511, 511, 428, 428, 511, 428, 511, 428, 0, 428, 511, 645, 645, 645, 645, 645, 645;
           429, 429, 429, 429, 429, 360, 429, 429, 429, 306, 360, 429, 360, 429, 360, 360, 0, 429, 542, 542, 542, 542, 542, 542;
           272, 272, 272, 272, 272, 553, 272, 272, 272, 553, 553, 272, 553, 272, 553, 553, 553, 0, 937, 937, 937, 937, 937, 937;
           675, 675, 675, 675, 675, 730, 675, 675, 675, 730, 730, 675, 730, 675, 730, 730, 730, 675, 660, 660, 660, 660, 660, 660;
           879, 879, 879, 879, 879, 952, 879, 879, 879, 952, 952, 879, 952, 879, 952, 952, 952, 879, 860, 0, 860, 860, 860, 860;
           715, 715, 715, 715, 715, 775, 715, 715, 715, 775, 775, 715, 775, 715, 775, 775, 775, 715, 700, 700, 700, 700, 700, 700;
           511, 511, 511, 511, 511, 553, 511, 511, 511, 553, 553, 511, 553, 511, 553, 553, 553, 511, 500, 500, 500, 0, 500, 500;
           675, 675, 675, 675, 675, 730, 675, 675, 675, 730, 730, 675, 730, 675, 730, 730, 730, 675, 660, 660, 660, 660, 0, 660;
           675, 675, 675, 675, 675, 730, 675, 675, 675, 730, 730, 675, 730, 675, 730, 730, 730, 675, 660, 660, 660, 660, 660, 0];

    pond_coefs = [1.54535294, 1.54535294, 1.11218247, 1.72777031, 2.75490793, ...
        2.14294615, 1.41994992, 1, 1.11218247, 1.41994992,...
           2.14294615, 2.35701223, 2.75490793, 2.21242079, 3.        ,...
           1.59083705, 1.41994992, 1.01231062, 1.41994992, 2.35701223,...
           1.11218247, 2.14294615, 1.72777031, 1.45204153]';
    
    pond_coefs_tens(1,:,:) = pond_coefs.*ones(1,n);
    pond_coefs_tens(2,:,:) = permute(pond_coefs_tens(1,:,:),[1 3 2]);
    pond_coefs_mat = squeeze(permute(max(pond_coefs_tens(1,:,:),pond_coefs_tens(2,:,:)),[2,3,1]));

    crec_coefs = [1.6, 1.6, 1.1469802107427398, 0.9, 1.2081587290019313, 0.9661783579052806, 1.0802156586966714, ...
        0.9, 1.1469802107427398, 1.0802156586966714, 0.9661783579052806, 0.9989307690507944,...
        1.2081587290019313, 0.9, 1.0730507219686063, 0.9, 1.0802156586966714, 1.3414585127763425, ...
        1.0802156586966714, 0.9989307690507944, 1.1469802107427398, 0.9661783579052806, 0.9, 1.1224800389355853];

    % for o=1:n
    %     for d=1:n
    %         demand(o,d) = crec_coefs(o)*crec_coefs(d)*demand(o,d);
    %     end
    % end

   % demand = demand.*pond_coefs_mat;
    %demand = max(demand,demand');
    
    %fixed cost for constructing links

    link_cost = 1e6.*distance./(365.25*25);
    group_1 = [3,8,9,10,14,16,18,21];
    group_2 = [1,2,7,15,17,19];
    group_3 = [4,5,6,11,12,13,20,22,23,24];

    % 1-1: 3* (1.5 el 1, 1 el 2, 0.5 el 3)
    % 1-2: 2.5*
    % 1-3: 2*
    % 2-2: 2*
    % 2-3: 1.5*
    % 3-3: 1*

    for i=1:n
        if ismember(i,group_1)
            ci = 1.5;
        end
        if ismember(i,group_2)
            ci = 1;
        end
        if ismember(i,group_3)
            ci = 0.5;
        end

        for j=i+1:n
            if ismember(j,group_1)
                cj = 1.5;
            end
            if ismember(j,group_2)
                cj = 1;
            end
            if ismember(j,group_3)
                cj = 0.5;
            end
            link_cost(i,j) = link_cost(i,j).*(ci+cj);
            link_cost(j,i) = link_cost(i,j);
        end

    end





    nodos_centricos = [3,8,9,10,14,16,17,18,19,21];
    % for i=1:(length(nodos_centricos)-1)
    %     for j=i+1:length(nodos_centricos)
    %         link_cost(nodos_centricos(i),nodos_centricos(j)) = link_cost(nodos_centricos(i),nodos_centricos(j)) + 1e6.*5./(365.25.*25);
    %         link_cost(nodos_centricos(j),nodos_centricos(i)) =link_cost(nodos_centricos(i),nodos_centricos(j));
    %     end
    % end
    
    %fixed cost for constructing stations
   
  

    river_plus = 100;
   % river_plus = 0;


    link_cost(12,[17,10,9,19,14,7,20]) =  link_cost(12,[17,10,9,19,14,7,20]) + river_plus;
    link_cost([17,10,9,19,14,7,20],12) = link_cost([17,10,9,19,14,7,20],12) + river_plus;

    link_cost(5,[1,2,10,9,21,14,17]) = link_cost(5,[1,2,10,9,21,14,17]) + river_plus;
    link_cost([1,2,10,9,21,14,17],5) = link_cost([1,2,10,9,21,14,17],5) + river_plus;

    link_cost(13,[17,19,7]) = link_cost(13,[17,19,7]) + river_plus;
    link_cost([17,19,7],13) = link_cost([17,19,7],13) + river_plus;

    link_cost(24,[1,2,23,3,8,10,17]) = link_cost(24,[1,2,23,3,8,10,17]) + river_plus;
    link_cost([1,2,23,3,8,10,17],24) = link_cost([1,2,23,3,8,10,17],24) + river_plus;

    link_cost(7,[21,14,11,22,6,20]) = link_cost(7,[21,14,11,22,6,20]) + river_plus;
    link_cost([21,14,11,22,6,20],7) = link_cost([21,14,11,22,6,20],7) + river_plus;

    link_cost(19,[9,3,21,18,14,11]) = link_cost(19,[9,3,21,18,14,11]) + river_plus;
    link_cost([9,3,21,18,14,11],19) = link_cost([9,3,21,18,14,11],19) + river_plus;

    link_cost(17,[9,21,14]) =  link_cost(17,[9,21,14]) + river_plus;
    link_cost([9,21,14],17) =  link_cost([9,21,14],17) + river_plus;

    link_cost(10,[9,20,21]) = link_cost(10,[9,20,21]) + river_plus;
    link_cost([9,20,21],10) = link_cost([9,20,21],10) + river_plus;

    link_cost(8,[2,3,9]) = link_cost(8,[2,3,9]) + river_plus;
    link_cost([2,3,9],8) = link_cost([2,3,9],8) + river_plus;

    link_cost(1,[4,23,2,18,3,9]) = link_cost(1,[4,23,2,18,3,9]) + river_plus;
    link_cost([4,23,2,18,3,9],1) = link_cost([4,23,2,18,3,9],1) + river_plus;

    % link_cost(12,[9,14,20]) = link_cost(12,[9,14,20]) + river_plus;
    % link_cost([9,14,20],12) = link_cost([9,14,20],12) + river_plus;
    % 
    % link_cost(5,[2,9,21,14]) = link_cost(5,[2,9,21,14]) + river_plus;
    % link_cost([2,9,21,14],5) = link_cost([2,9,21,14],5) + river_plus;
    % 
    % link_cost(24,[2,23,3,15]) = link_cost(24,[2,23,3,15]) + river_plus;
    % link_cost([2,23,3,15],24) = link_cost([2,23,3,15],24) + river_plus;

    link_capacity_slope = 0.3.*link_cost; 

    station_cost = 1e3.*population./(365.25*25);
    

    station_capacity_slope = 0.2.*station_cost;
    
    
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
    sigma = 0.25;
    a_max = 1e9;
    eps = 1e-3;
end

function distancia = haversine(lat1, lon1, lat2, lon2)
    % Convierte las coordenadas de grados a radianes
    lat1 = deg2rad(lat1);
    lon1 = deg2rad(lon1);
    lat2 = deg2rad(lat2);
    lon2 = deg2rad(lon2);

    % Diferencias en coordenadas
    dlat = lat2 - lat1;
    dlon = lon2 - lon1;

    % Fórmula haversine
    a = sin(dlat/2)^2 + cos(lat1) * cos(lat2) * sin(dlon/2)^2;
    c = 2 * atan2(sqrt(a), sqrt(1-a));

    % Radio de la Tierra en kilómetros (aproximado)
    radio_tierra = 6371;

    % Calcula la distancia
    distancia = radio_tierra * c;
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



function [] = plot_network(A,lam,beta)

    coordinates = readtable('./coordenadas_Sevilla.xlsx','Sheet',1);
    id_node = table2array(coordinates(1:24,1));
    coor_x_raw = table2array(coordinates(1:24,7));
    coor_y_raw = table2array(coordinates(1:24,3));
    
    figure; hold off; axis equal;
    %geoscatter(lat, lon, 40, 'r', 'filled');  
    %hold on;
    figure('Position', [100, 100, 1000, 700]);
    geoscatter(coor_x_raw,coor_y_raw, 40, 'r', 'filled');  
    hold on;


    [i_idx, j_idx] = find(A > 1e-2);     % solo enlaces con valor positivo
    linecolor = [0, 0, 0.8, 0.3];           % negro semitransparente
    
    hold on;
    for k = 1:length(i_idx)
        i = i_idx(k);
        j = j_idx(k);
    
        val = A(i,j);                 % valor del enlace
        maxVal = max(A(:));
        lw = 0.75 + 3*(val/maxVal);

    
        geoplot([coor_x_raw(i), coor_x_raw(j)], [coor_y_raw(i), coor_y_raw(j)], ...
                'Color', linecolor, 'LineWidth', lw);
    end

    
    for n = 1:numel(coor_x_raw)
        text(coor_x_raw(n), coor_y_raw(n), sprintf('%d', n), ...
             'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'bottom', ...
             'FontSize', 11, ...
             'Color', 'black');
    end
    
    geobasemap('topographic');
    figurename = sprintf('./figures/topologia_lam=%d_beta=%d.png',round(lam),beta);
    saveas(gcf, figurename);
end

