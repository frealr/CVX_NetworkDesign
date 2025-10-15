


clear all; close all; clc;

[n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    op_link_cost,congestion_coef_stations,...
    congestion_coef_links,alt_price,a_nom,tau,sigma,...
    a_max,candidasourcertes] = parameters_large_sevilla_network();

station_cost = station_cost.*ones(1,n);
station_capacity_slope = station_capacity_slope.*ones(1,n);

disp(n);
%%

% === 1) RUTA de tu GAMS ===
gamsHome = 'C:\GAMS\50';  % <-- CAMBIA esto

[basedir,~,~] = fileparts(mfilename('fullpath'));
basedir = fullfile(basedir, 'export_csv');   % subcarpeta de exportación
if ~exist(basedir,'dir'); mkdir(basedir); end


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
n=51;


%betas for lam 5
lams = [5];

betas = 1e2:1e2:1e3; %por ahora estas para las nuevas simulaciones
betas = [550,725:25:900];
betas = [810,820,830,840];


lams = [10];
betas = 200:50:500;
betas = [455:5:550];
betas = [492];


% repito simulaciones 9/10

lams = [5];
betas = [1e-5,1e-4,1e-3];
betas = [1e-3:1e-3:1e-2];
betas = [0.0051:0.0001:0.0059,0.0062:0.0002:0.007];


lams = [10];
betas = [5e-4:5e-4:5e-3];

betas = 2.6e-3:1e-4:3.7e-3;



for ll = 1:length(lams)
   lam = lams(ll);
   for bb=1:length(betas)
       beta = betas(bb);
       tic;
       [sprim,s,deltas,aprim,a,deltaa,f,fext,fij,comp_time] = compute_sim(niters,lam,beta,alfa,n);

       [obj_val,pax_obj,op_obj] = get_obj_val(alfa, op_link_cost,...
            congestion_coef_links, ...
            congestion_coef_stations,prices,alt_price,aprim,deltaa, ...
            sprim,deltas,fij,f,fext,demand);
       budget = get_budget(s,sprim,a,aprim,n,...
            station_cost,station_capacity_slope,link_cost,link_capacity_slope,lam)
       filename = sprintf('./sevilla_rebutal/beta=%d_lam=%d.mat',beta,lam);
        save(filename,'s','sprim','deltas', ...
        'a','aprim','deltaa','f','fext','fij','comp_time','budget', ...
        'pax_obj','op_obj','obj_val');
   end

end






%% Load obtained results
lams = [10];
betas = [1e-3,5e-3,1e-2,5e-2,1e-1,3e-1,5e-1,7e-1,1];
betas = 1e-2:1e-2:7e-2;
betas = 4e-2:1e-3:6e-2;
betas = 0.046;
betas = [5e-2:1e-2:1e-1,6e-2:1e-3:8.2e-2];
betas = [1e-2,2.5e-2,5e-2,7.5e-2,1e-1:1e-2:3e-1,0.13:1e-3:0.15]; 
%betas = [1e-2:1e-2:7e-2,4e-2:1e-3:6e-2]; %lam 15
%betas = [0,5e-2:1e-2:1e-1,6e-2:1e-3:8.2e-2]; %lam 10
betas = [1e-2,2.5e-2,5e-2,7.5e-2,1e-1:1e-2:3e-1,0.13:1e-3:0.15, 1e-1:2.5e-3:1.3e-1]; %lam 5
%betas = 1e-2:2.5e-4:1.3e-2;

% betas = [1e-3,1e-2,1e-1,1,1e1,1e2,1e3];
% betas = [5e2,1e3,1e4];
% betas = [6e2:1e2:1.5e3];
% 
% betas = [1.1e3:25:1.275e3];
% betas = [6e2:1e2:1.5e3];
% 
% betas = [1e2:1e2:7e2];
% 
% betas = [525,550,575,1e2:1e2:7e2];

betas = [1e-3,1e-2,1e-1,1,1e1,1e2,1e3,6e2:1e2:1.5e3,1.1e3:25:1.275e3]; %betas for lam 5 51 nodos
betas = [600,800,1000,1150,1175,1200,1250,1275]; %betas representar para lam 5

%betas = [6e2:1e2:1.5e3,1e2:1e2:5e2,525,550,575]; % betas for lam 10 51 nodos
betas = [300,400,500,575,600]; %representar para lam10, faltan algunas mas

%betas = [6e2:1e2:1.5e3];
%betas = [1.1e3:25:1.275e3];

betas = [1e-3,1e-2,1e-1,1,1e1,1e2,1e3];
betas = 1e2:1e2:1e3;
betas = [550,725:25:900];
betas = [800,810,820,830,840,850,875];
lams = [5];

betas = [400,500,550,600,700,725,750, 775, 800, 810, 820, 825, 830, 840, 850, 875, 900]; % seleccionar de aqui
betas = [400,600,750,800,810,820,840,850,875,900]; %seleccion para tabla lam 5

% 
% % betas = 200:50:500;
% lams = [10];
% 
% betas = [200:50:450,455:5:550]; %seleccionar de aqui
% betas = [250,350,450,475,485,492,500,510,520,525]; %seleccion para tabla lam 10

betas = [1e-2,1e-1,1,1e1,1e2,1e3];
betas = 1e-2;
betas = [1e-5,1e-4,1e-3];
betas = [1e-3:1e-3:1e-2];
betas = [0.0051:0.0001:0.0059,0.0062:0.0002:0.007];

lams = [5];
betas = [1e-5,1e-4,1e-3];
betas = [1e-3:1e-3:1e-2];
betas = [0.0051:0.0001:0.0059,0.0062:0.0002:0.007];
betas = [1e-3,4e-3,5.2e-3,5.5e-3,5.7e-3,6e-3,6.2e-3,6.6e-3,6.8e-3,7e-3];

lams = [10];
betas = [5e-4:5e-4:5e-3];
betas = 2.6e-3:1e-4:3.7e-3;
betas = [1e-3,2e-3,2.7e-3,2.9e-3,3.2e-3,3.3e-3,3.4e-3,3.5e-3,3.6e-3,3.7e-3];

for ll=1:length(lams)
    lam = lams(ll);
    for bb=1:length(betas)
        beta_or = betas(bb);
        [n,link_cost,station_cost,link_capacity_slope,...
        station_capacity_slope,demand,prices,...
        op_link_cost,congestion_coef_stations,...
        congestion_coef_links,alt_price,a_nom,tau,sigma,...
        a_max,candidasourcertes] = parameters_large_sevilla_network();
      
        filename = sprintf('./sevilla_rebutal/beta=%d_lam=%d.mat',beta_or,lam);
        load(filename);
        att_d(bb) = 100*sum(sum(f.*demand))/sum(sum(demand));
        bud(bb) = budget;
        att_dem = round(100*sum(sum(f.*demand))/sum(sum(demand)),1);
        n_links = sum(sum(aprim));
        n_nodes = sum(sprim > 0.1);
        served_demand = round(100*sum(sum(f.*demand))./sum(sum(demand)),1);
        %  disp(['beta = ',num2str(beta_or), ', lam = ', num2str(lam)]);
        %  disp(['att_dem = ',num2str(served_demand),' %']);
        % % % disp(['Pres. = ',num2str(budget)]);
        %  disp(['Arcos = ',num2str(sum(sum(aprim > 0.1)))]);
        % disp(['Nodos = ',num2str(sum(sprim > 0.1))]);
        % disp(['Cap. = ',num2str(sum(sum(aprim)))]);
        % disp(['distancia por pasajero red nueva = ',num2str(d_med(bb))]);
        % disp(['distancia por pasajero red actual = ',num2str(u_med(bb))]);
        % disp(['distancia por ruta = ',num2str(long_mean(bb))]);
        % disp(['tiempo computacional = ',num2str(comp_time)]);
        % disp('\n');
        disp([ num2str(lam),'&',...
            sprintf('%.2e',beta_or),'&', ...
            num2str(served_demand), ...
            '&',sprintf('%.3e',budget),'&',num2str(sum(sum(a > 0.1))), ...
            '&',num2str(sum(s > 0.1)), ...
            '&',sprintf('%.2f',n_links), ...
            '&',sprintf('%.2e',comp_time),'\\ \hline']);
    
    
        
    end
end

%% Plot
close all;
betas = [1e-3,4e-3,5.2e-3,6e-3]; %lam 5
betas = [1e-3,2e-3,2.7e-3,3.2e-3]; %lam 10
lam = 10;


    

for bb=1:length(betas)
    beta = betas(bb);
    filename = sprintf('./sevilla_rebutal/beta=%d_lam=%d.mat',beta,lam);
    load(filename);

    plot_network(a,'g',lam,beta);
    figurename = sprintf('./figures/topologia_lam=%d_beta=%d.png',round(lam),beta);
    saveas(gcf, figurename);
end

%%
lam = 0;
beta = 0;
plot_network(zeros(51),'g',lam,beta);
figurename = sprintf('./figures/topologia_51_example.png');
saveas(gcf, figurename);
%%
    close all;
    n_fs = 140; %number of demand centroids
    
    nlines = 8; %resolution for the grid
    sigmaX = 3;
    sigmaY= sigmaX;

    plot_network(a,'g');

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
    congestion_coef_stations,prices,alt_price,a_prim,delta_a, ...
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
                    pax_obj = pax_obj + 1e-3*(demand(o,d).*prices(i,j)).*fij(i,j,o,d);
                end
            end
        end
    end
    pax_obj = pax_obj + 1e-3*(sum(sum(demand.*(alt_price).*fext)));

    entro = max(f.*(log(f+eps)-1),0);
    entro_ext = max(fext.*(log(fext+eps)-1),0);
    pax_obj = pax_obj + 1e-3*(sum(sum(demand.*(entro))));
    pax_obj = pax_obj + 1e-3*(sum(sum(demand.*(entro_ext))));
    obj_val = (alfa/(1))*pax_obj + ((1-alfa)/(1))*op_obj;
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
    
    gmsFile  = 'C:\Users\freal\MATLAB\Projects\untitled\code\sevilla_net_49\new1.gms';
    comp_time = 0;

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

        results_file_ctime = readtable('./output_all.xlsx','Sheet','solver_time');
        comp_time = comp_time + table2array(results_file_ctime);
    
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
    congestion_coef_links,alt_price,a_nom,tau,sigma,...
    a_max,candidates] = parameters_large_sevilla_network()

    nlines = 8;
    n_fs = 140; %number of demand centroids
    
    sigmaX = 3;
    sigmaY= sigmaX;

    [n,demand,candidates,distance,lat,xpos,distance_alt] = set_network(n_fs,...
    nlines,sigmaX,sigmaY,'g');
    
    
    %fixed cost for constructing links

    link_cost = 2e7.*distance./(365.25*25);

    link_capacity_slope = 0.8.*link_cost; 

    station_cost = 1e7./(365.25*25);

    station_capacity_slope = 0.6.*station_cost;
    
    
    % Op Link Cost
    op_link_cost = 3.*distance;
    
    % Congestion Coefficients
    congestion_coef_stations = 0.1 .* ones(n, 1);
    congestion_coef_links = 0.1 .* ones(n);
    
    % Prices
    riding_cost = 0.083; %e/min
    train_speed = 40; %km/h
    fare = 0; %e

    prices = fare + (riding_cost*60/train_speed).*distance;

    %Alternative: car
    
    fixed_alt_price = 1.75; %e
    variable_alt_price = 0.12; %e/km
    alt_time_val = 0.05; %e/min
    
    average_parking_time = 10; %min
    parking_cost = 0.25; %e/min
    fixed_parking_cost = 1; %e
    alt_speed = 60; %km/h
    
    alt_price = variable_alt_price.*distance_alt ...
    + alt_time_val*(60/alt_speed).*distance_alt;
    
    a_nom = 607;             
    
    tau = 0.57;
    sigma = 0.25;
    a_max = 1e9;
    eps = 1e-3;
end

function [] = plot_network(a,cell_distribution,lam,beta)

    n_fs = 140; %number of demand centroids
    nlines = 8; %resolution for the grid
    sigmaX = 3;
    sigmaY= sigmaX;

    coordinates = readtable('./NODOSSEV_MT.xlsx','Sheet',1);
    id_node = table2array(coordinates(1:140,1));
    coor_x_raw = table2array(coordinates(1:140,2));
    coor_y_raw = table2array(coordinates(1:140,3));
    
    %real distance between nodes 1 and 60
    lat1 = 37.413969; lon1= -6.007699;
    lat60 = 37.361042; lon60 = -5.959065;
    real_d = haversine(lat1,lon1,lat60,lon60);
    
    %scaled distance between nodes 1 and 60
    scaled_d = sqrt( (coor_x_raw(id_node == 1) - coor_x_raw(id_node == 60))^2 ...
        +(coor_y_raw(id_node == 1) - coor_y_raw(id_node == 60) )^2   );
    
    scale = real_d / scaled_d; % scale of coordinates
    
    coor_x(id_node) = coor_x_raw*scale;
    coor_y(id_node) = coor_y_raw*scale;
    
    % constant values
    
    meters_per_deg_lat = 111.320; 
    meters_per_deg_lon = 111.320 * cosd(lat1);  
    
    % distance wrt reference node
    dx = coor_x - coor_x(id_node==1);
    dy = coor_y - coor_y(id_node==1);
    
    % Approximate lat/lon
    lat = lat1 + dy / meters_per_deg_lat;
    lon = lon1 + dx / meters_per_deg_lon;
    lon(18) = lon(18) + 0.0005;
    lat(18) = lat(18) + 0.0005;
    
    lon(34) = lon(34) + 0.0005;
    lat(34) = lat(34) + 0.0005;
    %
    
    demand_raw = readtable('./NODOSSEV_MT.xlsx','Sheet',2);
    id_node_x = table2array(demand_raw(1,2:141));
    id_node_y = table2array(demand_raw(2:141,1));
    
    demand_aux = table2array(demand_raw(2:141,2:141));
    
    demand(id_node_x,id_node_y) = demand_aux;
    
    %full adyacency matrix
    A_full = zeros(n_fs);
    G_full = graph(A_full);
    figure;
    
    scatter(coor_x,coor_y,'filled');
    hold on;
    xlims = xlim; ylims = ylim;
    
    ix = (1:nlines) - 0.5 - nlines/2;
    iy = (1:nlines) - 0.5 - nlines/2;
    
    densX = exp(-(ix.^2)/(2*sigmaX^2));
    densY = exp(-(iy.^2)/(2*sigmaY^2));
    
    weights = @(n,sigma)...
        exp(-(( (1:n) -0.5 - n/2  ).^2 ) / (2*sigma^2));
    normalize = @(v) v / sum(v);

    switch cell_distribution
        case 'u'
            % uniform split for cells
            wx = ones(1,nlines)/nlines;
            wy = ones(1,nlines)/nlines;
        otherwise
            % gaussian split for cells
            wx = normalize(1./densX);
            wy = normalize(1./densY);
    end
    
    widths = wx * (xlims(2) - xlims(1));
    heights = wy * (ylims(2) - ylims(1));
    
    xedges = [xlims(1),xlims(1) + cumsum(widths)];
    
    yedges =  [ylims(1),ylims(1) + cumsum(heights)];
    
    x_cent = (xedges(1:end-1)+xedges(2:end))/2;
    y_cent = (yedges(1:end-1)+yedges(2:end))/2;
    
    [Xc,Yc] = meshgrid(x_cent,y_cent);
    
    xc = Xc(:);
    yc = Yc(:);
    
    for xv=xedges
        xline(xv,':k','LineWidth',0.8);
    end
    for yv=yedges
        yline(yv,':k','LineWidth',0.8);
    end
    
    
    N = numel(coor_x);
    
    ix = discretize(coor_x, xedges);  
    iy = discretize(coor_y, yedges); 
    
    valid = ~isnan(ix) & ~isnan(iy);     
    cell_lin = nan(N,1);
    cell_lin(valid) = sub2ind([nlines nlines], iy(valid), ix(valid));   
    
    K = nlines*nlines;  
    
    
    demanda_total = sum(demand, 2) + sum(demand, 1)';
    
    centroides_x = nan(K,1);
    centroides_y = nan(K,1);
    
    for k = 1:K
        % nodes in cell k
        idx_nodos = find(cell_lin == k);
    
        if isempty(idx_nodos)
            centroides_x(k) = Xc(k);
            centroides_y(k) = Yc(k);
            continue;  % empty cell
        end
    
        % total demand per cell
        dt = demanda_total(idx_nodos);
    
        % largest demand node index per cell, allocate station there
        [~, idx_max] = max(dt);
        idx_nodo_max = idx_nodos(idx_max);
    
        centroides_x(k) = lat(idx_nodo_max);
        centroides_y(k) = lon(idx_nodo_max);
    end
    
    [I,J,V] = find(demand);                   
    ci = cell_lin(I);  cj = cell_lin(J);  
    
    mask = ~isnan(ci) & ~isnan(cj);      
    demand_cells = accumarray([ci(mask), cj(mask)], V(mask), [K K], @sum, 0);
    demand_cells(1:size(demand_cells,1)+1:end) = 0;
    
    
    nodos_inactivos = 0;
    for i=1:length(demand_cells)
        if sum(demand_cells(i,:)) == 0
            nodos_inactivos = nodos_inactivos + 1;
        end
    end
    
    n = size(Xc, 1);            
    K = n * n;
    
    xpos = centroides_x;
    ypos = centroides_y;   

    
    % aggregate demand per cell
    demanda_por_celda = accumarray(cell_lin(valid), demanda_total(valid), [K 1], @sum, 0);
    dtotal = demand_cells + demand_cells';
    zonas_con_demanda = find(sum(dtotal) > 0);
    
    % potential link matrix
    A_sub = a;
    demand_cells_sub = demand_cells(zonas_con_demanda,zonas_con_demanda);
    
    % coordinates of cells with demand
    x_sub = xpos(zonas_con_demanda);
    y_sub = ypos(zonas_con_demanda);
    
    figure; hold off; axis equal;
    %geoscatter(lat, lon, 40, 'r', 'filled');  
    %hold on;
    figure('Position', [100, 100, 1000, 700]);
    geoscatter(xpos(xpos > 30), ypos(ypos < -3), 40, 'r', 'filled');  
    hold on;



    [i_idx, j_idx] = find(A_sub > 1e-4);     % solo enlaces con valor positivo
    linecolor = [0, 0, 0.8, 0.3];           % negro semitransparente
    
    hold on;
    for k = 1:length(i_idx)
        i = i_idx(k);
        j = j_idx(k);
    
        val = A_sub(i,j);                 % valor del enlace
        maxVal = max(A_sub(:));
        lw = 0.75 + 3*(val/maxVal);

    
        geoplot([x_sub(i), x_sub(j)], [y_sub(i), y_sub(j)], ...
                'Color', linecolor, 'LineWidth', lw);
    end

    
    for n = 1:numel(x_sub)
        text(x_sub(n), y_sub(n), sprintf('%d', n), ...
             'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'bottom', ...
             'FontSize', 11, ...
             'Color', 'black');
    end
    
    geobasemap('topographic');
    figurename = sprintf('./figures/topologia_lam=%d_beta=%d.png',round(lam),beta);
    saveas(gcf, figurename);
end


function [n,demand_cells_sub,candidates,distance,lat,xpos,distance_alt] = set_network(n_fs,...
    nlines,sigmaX,sigmaY,cell_distribution)

    coordinates = readtable('./NODOSSEV_MT.xlsx','Sheet',1);
    id_node = table2array(coordinates(1:140,1));
    coor_x_raw = table2array(coordinates(1:140,2));
    coor_y_raw = table2array(coordinates(1:140,3));
    
    %real distance between nodes 1 and 60
    lat1 = 37.413969; lon1= -6.007699;
    lat60 = 37.361042; lon60 = -5.959065;
    real_d = haversine(lat1,lon1,lat60,lon60);
    
    %scaled distance between nodes 1 and 60
    scaled_d = sqrt( (coor_x_raw(id_node == 1) - coor_x_raw(id_node == 60))^2 ...
        +(coor_y_raw(id_node == 1) - coor_y_raw(id_node == 60) )^2   );
    
    scale = real_d / scaled_d; % scale of coordinates
    
    coor_x(id_node) = coor_x_raw*scale;
    coor_y(id_node) = coor_y_raw*scale;
    
    % constant values
    
    meters_per_deg_lat = 111.320; 
    meters_per_deg_lon = 111.320 * cosd(lat1);  
    
    % distance wrt reference node
    dx = coor_x - coor_x(id_node==1);
    dy = coor_y - coor_y(id_node==1);
    
    % Approximate lat/lon
    lat = lat1 + dy / meters_per_deg_lat;
    lon = lon1 + dx / meters_per_deg_lon;
    lon(18) = lon(18) + 0.0005;
    lat(18) = lat(18) + 0.0005;
    
    lon(34) = lon(34) + 0.0005;
    lat(34) = lat(34) + 0.0005;
    %
    
    demand_raw = readtable('./NODOSSEV_MT.xlsx','Sheet',2);
    id_node_x = table2array(demand_raw(1,2:141));
    id_node_y = table2array(demand_raw(2:141,1));
    
    demand_aux = table2array(demand_raw(2:141,2:141));
    
    demand(id_node_x,id_node_y) = demand_aux;
    
    %full adyacency matrix
    A_full = zeros(n_fs);
    G_full = graph(A_full);
    figure;
    
    scatter(coor_x,coor_y,'filled');
    hold on;
    xlims = xlim; ylims = ylim;
    
    ix = (1:nlines) - 0.5 - nlines/2;
    iy = (1:nlines) - 0.5 - nlines/2;
    
    densX = exp(-(ix.^2)/(2*sigmaX^2));
    densY = exp(-(iy.^2)/(2*sigmaY^2));
    
    weights = @(n,sigma)...
        exp(-(( (1:n) -0.5 - n/2  ).^2 ) / (2*sigma^2));
    normalize = @(v) v / sum(v);

    switch cell_distribution
        case 'u'
            % uniform split for cells
            wx = ones(1,nlines)/nlines;
            wy = ones(1,nlines)/nlines;
        otherwise
            % gaussian split for cells
            wx = normalize(1./densX);
            wy = normalize(1./densY);
    end
    
    widths = wx * (xlims(2) - xlims(1));
    heights = wy * (ylims(2) - ylims(1));
    
    xedges = [xlims(1),xlims(1) + cumsum(widths)];
    
    yedges =  [ylims(1),ylims(1) + cumsum(heights)];
    
    x_cent = (xedges(1:end-1)+xedges(2:end))/2;
    y_cent = (yedges(1:end-1)+yedges(2:end))/2;
    
    [Xc,Yc] = meshgrid(x_cent,y_cent);
    
    xc = Xc(:);
    yc = Yc(:);
    
    for xv=xedges
        xline(xv,':k','LineWidth',0.8);
    end
    for yv=yedges
        yline(yv,':k','LineWidth',0.8);
    end
    
    
    N = numel(coor_x);
    
    ix = discretize(coor_x, xedges);  
    iy = discretize(coor_y, yedges); 
    
    valid = ~isnan(ix) & ~isnan(iy);     
    cell_lin = nan(N,1);
    cell_lin(valid) = sub2ind([nlines nlines], iy(valid), ix(valid));   
    
    K = nlines*nlines;  
    
    
    demanda_total = sum(demand, 2) + sum(demand, 1)';
    
    centroides_x = nan(K,1);
    centroides_y = nan(K,1);
    
    for k = 1:K
        % nodes in cell k
        idx_nodos = find(cell_lin == k);
    
        if isempty(idx_nodos)
            centroides_x(k) = Xc(k);
            centroides_y(k) = Yc(k);
            continue;  % empty cell
        end
    
        % total demand per cell
        dt = demanda_total(idx_nodos);
    
        % largest demand node index per cell, allocate station there
        [~, idx_max] = max(dt);
        idx_nodo_max = idx_nodos(idx_max);
    
        centroides_x(k) = lat(idx_nodo_max);
        centroides_y(k) = lon(idx_nodo_max);
    end
    
    [I,J,V] = find(demand);                   
    ci = cell_lin(I);  cj = cell_lin(J);  
    
    mask = ~isnan(ci) & ~isnan(cj);      
    demand_cells = accumarray([ci(mask), cj(mask)], V(mask), [K K], @sum, 0);
    demand_cells(1:size(demand_cells,1)+1:end) = 0;
    
    
    nodos_inactivos = 0;
    for i=1:length(demand_cells)
        if sum(demand_cells(i,:)) == 0
            nodos_inactivos = nodos_inactivos + 1;
        end
    end
    
    n = size(Xc, 1);            
    K = n * n;
    
    % Initialize candidate links matrix
    A = zeros(K, K);
    
    % Connection probability
    p = 0.8;
    rng(123);
    
    for row = 1:n
        for col = 1:n
            i = sub2ind([n n], row, col); 
    
            neighbors = [];
    
            if row > 1     
                j = sub2ind([n n], row-1, col);
                if j > i && rand < p
                    A(i, j) = 1;
                    A(j, i) = 1; 
                end
            end
            if row < n    
                j = sub2ind([n n], row+1, col);
                if j > i && rand < p
                    A(i, j) = 1;
                    A(j, i) = 1;
                end
            end
            if col > 1  
                j = sub2ind([n n], row, col-1);
                if j > i && rand < p
                    A(i, j) = 1;
                    A(j, i) = 1;
                end
            end
            if col < n 
                j = sub2ind([n n], row, col+1);
                if j > i && rand < p
                    A(i, j) = 1;
                    A(j, i) = 1;
                end
            end
        end
    end
    
    xpos = centroides_x;
    ypos = centroides_y;   
    G = graph(A);
    
    % aggregate demand per cell
    demanda_por_celda = accumarray(cell_lin(valid), demanda_total(valid), [K 1], @sum, 0);
    dtotal = demand_cells + demand_cells';
    zonas_con_demanda = find(sum(dtotal) > 0);
    G_sub = subgraph(G, zonas_con_demanda);
    
    % potential link matrix
    A_sub = full(adjacency(G_sub));
    demand_cells_sub = demand_cells(zonas_con_demanda,zonas_con_demanda);
    
    % coordinates of cells with demand
    x_sub = xpos(zonas_con_demanda);
    y_sub = ypos(zonas_con_demanda);
    
    figure; hold off; axis equal;
    %geoscatter(lat, lon, 40, 'r', 'filled');  
    %hold on;
    geoscatter(xpos(xpos > 30), ypos(ypos < -3), 40, 'r', 'filled');  
    hold on;
    
    
    [i_idx, j_idx] = find(A_sub); 
    linecolor = [0, 0, 0, 0.2];
    for k = 1:length(i_idx)
        i = i_idx(k);
        j = j_idx(k);
        geoplot([x_sub(i), x_sub(j)], [y_sub(i), y_sub(j)], 'k-', 'LineWidth', 1.5,'Color',linecolor); hold on;
    end
    
    for n = 1:numel(x_sub)
        text(x_sub(n), y_sub(n), sprintf('%d', n), ...
             'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'bottom', ...
             'FontSize', 11, ...
             'Color', 'black');
    end
    
    geobasemap('topographic');
    n = length(x_sub);
    
    coor_x_sub = centroides_x(zonas_con_demanda);
    coor_y_sub = centroides_y(zonas_con_demanda);
    
    K = size(A_sub, 1);
    candidates = cell(K, 1);
    
    for i = 1:K
        candidates{i} = find(A_sub(i, :) == 1);
    end
    distance = 1e6.*ones(n);
    for i=1:n
        distance(i,i) = 0;
        cand = candidates(i);
        cand = cand{1};
        cand = cand(cand > i);
        for j=i+1:n
            if sum(j == cand) > 0
                distance(i,j) = haversine(coor_x_sub(i), coor_y_sub(i), coor_x_sub(j), coor_y_sub(j));
                distance(j,i) = distance(i,j);
            end
        end
    end
    distance_alt = 1e6.*ones(n);
    for i=1:n
        distance(i,i) = 0;
        for j=i+1:n
            distance_alt(i,j) = haversine(coor_x_sub(i), coor_y_sub(i), coor_x_sub(j), coor_y_sub(j));
            distance_alt(j,i) = distance_alt(i,j);
        end
    end

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
