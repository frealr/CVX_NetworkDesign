close all;
clear all;
clc;
cvx_begin
cvx_end;

%Problema del operador:

%op_obj = operation_costs. Problema del operador.

%Problema de los pasajeros:
%pax_obj = (distances + prices) + entropies + congestion


cvx_solver_settings -clearall
cvx_solver mosek
cvx_precision high
cvx_save_prefs 


list_of_parameters = {'demand','alt_time','prices','alt_price','fixed_costs',...
    'linear_costs','op_cost'};

eps = 1e-3;
alfa = 0.6;
dm_pax = 1.2;
dm_op = 0.008;

betas = [1.25,1.5];
lams = [5,6];

for lind=1:length(list_of_parameters)
    for ll=1:length(lams)
        for bb=1:length(betas)
            lam = lams(ll);
            beta = betas(bb);
            param_name = string(list_of_parameters(lind));
            if ~((beta == 1.5) & (lam == 6))
                disp(['simulation for parameter ',param_name,', beta = ',...
                    num2str(beta),', lambda = ',num2str(lam)]);
                filename = sprintf('./uthopy_results/sol_beta=%d_lam=%d.mat',beta,lam);
                load(filename);
                nom_bud = budget;
                
                sens_coefs = 0.8:0.05:1.2;
                beta_or = beta;
                
                for ss=1:length(sens_coefs)
                    sens_coef = sens_coefs(ss);
                    conv = 0;
                    beta = beta_or;
                    beta_max = 5;
                    beta_min = 0;
                    counter = 0;
                    disp(['coef = ',num2str(sens_coef)]);
                    
                    while ~conv
                        [s,s_prim,delta_s,a,a_prim,delta_a,f,fext,fij,comp_time,budget,obj_val,...
                        pax_obj,op_obj] = compute_sim_sens_coef(lam,beta,alfa,dm_pax,dm_op,sens_coef,param_name);
                        if ((abs((budget - nom_bud))/(nom_bud) <= 0.01) | (counter > 15)) & (budget >= nom_bud) %si no converge lo bajas. %si se va mucho, haz 2 scripts, uno que se quede
                            %siempre por debajo del presupuesto y otro por encima, y luego
                            %haces la combinacion lineal.
                            conv = 1;
                            disp('converged');
                        end
                        if (budget < nom_bud)
                            beta_max = beta;
                        elseif (budget > nom_bud)
                            beta_min = beta;
                        end
                        beta = (beta_max+beta_min)/2;
                        disp(['budget = ',num2str(budget)]);
                        disp(['counter = ', num2str(counter)]);
                        counter = counter + 1;
                    
                    end
                    filename = sprintf('./uthopy_results/sol_beta=%d_lam=%d_sens_analysis_%s_coef=%d_moreb.mat',beta_or,lam,param_name,sens_coef);
                    save(filename,'s','s_prim','delta_s', ...
                        'a','a_prim','delta_a','f','fext','fij','obj_val','pax_obj','op_obj','comp_time','budget');
                end
            end
        end
    end
end


%% 

% 
% eps = 1e-3;
% alfa = 0.6;
% lam = 5;
% beta = 1.25;
% dm_pax = 1.2;
% dm_op = 0.008;
% 
% sens_coefs = 0.8:0.05:1.2;
% beta_or = beta;
% 
% for ss=1:length(sens_coefs)
%     sens_coef = sens_coefs(ss);
%     conv = 0;
%     beta = beta_or;
%     disp(['coef = ',num2str(sens_coef)]);
%     filename = sprintf('./uthopy_results/sol_beta=%d_lam=%d_sens_analysis_op_cost_coef=%d_lessb.mat',beta_or,lam,sens_coef);
%     load(filename);
%     disp(['obj_val = ',num2str(obj_val),', pax_obj = ',num2str(pax_obj), ...
%         ', op_obj = ',num2str(op_obj), ', budget = ',num2str(budget)]);
% 
%     
% end
% 


%% Functions

function budget = get_budget(s,s_prim,a,a_prim,n,...
    station_cost,station_capacity_slope,link_cost,link_capacity_slope,lam)
    budget = 0;
    for i=1:n
        if s(i) > 1e-3
            budget = budget + lam*station_cost(i) + ...
                station_capacity_slope(i)*s_prim(i);
        end
        for j=1:n
            if a(i,j) > 1e-3
                budget = budget + lam*link_cost(i,j) + ...
                    link_capacity_slope(i,j) * a_prim(i,j);
            end
        end
    end
end

function [obj_val,pax_obj,op_obj] = get_obj_val(alfa, op_link_cost,...
    congestion_coef_links, ...
    congestion_coef_stations,travel_time,prices,alt_time,alt_price,a_prim,delta_a, ...
    s_prim,delta_s,fij,f,fext,demand,dm_pax,dm_op)
    n = 9;
    pax_obj = 0;
    op_obj = 0;
    eps = 1e-3;
    op_obj = op_obj + 1e-6*(sum(sum(op_link_cost.*a_prim))); %operational costs
    for i=1:n
        if s_prim(i) > eps
            pax_obj = pax_obj + 1e-6*inv_pos(congestion_coef_stations(i)*delta_s(i) + eps);
        end

        for j=1:n
            if a_prim(i,j) > eps
                pax_obj = pax_obj + 1e-6*inv_pos(congestion_coef_links(i,j)*delta_a(i,j) + eps);
            end
        end
    end
    for o=1:n
        for d=1:n
            for i=1:n
                for j=1:n
                    pax_obj = pax_obj + 1e-6*(demand(o,d).*(travel_time(i,j)+prices(i,j)).*fij(i,j,o,d));
                end
            end
        end
    end
    pax_obj = pax_obj + 1e-6*(sum(sum(demand.*(alt_time+alt_price).*fext)));
    pax_obj = pax_obj + 1e-6*(sum(sum(demand.*(-entr(f) - f))));
    pax_obj = pax_obj + 1e-6*(sum(sum(demand.*(-entr(fext) - fext))));
    obj_val = pax_obj*alfa/dm_pax + op_obj*(1-alfa)/dm_op;
end


function [s,s_prim,delta_s,a,a_prim,delta_a,f,fext,fij,comp_time,budget,obj_val,...
    pax_obj,op_obj] = compute_sim(lam,beta_or,alfa,dm_pax,dm_op)

    [n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    load_factor,op_link_cost,congestion_coef_stations,...
    congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,eta,...
    a_max,candidates] = parameters_9node_network();

    niters = 15;           
    eps = 1e-3;
    a_prev = 1e4.*ones(n);
    s_prev= 1e4.*ones(n,1);
    disp(['beta = ',num2str(beta_or),', lam = ',num2str(lam)]);
    tic;
    for iter=1:niters
        cvx_begin quiet
            variable s(n)
            variable s_prim(n)
            variable delta_s(n)
            variable a(n,n)
            variable a_prim(n,n)
            variable delta_a(n,n)
            variable f(n,n)
            variable fext(n,n)
            variable fij(n,n,n,n)
            op_obj = 0;
            pax_obj = 0;
            bud_obj = 0;
            bud_obj = bud_obj + 1e-6*sum(station_capacity_slope'.*s_prim);
            bud_obj = bud_obj + 1e-6*(sum(sum(link_capacity_slope.*a_prim)));  %linear construction costs
            op_obj = op_obj + 1e-6*(sum(sum(op_link_cost.*a_prim))); %operation costs
            if iter < niters
                pax_obj = pax_obj + 1e-6*(sum(sum(inv_pos(congestion_coef_links.*delta_a + eps)))) + 1e-6*(sum(inv_pos(congestion_coef_stations'.*delta_s + eps))); %congestion costs
                bud_obj = bud_obj + 1e-6*lam*sum(sum((link_cost.*a_prim.*(1./(a_prev+eps))))) + 1e-6*lam*sum((station_cost'.*s_prim.*(1./(s_prev+eps)))); %fixed construction costs
                
            end
    
            for o=1:n
                for d=1:n
                    pax_obj = pax_obj + 1e-6*(demand(o,d).*sum(sum((travel_time+prices).*fij(:,:,o,d)))); 
                end
            end
            pax_obj = pax_obj + 1e-6*(sum(sum(demand.*(alt_time+alt_price).*fext)));
            pax_obj = pax_obj + 1e-6*(sum(sum(demand.*(-entr(f) - f))));
            pax_obj = pax_obj + 1e-6*(sum(sum(demand.*(-entr(fext) - fext))));
    
    
            if iter == niters
                for i=1:n
                    if s_prev(i) >= eps
                        pax_obj = pax_obj + 1e-6*(sum(inv_pos(congestion_coef_stations(i).*delta_s(i) + eps)));
                    end
                    for j=1:n
                        if a_prev(i,j) >= eps
                            pax_obj = pax_obj + 1e-6*(sum(sum(inv_pos(congestion_coef_links(i,j).*delta_a(i,j) + eps))));
                         end
                     end
                 end
            end
            obj = beta_or*bud_obj + (alfa/(dm_pax))*pax_obj + ((1-alfa)/(dm_op))*op_obj;
            minimize obj
            % constraints
            s >= 0;
            s_prim >= 0;
            delta_s >= 0;
            a >= 0;
            a_prim >= 0;
            delta_a >= 0;
            f >= 0;
            f <= 1;
            fij >= 0;
            fij <= 1;
            fext >= 0;
            fext <= 1;
            s_prim == s + delta_s;
            a_prim == a + delta_a;
    
            for i=1:n
                for j=1:n
                    squeeze(sum(sum(squeeze(permute(fij(i,j,:,:),[3 4 1 2]).*demand)))) <= tau.*a(i,j).*a_nom; %multiplicar por demanda
                end
            end
            for i=1:n
                eta*sum(a(:,i)) <= load_factor(i).*s(i)
            end
            sum(sum(a)) <= a_max;
            for o=1:n
                for d=1:n
                    sum(fij(o,:,o,d)) == f(o,d);
                end
            end
            
            for o=1:n
                squeeze(sum(fij(o,:,o,[1:(o-1),(o+1):n]),2)) - squeeze(sum(permute(fij(:,o,o,[1:(o-1),(o+1):n]),[2,1,3,4]),2)) == transpose(1 - fext(o,[1:(o-1),(o+1):n])); 
            end
            for d=1:n
                squeeze(sum(fij(d,:,[1:(d-1),(d+1):n],d),2)) - squeeze(sum(permute(fij(:,d,[1:(d-1),(d+1):n],d),[2,1,3,4]),2)) == -1 + fext([1:(d-1),(d+1):n],d);
            end
            for i=1:n
                fij(i,i,:,:) == 0;
                squeeze(sum(fij(i,:,[1:(i-1),(i+1):n],[1:(i-1),(i+1):n]),2)) - squeeze(sum(permute(fij(:,i,[1:(i-1),(i+1):n],[1:(i-1),(i+1):n]),[2,1,3,4]),2)) == 0;
            end
            for o=1:n
                fext(o,o) == 0;
                f(o,o) == 0;
                fij(:,o,o,:) == 0;
            end
            for d=1:n
                fij(d,:,:,d) == 0;
            end
            for o=1:n
                for d=1:n
                    if o ~= d
                        f(o,d) + fext(o,d) == 1;
                    end
                end
            end
            for i=1:n
                for j=1:n
                    if ~ismember(j,candidates{i})
                        a_prim(i,j) == 0;
                    end
                end
            end
            if iter == niters
                for i=1:n
                    if s_prev(i) <= 0.1
                        s_prim(i) == 0;
                    end
                    for j=1:n
                        if a_prev(i,j) <= 0.1
                            a_prim(i,j) == 0;
                         end
                     end
                 end
             end
    
        cvx_end
        a_prev = a_prim;
        s_prev = s_prim;
    end
    comp_time = toc;
    
    
    budget = get_budget(s,s_prim,a,a_prim,n,...
        station_cost,station_capacity_slope,link_cost,link_capacity_slope,lam);
    
    [obj_val,pax_obj,op_obj] = get_obj_val(alfa, op_link_cost,...
        congestion_coef_links, ...
        congestion_coef_stations,travel_time,prices,alt_time,alt_price,a_prim,delta_a, ...
        s_prim,delta_s,fij,f,fext,demand,dm_pax,dm_op);
    filename = sprintf('./uthopy_results/sol_beta=%d_lam=%d_uthopy_alfa=%d.mat',beta_or,lam,alfa);
    %save(filename,'s','s_prim','delta_s', ...
     %   'a','a_prim','delta_a','f','fext','fij','obj_val','pax_obj','op_obj','comp_time','budget');

end



function [s,s_prim,delta_s,a,a_prim,delta_a,f,fext,fij,comp_time,budget,obj_val,...
    pax_obj,op_obj] = compute_sim_sens_coef(lam,beta_or,alfa,dm_pax,dm_op,sens_coef,parameter_name)

    [n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    load_factor,op_link_cost,congestion_coef_stations,...
    congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,eta,...
    a_max,candidates] = parameters_9node_network();

    if parameter_name == 'alt_time'
        alt_time = sens_coef.*alt_time;
    elseif parameter_name == 'demand'
        demand = sens_coef.*demand;
    elseif parameter_name == 'prices'
        prices = sens_coef.*prices;
    elseif parameter_name == 'alt_price'
        alt_price = sens_coef.*alt_price;
    elseif parameter_name == 'fixed_costs'
        link_cost = sens_coef.*link_cost;
        station_cost = sens_coef.*station_cost;
    elseif parameter_name == 'linear_costs'
        link_capacity_slope = sens_coef.*link_capacity_slope;
        station_capacity_slope = sens_coef.*station_capacity_slope;
    elseif parameter_name == 'op_cost'
        op_link_cost = sens_coef.*op_link_cost;
    end

    niters = 15;           
    eps = 1e-3;
    a_prev = 1e4.*ones(n);
    s_prev= 1e4.*ones(n,1);
    disp(['beta = ',num2str(beta_or),', lam = ',num2str(lam)]);
    tic;
    for iter=1:niters
        cvx_begin 
            variable s(n)
            variable s_prim(n)
            variable delta_s(n)
            variable a(n,n)
            variable a_prim(n,n)
            variable delta_a(n,n)
            variable f(n,n)
            variable fext(n,n)
            variable fij(n,n,n,n)
            op_obj = 0;
            pax_obj = 0;
            bud_obj = 0;
            bud_obj = bud_obj + 1e-6*sum(station_capacity_slope'.*s_prim);
            bud_obj = bud_obj + 1e-6*(sum(sum(link_capacity_slope.*a_prim)));  %linear construction costs
            op_obj = op_obj + 1e-6*(sum(sum(op_link_cost.*a_prim))); %operation costs
            if iter < niters
                pax_obj = pax_obj + 1e-6*(sum(sum(inv_pos(congestion_coef_links.*delta_a + eps)))) + 1e-6*(sum(inv_pos(congestion_coef_stations'.*delta_s + eps))); %congestion costs
                bud_obj = bud_obj + 1e-6*lam*sum(sum((link_cost.*a_prim.*(1./(a_prev+eps))))) + 1e-6*lam*sum((station_cost'.*s_prim.*(1./(s_prev+eps)))); %fixed construction costs
                
            end
    
            for o=1:n
                for d=1:n
                    pax_obj = pax_obj + 1e-6*(demand(o,d).*sum(sum((travel_time+prices).*fij(:,:,o,d)))); 
                end
            end
            pax_obj = pax_obj + 1e-6*(sum(sum(demand.*(alt_time+alt_price).*fext)));
            pax_obj = pax_obj + 1e-6*(sum(sum(demand.*(-entr(f) - f))));
            pax_obj = pax_obj + 1e-6*(sum(sum(demand.*(-entr(fext) - fext))));
    
    
            if iter == niters
                for i=1:n
                    if s_prev(i) >= eps
                        pax_obj = pax_obj + 1e-6*(sum(inv_pos(congestion_coef_stations(i).*delta_s(i) + eps)));
                    end
                    for j=1:n
                        if a_prev(i,j) >= eps
                            pax_obj = pax_obj + 1e-6*(sum(sum(inv_pos(congestion_coef_links(i,j).*delta_a(i,j) + eps))));
                         end
                     end
                 end
            end
            obj = beta_or*bud_obj + (alfa/(dm_pax))*pax_obj + ((1-alfa)/(dm_op))*op_obj;
            minimize obj
            % constraints
            s >= 0;
            s_prim >= 0;
            delta_s >= 0;
            a >= 0;
            a_prim >= 0;
            delta_a >= 0;
            f >= 0;
            f <= 1;
            fij >= 0;
            fij <= 1;
            fext >= 0;
            fext <= 1;
            s_prim == s + delta_s;
            a_prim == a + delta_a;
    
            for i=1:n
                for j=1:n
                    squeeze(sum(sum(squeeze(permute(fij(i,j,:,:),[3 4 1 2]).*demand)))) <= tau.*a(i,j).*a_nom; %multiplicar por demanda
                end
            end
            for i=1:n
                eta*sum(a(:,i)) <= load_factor(i).*s(i)
            end
            sum(sum(a)) <= a_max;
            for o=1:n
                for d=1:n
                    sum(fij(o,:,o,d)) == f(o,d);
                end
            end
            
            for o=1:n
                squeeze(sum(fij(o,:,o,[1:(o-1),(o+1):n]),2)) - squeeze(sum(permute(fij(:,o,o,[1:(o-1),(o+1):n]),[2,1,3,4]),2)) == transpose(1 - fext(o,[1:(o-1),(o+1):n])); 
            end
            for d=1:n
                squeeze(sum(fij(d,:,[1:(d-1),(d+1):n],d),2)) - squeeze(sum(permute(fij(:,d,[1:(d-1),(d+1):n],d),[2,1,3,4]),2)) == -1 + fext([1:(d-1),(d+1):n],d);
            end
            for i=1:n
                fij(i,i,:,:) == 0;
                squeeze(sum(fij(i,:,[1:(i-1),(i+1):n],[1:(i-1),(i+1):n]),2)) - squeeze(sum(permute(fij(:,i,[1:(i-1),(i+1):n],[1:(i-1),(i+1):n]),[2,1,3,4]),2)) == 0;
            end
            for o=1:n
                fext(o,o) == 0;
                f(o,o) == 0;
                fij(:,o,o,:) == 0;
            end
            for d=1:n
                fij(d,:,:,d) == 0;
            end
            for o=1:n
                for d=1:n
                    if o ~= d
                        f(o,d) + fext(o,d) == 1;
                    end
                end
            end
            for i=1:n
                for j=1:n
                    if ~ismember(j,candidates{i})
                        a_prim(i,j) == 0;
                    end
                end
            end
            if iter == niters
                for i=1:n
                    if s_prev(i) <= 0.1
                        s_prim(i) == 0;
                    end
                    for j=1:n
                        if a_prev(i,j) <= 0.1
                            a_prim(i,j) == 0;
                         end
                     end
                 end
             end
    
        cvx_end
        a_prev = a_prim;
        s_prev = s_prim;
    end
    comp_time = toc;
    
    if parameter_name == 'alt_time'
        alt_time = alt_time./sens_coef;
    elseif parameter_name == 'demand'
        demand = demand./sens_coef;
    elseif parameter_name == 'prices'
        prices = prices./sens_coef;
    elseif parameter_name == 'alt_price'
        alt_price = alt_price./sens_coef;
%     elseif parameter_name == 'fixed_costs'
%         link_cost = link_cost./sens_coef;
%         station_cost = station_cost./sens_coef;
%     elseif parameter_name == 'linear_costs'
%         link_capacity_slope = link_capacity_slope./sens_coef;
%         station_capacity_slope = station_capacity_slope./sens_coef;
    elseif parameter_name == 'op_cost'
        op_link_cost = op_link_cost./sens_coef;
    end
    
    budget = get_budget(s,s_prim,a,a_prim,n,...
        station_cost,station_capacity_slope,link_cost,link_capacity_slope,lam);
    
    [obj_val,pax_obj,op_obj] = get_obj_val(alfa, op_link_cost,...
        congestion_coef_links, ...
        congestion_coef_stations,travel_time,prices,alt_time,alt_price,a_prim,delta_a, ...
        s_prim,delta_s,fij,f,fext,demand,dm_pax,dm_op);


end



function [s,s_prim,delta_s,a,a_prim,delta_a,f,fext,fij,comp_time,budget,obj_val,...
    pax_obj,op_obj] = compute_sim_MIP(lam,nom_bud,alfa,dm_pax,dm_op)

    [n,link_cost,station_cost,link_capacity_slope,...
    station_capacity_slope,demand,prices,...
    load_factor,op_link_cost,congestion_coef_stations,...
    congestion_coef_links,travel_time,alt_time,alt_price,a_nom,tau,eta,...
    a_max,candidates]  = parameters_9node_network();

    M = 1e6;
     
    eps = 1e-3;
    tic;
    
    cvx_begin
        variable s(n)
        variable s_prim(n)
        variable s_bin(n) binary
        variable delta_s(n)
        variable a(n,n)
        variable a_prim(n,n)
        variable a_bin(n,n) binary
        variable delta_a(n,n)
        variable f(n,n)
        variable fext(n,n)
        variable fij(n,n,n,n)
        op_obj = 0;
        pax_obj = 0;
        op_obj = op_obj + 1e-6*(sum(sum(op_link_cost.*a_prim))); %operation costs
        pax_obj = pax_obj + 1e-6*(sum(sum(inv_pos(M.*(1-a_bin) + congestion_coef_links.*delta_a + eps)))) + 1e-6*(sum(inv_pos( M.*(1-s_bin) + congestion_coef_stations'.*delta_s + eps))); %congestion costs
        for o=1:n
            for d=1:n
                pax_obj = pax_obj + 1e-6*(demand(o,d).*sum(sum((travel_time+prices).*fij(:,:,o,d)))); 
            end
        end
        pax_obj = pax_obj + 1e-6*(sum(sum(demand.*(alt_time+alt_price).*fext)));
        pax_obj = pax_obj + 1e-6*(sum(sum(demand.*(-entr(f) - f))));
        pax_obj = pax_obj + 1e-6*(sum(sum(demand.*(-entr(fext) - fext))));
        obj = (alfa/(dm_pax))*pax_obj + ((1-alfa)/(dm_op))*op_obj;
        minimize obj
        % constraints
        s >= 0;
        s_prim >= 0;
        delta_s >= 0;
        a >= 0;
        a_prim >= 0;
        delta_a >= 0;
        f >= 0;
        f <= 1;
        fij >= 0;
        fij <= 1;
        fext >= 0;
        fext <= 1;
        s_prim == s + delta_s;
        a_prim == a + delta_a;
        M.*s_bin >= s_prim;
        M.*a_bin >= a_prim;
        bg = 0;
        bg = bg + sum(station_capacity_slope'.*s_prim);
        bg = bg + sum(sum(link_capacity_slope.*a_prim));
        bg = bg + lam*sum(sum((link_cost.*a_bin)));
        bg = bg + lam*sum(station_cost'.*s_bin);
        bg <= nom_bud;

        for i=1:n
            for j=1:n
                squeeze(sum(sum(squeeze(permute(fij(i,j,:,:),[3 4 1 2]).*demand)))) <= tau.*a(i,j).*a_nom; 
            end
        end
        for i=1:n
            eta*sum(a(:,i)) <= load_factor(i).*s(i)
        end
        sum(sum(a)) <= a_max;
        for o=1:n
            for d=1:n
                sum(fij(o,:,o,d)) == f(o,d);
            end
        end
        
        for o=1:n
            squeeze(sum(fij(o,:,o,[1:(o-1),(o+1):n]),2)) - squeeze(sum(permute(fij(:,o,o,[1:(o-1),(o+1):n]),[2,1,3,4]),2)) == transpose(1 - fext(o,[1:(o-1),(o+1):n])); 
        end
        for d=1:n
            squeeze(sum(fij(d,:,[1:(d-1),(d+1):n],d),2)) - squeeze(sum(permute(fij(:,d,[1:(d-1),(d+1):n],d),[2,1,3,4]),2)) == -1 + fext([1:(d-1),(d+1):n],d);
        end
        for i=1:n
            fij(i,i,:,:) == 0;
            squeeze(sum(fij(i,:,[1:(i-1),(i+1):n],[1:(i-1),(i+1):n]),2)) - squeeze(sum(permute(fij(:,i,[1:(i-1),(i+1):n],[1:(i-1),(i+1):n]),[2,1,3,4]),2)) == 0;
        end
        for o=1:n
            fext(o,o) == 0;
            f(o,o) == 0;
            fij(:,o,o,:) == 0;
        end
        for d=1:n
            fij(d,:,:,d) == 0;
        end
        for o=1:n
            for d=1:n
                if o ~= d
                    f(o,d) + fext(o,d) == 1;
                end
            end
        end
        for i=1:n
            for j=1:n
                if ~ismember(j,candidates{i})
                    a_prim(i,j) == 0;
                end
            end
        end

    cvx_end
    a_prev = a_prim;
    s_prev = s_prim;
    comp_time = toc;
        
    budget = get_budget(s,s_prim,a,a_prim,n,...
        station_cost,station_capacity_slope,link_cost,link_capacity_slope,lam);
    
    [obj_val,pax_obj,op_obj] = get_obj_val(alfa, op_link_cost,...
        congestion_coef_links, ...
        congestion_coef_stations,travel_time,prices,alt_time,alt_price,a_prim,delta_a, ...
        s_prim,delta_s,fij,f,fext,demand,dm_pax,dm_op);
    filename = sprintf('./uthopy_results/sol_MIP_budget=%d_lam=%d_uthopy_alfa=%d.mat',nom_bud,lam,alfa);
    save(filename,'s','s_prim','delta_s', ...
        'a','a_prim','delta_a','f','fext','fij','obj_val','pax_obj','op_obj','comp_time','budget');

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