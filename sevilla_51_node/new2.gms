


$setglobal TXTDIR "C:\Users\freal\MATLAB\Projects\untitled\code\sevilla_net_49\export_txt"




Set i "nodos" /i1*i51/;
Alias (i,j,o,d);
* === Auxiliares (decláralos UNA sola vez en tu código) ===
*


Parameter demand(o,d)

/
$include "%TXTDIR%\demand.txt"
/;


* --- (o,d) ---
Parameter alt_price(o,d)
/
$include "%TXTDIR%\alt_price.txt"
/;

* --- (i,j) ---
Parameter link_cost(i,j)
/
$include "%TXTDIR%\link_cost.txt"
/;

Parameter link_capacity_slope(i,j)
/
$include "%TXTDIR%\link_capacity_slope.txt"
/;

Parameter prices(i,j)
/
$include "%TXTDIR%\prices.txt"
/;

Parameter op_link_cost(i,j)
/
$include "%TXTDIR%\op_link_cost.txt"
/;

Parameter congestion_coefs_links(i,j)
/
$include "%TXTDIR%\congestion_coefs_links.txt"
/;

Parameter candidates(i,j)
/
$include "%TXTDIR%\candidates.txt"
/;



* --- (i) ---
Parameter station_cost(i)
/
$include "%TXTDIR%\station_cost.txt"
/;

Parameter station_capacity_slope(i)
/
$include "%TXTDIR%\station_capacity_slope.txt"
/;

Parameter congestion_coefs_stations(i)
/
$include "%TXTDIR%\congestion_coefs_stations.txt"
/;


Parameter a_prev(i,j)
/
$include "%TXTDIR%\a_prev.txt"
/;

Parameter s_prev(i)
/
$include "%TXTDIR%\s_prev.txt"
/;


Scalars tau, sigma, a_nom, a_max;
tau = 0.57;
sigma = 0.25;
a_nom = 588;
a_max = 1e5;

Scalar n; n = card(i);

Scalar epsi; epsi = 1e-3;

Scalar lam /
$include "%TXTDIR%\lam.txt"
/;

Scalar alfa /
$include "%TXTDIR%\alfa.txt"
/;

Scalar beta /
$include "%TXTDIR%\beta.txt"
/;

Scalar iter /
$include "%TXTDIR%\current_iter.txt"
/;

Scalar niters /
$include "%TXTDIR%\niters.txt"
/;

Scalar logit_coef; logit_coef = 0.2;



display link_cost,link_capacity_slope,prices,congestion_coefs_links,candidates,station_cost,station_capacity_slope,congestion_coefs_stations,a_max;

