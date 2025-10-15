* ===============================================================
*  netdesign_demo.gms — Demo runnable con datos de juguete (MOSEK, sin SOS2, con set seg)
* ===============================================================

* ---------- Conjuntos ----------
$include "C:\Users\freal\MATLAB\Projects\untitled\code\3node\param_definition.gms";

*dm_pax = 1e4;
*dm_op = 1e2;


dm_pax = 1e4;
dm_op = 1e2;
* ===============================================================
*  netdesign_demo.gms — Demo runnable con datos de juguete (MOSEK, sin SOS2, con set seg)
* ===============================================================

* ---------- Variables ----------
Variables
    obj, op_obj, pax_obj, xlen(o,d), budvar
;

Positive Variables
    s(i), s_prim(i), delta_s(i)
    a(i,i), a_prim(i,i), delta_a(i,i)
    f(i,i), fext(i,i)
    zij(i,i,i,i)
    ghat(o,d)
    delta(seg,o,d)
;

Binary Variables
    s_bin(i)
    a_bin(i,i)
    z(i,i)
    fij(i,i,i,i)
    alfaseg(seg,o,d)         
;


* Cotas superiores
f.up(i,i)       = 1;
fext.up(i,i)    = 1;
zij.up(i,i,o,d) = 1;

* ---------- Ecuaciones ----------
Equations
    def_op_obj, def_pax_obj, def_obj, def_bud
    sprim_link, aprim_link
    bigM_s, bigM_a
    bud_avail
    link_cap
    station_in_cap, station_out_cap
    tot_links
    fij_od_flow, fij_do_flow, fij_trans
    diag_zero_fo, fij_zero_diag, fij_zero_oo, fij_zero_dd
    split_fd, f_leq_Mz
    cand_zero
    zij_bigM_lb, zij_leq_f, zij_leq_Mfij
    def_xlen, alfaseg_sum, delta_bound,delta_bound_last, xlen_interp, slice_bound, ghat_interp, f_leq_ghat
;


def_bud.. budvar =e= sum(i, station_capacity_slope(i)*s_prim(i))
         + sum((i,j), link_capacity_slope(i,j)*a_prim(i,j))
         + lam*sum((i,j),
                        link_cost(i,j)*a_bin(i,j))
         + lam*sum(i,
                        station_cost(i)*s_bin(i));

* ---------- Objetivo ----------
def_op_obj..
    op_obj =e=  sum((i,j), op_link_cost(i,j) * a_prim(i,j));

def_pax_obj..
    pax_obj =e=  sum((o,d), demand(o,d) * fext(o,d));

def_obj..
    obj =e=  ( (alfa/dm_pax)*pax_obj + ((1-alfa)/dm_op)*op_obj );
    


* ---------- Restricciones de infraestructura ----------
sprim_link(i)..          s_prim(i)   =e= s(i) + delta_s(i);
aprim_link(i,j)..        a_prim(i,j) =e= a(i,j) + delta_a(i,j);

bigM_s(i)..              s_prim(i)   =l= M * s_bin(i);
bigM_a(i,j)..            a_prim(i,j) =l= M * a_bin(i,j);

bud_avail..
    sum(i, station_capacity_slope(i)*s_prim(i))
         + sum((i,j), link_capacity_slope(i,j)*a_prim(i,j))
         + lam*sum((i,j),
                        link_cost(i,j)*a_bin(i,j))
         + lam*sum(i,
                        station_cost(i)*s_bin(i)) =l= budget;

link_cap(i,j)..
     sum((o,d), zij(i,j,o,d) * demand(o,d))*tau =l= a(i,j) * a_nom;

station_in_cap(i)..   sigma * sum(j, a(j,i)) =l= s(i);
station_out_cap(i)..  sigma * sum(j, a(i,j)) =l= s(i);

tot_links.. sum((i,j), a(i,j)) =l= a_max;

* ---------- Restricciones de flujo ----------
fij_od_flow(o,d)$(ord(o)<>ord(d))..
    sum(j$(ord(j)<>ord(o)), fij(o,j,o,d))
  - sum(i$(ord(i)<>ord(o)), fij(i,o,o,d))
  =e= z(o,d);

fij_do_flow(o,d)$(ord(o)<>ord(d))..
    sum(j$(ord(j)<>ord(d)), fij(d,j,o,d))
  - sum(i$(ord(i)<>ord(d)), fij(i,d,o,d))
  =e= -z(o,d);

fij_trans(i,o,d)$(ord(i)<>ord(o) and ord(i)<>ord(d))..
    sum(j$(ord(j)<>ord(i)), fij(i,j,o,d))
  - sum(ii$(ord(ii)<>ord(i)), fij(ii,i,o,d))
  =e= 0;

diag_zero_fo(o)..            f(o,o) + fext(o,o) =e= 0;
fij_zero_diag(i,o,d)..       fij(i,i,o,d) =e= 0;
fij_zero_oo(o,i,d)..         fij(i,o,o,d) =e= 0;
fij_zero_dd(d,i,o)..         fij(d,i,o,d) =e= 0;

split_fd(o,d)$(ord(o)<>ord(d))..  f(o,d) + fext(o,d) =e= 1;
f_leq_Mz(o,d)$(ord(o)<>ord(d))..  f(o,d) =l=  z(o,d);

cand_zero(i,j)$(candidates(i,j)=0)..  a_prim(i,j) =e= 0;

zij_bigM_lb(i,j,o,d)..  zij(i,j,o,d) =g= f(o,d) - (1 - fij(i,j,o,d));
zij_leq_f(i,j,o,d)..    zij(i,j,o,d) =l= f(o,d);
zij_leq_Mfij(i,j,o,d).. zij(i,j,o,d) =l= fij(i,j,o,d);

* ---------- PWL con convex combination + binarios ----------
def_xlen(o,d)$(ord(o)<>ord(d))..
    xlen(o,d) =e= sum((i,j), fij(i,j,o,d) * logit_coef*(-travel_time(i,j) - prices(i,j)));
    
alfaseg_sum(o,d)$(ord(o)<>ord(d)).. sum(seg,alfaseg(seg,o,d)) =e= 1;

delta_bound(seg,o,d)$((ord(o)<>ord(d)) and (ord(seg) lt card(seg))).. delta(seg,o,d) =l= alfaseg(seg,o,d)*(b(seg+1,o,d)-b(seg,o,d));

delta_bound_last(seg,o,d)$(ord(seg)=card(seg))..  delta(seg,o,d) =l= alfaseg(seg,o,d)*(-b(seg,o,d));

xlen_interp(o,d)$(ord(o)<>ord(d)).. xlen(o,d) =e= sum(seg,alfaseg(seg,o,d)*b(seg,o,d) + delta(seg,o,d));

slice_bound(seg,o,d)$((ord(seg) lt card(seg)) and (ord(o)<>ord(d)) ).. b(seg,o,d) + delta(seg,o,d) =l= b(seg+1,o,d);

ghat_interp(o,d)$(ord(o)<>ord(d)).. ghat(o,d) =e= sum(seg,alfaseg(seg,o,d)*bord(seg,o,d) + mreg(seg,o,d)*delta(seg,o,d));   

f_leq_ghat(o,d)$(ord(o)<>ord(d))..     f(o,d) =l= ghat(o,d);


* ---------- Modelo y resolución ----------
Model netdesign /
*    def_bud
    def_op_obj, def_pax_obj, def_obj
    sprim_link, aprim_link
    bigM_s, bigM_a
    bud_avail
    link_cap
    station_in_cap, station_out_cap
    tot_links
    fij_od_flow, fij_do_flow, fij_trans
    diag_zero_fo, fij_zero_diag, fij_zero_oo, fij_zero_dd
    split_fd, f_leq_Mz
    cand_zero
    zij_bigM_lb, zij_leq_f, zij_leq_Mfij 
    def_xlen
    alfaseg_sum
    delta_bound
    delta_bound_last
    slice_bound
    xlen_interp
    ghat_interp, f_leq_ghat
/;


option threads = 60;
option mip     = mosek;

option reslim = 1800;

Parameter mipgap;

Solve netdesign using mip minimizing obj;

* ---------- Salidas básicas ----------
*Display f.l, fext.l, fij.l, xlen.l, ghat.l,
*        obj.l, op_obj.l, pax_obj.l,
*        s.l, s_bin.l, a.l, a_bin.l, z.l, b;

mipgap = abs(netdesign.objval - netdesign.objest) / (1e-9 + abs(netdesign.objval));


Parameter solverTime;

solverTime = netdesign.resusd;

Display b,budget,netdesign.objval, netdesign.objest, mipgap,solverTime;


Parameter fnew(o,i,j,d);

fnew(o,i,j,d)=fij.l(i,j,o,d);
execute_unload "resultado.gdx";


file fijx /'fij_long.csv'/; put fijx;
put 'i,j,o,d,value' /;
loop((i,j,o,d),
    put i.tl ',' j.tl ',' o.tl ',' d.tl ',' fij.l(i,j,o,d):0:15 / );
putclose fijx;

EmbeddedCode Connect:
- GAMSReader:
    symbols:
      - name: a_prim
- Projection:
    name: a_prim.l(i,j)
    newName: a_prim_level(i,j)
- ExcelWriter:
    file: output_a.xlsx
    symbols:
      - name: a_prim_level
endEmbeddedCode

EmbeddedCode Connect:
- GAMSReader:
    symbols:
      - name: s_prim
- Projection:
    name: s_prim.l(i)
    newName: s_prim_level(i)
- ExcelWriter:
    file: output_s.xlsx
    symbols:
      - name: s_prim_level
endEmbeddedCode

EmbeddedCode Connect:
- GAMSReader:
    symbols:
      - name: s_prim
      - name: s
      - name: delta_s
      - name: a_prim
      - name: a
      - name: delta_a
      - name: f
      - name: fext
      - name: fnew
      - name: mipgap
      - name: solverTime
- Projection:
    name: s_prim.l(i)
    newName: sprim_level(i)
- Projection:
    name: s.l(i)
    newName: s_level(i)
- Projection:
    name: delta_s.l(i)
    newName: deltas_level(i)
- Projection:
    name: a_prim.l(i,j)
    newName: aprim_level(i,j)
- Projection:
    name: a.l(i,j)
    newName: a_level(i,j)
- Projection:
    name: delta_a.l(i,j)
    newName: deltaa_level(i,j)
- Projection:
    name: f.l(o,d)
    newName: f_level(o,d)
- Projection:
    name: fext.l(o,d)
    newName: fext_level(o,d)
- Projection:
    name: fnew(o,i,j,d)
    newName: fij_level(o,i,j,d)
- Projection:
    name: mipgap
    newName: mip_opt_gap
- Projection:
    name: solverTime
    newName: solver_time
- ExcelWriter:
    file: output_all.xlsx
    symbols:
      - name: sprim_level
      - name: s_level
      - name: deltas_level
      - name: aprim_level
      - name: a_level
      - name: deltaa_level
      - name: f_level
      - name: fext_level
      - name: fij_level
      - name: mip_opt_gap
      - name: solver_time
endEmbeddedCode