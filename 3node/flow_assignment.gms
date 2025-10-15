* ===============================================================
*  netdesign_demo.gms — Demo runnable con datos de juguete (MOSEK, sin SOS2, con set seg)
* ===============================================================

* ---------- Conjuntos ----------
$include "C:\Users\freal\MATLAB\Projects\untitled\code\3node\param_definition.gms";

*dm_pax = 1e4;
*dm_op = 1e2;



Parameter a_prim_prev(i,j)
/
$include "%TXTDIR%\a_prim_mipreg.txt"
/;



dm_pax = 1e4;
dm_op = 1e2;

* ---------- Variables ----------
Variables
    obj, pax_obj, xlen(o,d)
;

Positive Variables
    f(i,i), fext(i,i)
    zij(i,i,i,i)
    ghat(o,d)
    delta(seg,o,d)
;

Binary Variables
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
    def_pax_obj, def_obj
    link_cap
    fij_od_flow, fij_do_flow, fij_trans
    diag_zero_fo, fij_zero_diag, fij_zero_oo, fij_zero_dd
    split_fd, f_leq_Mz
    zij_bigM_lb, zij_leq_f, zij_leq_Mfij
    def_xlen, alfaseg_sum, delta_bound,delta_bound_last, xlen_interp, slice_bound, ghat_interp, f_leq_ghat
;




* ---------- Objetivo ----------
def_pax_obj..
    pax_obj =e=  sum((o,d), demand(o,d) * fext(o,d));

def_obj..
    obj =e=  pax_obj;
    
* ---------- Restricciones de infraestructura ----------
link_cap(i,j)..
     sum((o,d), zij(i,j,o,d) * demand(o,d))*tau =l= a_prim_prev(i,j) * a_nom;

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
    def_pax_obj, def_obj
    fij_od_flow, fij_do_flow, fij_trans
    diag_zero_fo, fij_zero_diag, fij_zero_oo, fij_zero_dd
    link_cap
    split_fd, f_leq_Mz
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

Display b,budget,netdesign.objval, netdesign.objest, mipgap;


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
      - name: f
      - name: fext
      - name: fnew
- Projection:
    name: f.l(o,d)
    newName: f_level(o,d)
- Projection:
    name: fext.l(o,d)
    newName: fext_level(o,d)
- Projection:
    name: fnew(o,i,j,d)
    newName: fij_level(o,i,j,d)
- ExcelWriter:
    file: output_all.xlsx
    symbols:
      - name: f_level
      - name: fext_level
      - name: fij_level
endEmbeddedCode