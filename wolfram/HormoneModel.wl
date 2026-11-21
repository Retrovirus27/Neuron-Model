 (* ::Package:: *)

(* HormoneModel.wl - modelo hormonal desacoplado de NeuronModel.nb

   Uso en Jupyter:  Get["~/Downloads/HormoneModel.wl"]  *)




(* ::Section:: *)
(* Modelo hormonal (desacoplado del modelo neuronal) *)


(* ::Section:: *)
(* 1. Funciones de tiempo *)

gen[x_] := x*60*60;            (* horas -> segundos *)
num[x_] := x/(60*60);          (* tasa por hora -> por segundo *)
break[x_] := (x*60*60) + 2340; (* fin del desayuno: dura 2340 s *)
lun[x_]   := (x*60*60) + 1676; (* fin del almuerzo: dura 1676 s *)
din[x_]   := (x*60*60) + 1400; (* fin de la cena: dura 1400 s *)


(* ::Section:: *)
(* 2. Parámetros *)

(* Glucosa y ghrelina *)
{k12, k2x, sigma, khx, beta, beta2, gamma} =
  {num[2.30], num[0.71], 0.20, num[2.354], num[1591], num[1238], 1.29};

(* Insulina *)
{theta, delta, s, kix, rho} = {num[528.63], 0.38, 4, num[1.03], num[90]};

(* Leptina *)
{epsilon, FM, BV} = {num[6231.39], 60000, 5000};
{R, eta, l, EC50LEP} = {num[0.11], 24, 1, 661};


(* ::Section:: *)
(* 3. Entradas: patrones de alimentación *)


(* ::Section:: *)
(* Alimentación con horario restringido (TRE): 3 comidas a las 8, 13 y 17 h *)

TRE[t_, r_] := Which[
   gen[8] < t < break[8], r, gen[13] < t < lun[13], r, gen[17] < t < din[17], r,
   gen[32] < t < break[32], r, gen[37] < t < lun[37], r, gen[41] < t < din[41], r,
   gen[56] < t < break[56], r, gen[61] < t < lun[61], r, gen[65] < t < din[65], r,
   True, 0];


(* ::Section:: *)
(* Dieta común (CD): 3 comidas (`r`) + 2 meriendas (`s2`) a las 16 y 20.5 h *)

CD[t_, r_, s2_] := Which[
   gen[8] < t < break[8], r, gen[13] < t < lun[13], r, gen[16] < t < din[16], s2,
   gen[19] < t < din[19], r, gen[20.5] < t < din[20.5], s2,
   gen[32] < t < break[32], r, gen[37] < t < lun[37], r, gen[40] < t < din[40], s2,
   gen[43] < t < din[43], r, gen[44.5] < t < din[44.5], s2,
   gen[56] < t < break[56], r, gen[61] < t < lun[61], r, gen[64] < t < din[64], s2,
   gen[67] < t < din[67], r, gen[68.5] < t < din[68.5], s2,
   True, 0];


(* ::Section:: *)
(* Ayuno intermitente (IF): 2 comidas, sin desayuno (13 y 17 h) *)

IFI[t_, r_] := Which[gen[13] < t < lun[13], r, gen[17] < t < din[17], r, True, 0];


(* ::Section:: *)
(* 4. Ecuaciones del modelo hormonal *)


(* ::Section:: *)
(* Ec. 1: Glucosa, compartimento 1 (absorción) *)

dGL1[gl1_, gl2_] := -(gl1 k12)/(1 + sigma gl2);


(* ::Section:: *)
(* Ec. 2: Glucosa, compartimento 2 (plasma) *)

dGL2[gl1_, gl2_] := (k12/(1 + sigma gl2)) gl1 - k2x gl2;


(* ::Section:: *)
(* Ec. 3: Ghrelina *)

betaT[t_] := If[gen[0] <= t < gen[3] || gen[6] <= t < gen[27] ||
     gen[30] <= t < gen[51] || gen[54] <= t < gen[75] || gen[78] <= t < gen[99] ||
     gen[102] <= t < gen[123] || gen[126] <= t < gen[147] || gen[150] <= t < gen[171],
    beta, beta2];

dGHR[ghr_, gl2_, t_] := betaT[t]/(1 + gamma gl2) - khx ghr;


(* ::Section:: *)
(* Ec. 4: Insulina *)

dINS[gl2_, ins_] := rho + (theta gl2^s)/(delta^s + gl2^s) - kix ins;


(* ::Section:: *)
(* Ec. 5: Leptina *)

dLEP[lep_, ins_] := (FM (R eta)/(1 + (ins/EC50LEP)^l) - epsilon lep)/BV;


(* ::Section:: *)
(* 5. Sistema completo *)

hormoneModel[food_] := {
   GL1'[t] == dGL1[GL1[t], GL2[t]] + food[t],
   GL2'[t] == dGL2[GL1[t], GL2[t]],
   GHR'[t] == dGHR[GHR[t], GL2[t], t],
   INS'[t] == dINS[GL2[t], INS[t]],
   LEP'[t] == dLEP[LEP[t], INS[t]],
   GL1[0] == 0, GL2[0] == 0.01, GHR[0] == 670, INS[0] == 79, LEP[0] == 22};

vars = {GL1, GL2, GHR, INS, LEP};


(* ::Section:: *)
(* Resolver *)

solveHormones[food_, hours_: 72] :=
  First@NDSolve[hormoneModel[food], vars, {t, 0, gen[hours]}, MaxStepSize -> 60];

solTRE = solveHormones[TRE[#, num[2.1]] &];
solCD  = solveHormones[CD[#, num[2.1], num[0.8]] &];
solIF  = solveHormones[IFI[#, num[2.1]] &];


(* ::Section:: *)
(* 6. Datos experimentales (del notebook original) *)

datainsulin = {{0,120.83},{1,113.24},{2,106.86},{3,103.65},{4,90.16},{5,79.23},{6,79.44},{7,89.17},{8,89.17},{8.5,71.34},{9.5,391.79},{10,518.51},{10.5,502.54},{11,373.14},{11.5,331.01},{12,214.13},{12.5,168.69},{13,384.29},{13.5,341.9},{14.5,248.41},{15,244.37},{15.5,154.12},{16,199.26},{16.5,155.32},{17,170.29},{17.5,126.07},{18,144.55},{18.5,280.14},{19,382.17},{19.5,393.84},{20,349.66},{20.5,292.61},{21,222.61},{21.5,188.44},{22,159.02},{23,129.76}};

dataleptin = {{0,22.8286},{1,23.4473},{2,23.503},{3,23.5769},{4,23.9776},{5,23.9425},{6,23.2717},{7,22.6371},{8,22.5111},{8.5,21.0934},{9,19.4697},{9.5,17.1366},{10,16.0502},{10.5,15.6503},{11,15.9026},{11.5,16.3265},{12,16.6017},{12.5,17.3575},{13,17.0149},{13.5,17.9537},{14,16.3524},{14.5,17.7032},{15,17.7882},{15.5,18.0114},{16,18.506},{16.5,19.7699},{17,19.8795},{17.5,20.316},{18,20.7525},{18.5,21.298},{19,20.6447},{19.5,19.6463},{20,19.774},{20.5,20.0834},{21,20.2112},{22,21.1926},{23,21.393}};

dataghrelin = {{0,625.29},{1,650.48},{2,626.89},{3,613.66},{4,542.79},{5,504.43},{6,501.48},{7,600.55},{8,704.07},{8.5,525.29},{9,321.37},{9.5,343.57},{10,352.47},{10.5,345.11},{11,389.49},{11.5,450.11},{12,648.17},{12.5,689.58},{13,432.47},{13.5,451.71},{14,422.18},{14.5,408.92},{15,502.05},{15.5,484.35},{16,509.51},{16.5,540.58},{17,691.34},{17.5,704.68},{18,639.68},{18.5,514.1},{19,515.61},{19.5,450.62},{20,492.03},{20.5,475.81},{21,602.93},{22,608.91},{23,608.97}};

dataleptininpM  = dataleptin  /. {x_, y_} :> {x, y/2.0043};
dataghrelininpM = dataghrelin /. {x_, y_} :> {x, y/3.3709};


(* ::Section:: *)
(* 7. Gráficas *)

plotHormones[sol_, food_, hours_, title_, data_: False, meals_: {}] := Module[{opts, row},
  opts = {Frame -> True, ImageSize -> Medium, PlotStyle -> Black,
    LabelStyle -> {10, Black, Bold}, PlotTheme -> "Scientific",
    GridLines -> {meals, None}, GridLinesStyle -> Directive[Gray, Thick]};
  row[expr_, label_, range_, ar_, d_] := Show[
    Plot[Evaluate[expr /. sol], {h, 0, hours}, PlotRange -> range,
     FrameLabel -> {None, label}, AspectRatio -> ar, Evaluate[Sequence @@ opts]],
    If[d === None || data === False, Graphics[{}],
     ListPlot[d, PlotStyle -> {Gray, PointSize[Medium]}]]];
  Column[{
    Plot[food[gen[h]]/num[1], {h, 0, hours}, PlotRange -> {0, 2.5}, Filling -> Bottom,
     Frame -> True, FrameLabel -> {None, "Food Intake (L/h)"}, PlotLabel -> title,
     AspectRatio -> 1/4, ImageSize -> Medium, LabelStyle -> {10, Black, Bold},
     PlotPoints -> 400],
    row[GL1[gen[h]], "GI\:2081 (L)", {0, 0.8}, 1/7, None],
    row[GL2[gen[h]], "GI\:2082 (L)", {0, 0.8}, 1/7, None],
    row[GHR[gen[h]]/3.3709, "Ghrelin (pM)", {0, 300}, 1/6, dataghrelininpM],
    row[INS[gen[h]], "Insulin (pM)", {0, 600}, 1/6, datainsulin],
    row[LEP[gen[h]]/2.0043, "Leptin (pM)", {0, 20}, 1/8, dataleptininpM]}]];


(* ::Section:: *)
(* 24 h *)


(* ::Section:: *)
(* 72 h *)


(* ::Section:: *)
(* 8. Exportar una figura (opcional) *)


(* ::Section:: *)
(* Anexo: salida hacia las neuronas (no forma parte del modelo hormonal) *)

{tauactINSR, tauinacINSR, EC50INSR} = {300, 300, 661};
{tauactGHSR, tauinacGHSR, EC50GHSR} = {300, 300, 500};
{tauactLEPR, tauinacLEPR, EC50LEPR} = {300, 300, 16};

receptorModel[sol_] := Module[{ins, ghr, lep},
  {ins, ghr, lep} = {INS, GHR, LEP} /. sol;
  First@NDSolve[{
     AINSR'[t] == (ins[t]/(EC50INSR + ins[t]))/tauactINSR - AINSR[t]/tauinacINSR,
     AGHSR'[t] == (ghr[t]/(EC50GHSR + ghr[t]))/tauactGHSR - AGHSR[t]/tauinacGHSR,
     ALEPR'[t] == (lep[t]/(EC50LEPR + lep[t]))/tauactLEPR - ALEPR[t]/tauinacLEPR,
     AINSR[0] == 0.11, AGHSR[0] == 0.6, ALEPR[0] == 0.6},
    {AINSR, AGHSR, ALEPR}, {t, 0, gen[72]}, MaxStepSize -> 60]];

(* recTRE = receptorModel[solTRE];
   Plot[Evaluate[{AINSR[gen[h]], AGHSR[gen[h]], ALEPR[gen[h]]} /. recTRE], {h, 0, 24},
     PlotLegends -> {"INSR", "GHSR", "LEPR"}] *)
