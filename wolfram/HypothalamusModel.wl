(* ::Package:: *)
(* ============================================================================
   HypothalamusModel.wl  -  Modelo de referencia del sistema hipotalamico
   Base: Tamayo-Molina et al. (2024) Sci Rep 14:28897, codigo reference/NeuronModel.nb

   CAMBIOS respecto al codigo publicado (ver docs/CHANGELOG.md, CH-20261120-04):
     1. Tiempo en HORAS (antes segundos). Tabla de conversion: docs/UNIDADES.md
        y la variable unitTable de este archivo.
     2. Correccion: corriente de leptina en POMC con el voltaje de POMC
        (antes ILEPR[A[t],...]).
     3. Una sola funcion modelEqs[diet] (antes 4 copias del modelo).
     4. Horarios periodicos con Mod[t, 24]: comidas y beta de la ghrelina.
        Se puede simular cualquier numero de dias.
        (En el original, IFI solo tenia comidas el dia 1.)
     5. Nombres seguros: sin I como argumento, sin h[h_,...], letras griegas en ASCII.
   Se MANTIENEN las formas ajustadas a los datos humanos (no las del texto):
     leptina ^(+l), beta2 de la ghrelina entre 3 y 6 h, fase circadiana -Pi/3.

   Uso:  SetDirectory["~/Documents/22_Neuron_model"];
         Get["wolfram/HypothalamusModel.wl"]
   ========================================================================== *)

(* ::Section:: *)
(* 0. Unidades *)

sec = 1/3600.;          (* 1 segundo expresado en horas *)
Cm  = sec;              (* capacitancia efectiva: el original usaba C = 1 con t en s *)

(* Fondo blanco en todas las graficas (JupyterLab con tema oscuro las mostraba
   con fondo transparente y las lineas negras no se veian) *)
SetOptions[{Plot, ListPlot, ListLinePlot, ParametricPlot, StreamPlot, Graphics,
   GraphicsRow, GraphicsColumn, GraphicsGrid}, Background -> White];

(* ::Section:: *)
(* 1. Parametros (t en horas) *)

(* --- Corrientes principales --- *)
{gL, EL, gNa, ENa, sNa, v0Na, gK, EK, sK, v0K} = {0.1, -60, 3, 50, 0.25, -25, 4, -90, 0.25, -25};
tauK = 2 sec;

(* --- Ionotropicas --- *)
{gGLU, EGLU, gGABA, EGABA, sSyn, vSpike} = {0.1, 50, 0.1, -65, 1, -20};
tauGLU = 30 sec;  tauGABA = 30 sec;

(* --- Orexina --- *)
{gOX, EOX, Mmax} = {0.2, 50, 1};
tauOX = 300 sec;  tauInc = 27000 sec;  tauDec = 3312 sec;

(* --- MC4R (modelo de dos estados) --- *)
{gMC4R, EMC4R, basal, m} = {0.1, -38, 0.5, 1};
tauAcoMSH = 400 sec; tauEliMSH = 300 sec; tauAcoAgRP = 300 sec; tauEliAgRP = 300 sec;
tauActMC4R = 734.4 sec; tauInacMC4R = 2030.4 sec;
maxMSH = 1.*^-9; maxAgRP = 1.*^-8;
KK[kapp_, mx_] := kapp basal/mx;
kk[kapp_, mx_] := KK[kapp, mx] ((1 - basal)/basal) mx/(1 - mx);
TA = KK[3.6*^-9, 0.1];  ZA = kk[3.6*^-9, 0.1];
Ta = KK[0.16*^-9, 0.9]; Za = kk[0.16*^-9, 0.9];

(* --- Hormonas (tasas en 1/h: ya no hace falta num[]) --- *)
{k12, k2x, sigma, khx, beta, beta2, gamma} = {2.30, 0.71, 0.20, 2.354, 1591, 1238, 1.29};
{theta, delta, s, kix, rho} = {528.63, 0.38, 4, 1.03, 90};
{epsilon, FM, BV} = {6231.39, 60000, 5000};
{R0, eta, l, EC50LEP} = {0.11, 24, 1, 661};

(* --- Receptores hormonales --- *)
{gINSR, EC50INSR, gGHSR, EC50GHSR, gLEPR, EC50LEPR} = {0.1, 661, 0.1, 500, 0.045, 16};
tauActINSR = 300 sec; tauInacINSR = 300 sec;
tauActGHSR = 300 sec; tauInacGHSR = 300 sec;
tauActLEPR = 300 sec; tauInacLEPR = 300 sec;
(* Potenciales de inversion de los receptores (iguales al original) *)
{EGHSRagrp, ELEPRagrp, EINSRagrp, ELEPRpomc, EINSRpomc} = {-10, -90, -90, -20, -10};

(* --- Ritmo circadiano --- *)
wc = 0.000073*3600;      (* rad/h, = 0.2628; el original usaba 0.000073 rad/s *)
Ac = 0.045;

(* --- Conversiones para graficar en pM --- *)
ghrToPM = 1/3.3709;  lepToPM = 1/2.0043;

(* ::Section:: *)
(* 2. Funciones del modelo *)

sig[v_, slope_, v0_] := 1/(1 + Exp[-slope (v - v0)]);
Ipc[v_, aK_] := gL (v - EL) + gNa sig[v, sNa, v0Na] (v - ENa) + gK aK (v - EK);
Icurr[g_, a_, v_, erev_] := g a (v - erev);

Icirc[t_] := Ac (0.97 Sin[wc t - Pi/3] + 0.22 Sin[2 wc t - Pi/3] + 0.07 Sin[3 wc t - Pi/3] +
     0.03 Sin[4 wc t - Pi/3] + 0.001 Sin[5 wc t - Pi/3] + 1);

RaMC4R[msh_, agrp_] := (msh + Ta (1 + agrp/TA))/
   (Ta ((1 + m) + ((1/TA) + (m/ZA)) agrp) + (1 + (m Ta/Za)) msh);

hill1[c_, ec50_] := c/(ec50 + c);

(* Ghrelina: produccion basal beta2 entre las 3 y las 6 h de cada dia *)
betaT[t_?NumericQ] := If[3 <= Mod[t, 24] < 6, beta2, beta];

(* ::Section:: *)
(* 3. Dietas (periodicas cada 24 h) *)
(* {hora de inicio, duracion en h, tasa de ingesta en L/h} *)

diets = <|
   "TRE" -> {{8, 2340 sec, 2.1}, {13, 1676 sec, 2.1}, {17, 1400 sec, 2.1}},
   "SM"  -> {{8, 2340 sec, 1.1}, {13, 1676 sec, 1.1}, {17, 1400 sec, 1.1}},
   "CD"  -> {{8, 2340 sec, 2.1}, {13, 1676 sec, 2.1}, {16, 1400 sec, 0.8},
             {19, 1400 sec, 2.1}, {20.5, 1400 sec, 0.8}},
   "IF"  -> {{13, 1676 sec, 2.1}, {17, 1400 sec, 2.1}}|>;

food[diet_List][t_?NumericQ] := With[{hd = Mod[t, 24]},
   Total[If[#[[1]] < hd < #[[1]] + #[[2]], #[[3]], 0] & /@ diet]];

(* ::Section:: *)
(* 4. Ecuaciones *)

hormoneVars = {GL1, GL2, GHR, INS, LEP};
receptorVars = {AINSR, AGHSR, ALEPR};
neuronVars = {V, A, P, L, aKV, aKA, aKP, aKL, aGLUP, aGLUV, GABAA, GABAL,
   aOX, MOX, aMSH, AgRP, AMC4R};
vars = Join[neuronVars, hormoneVars, receptorVars];

hormoneEqs[diet_List] := {
   GL1'[t] == -(k12 GL1[t])/(1 + sigma GL2[t]) + food[diet][t],
   GL2'[t] == (k12/(1 + sigma GL2[t])) GL1[t] - k2x GL2[t],
   GHR'[t] == betaT[t]/(1 + gamma GL2[t]) - khx GHR[t],
   INS'[t] == rho + (theta GL2[t]^s)/(delta^s + GL2[t]^s) - kix INS[t],
   LEP'[t] == (FM R0 eta/(1 + (INS[t]/EC50LEP)^l) - epsilon LEP[t])/BV};

receptorEqs = {
   AINSR'[t] == hill1[INS[t], EC50INSR]/tauActINSR - AINSR[t]/tauInacINSR,
   AGHSR'[t] == hill1[GHR[t], EC50GHSR]/tauActGHSR - AGHSR[t]/tauInacGHSR,
   ALEPR'[t] == hill1[LEP[t], EC50LEPR]/tauActLEPR - ALEPR[t]/tauInacLEPR};

neuronEqs = {
   (* PVH *)
   V'[t] == (-Ipc[V[t], aKV[t]] - Icurr[gGLU, aGLUP[t], V[t], EGLU] -
       Icurr[gGABA, GABAA[t], V[t], EGABA] - Icurr[gGABA, GABAL[t], V[t], EGABA] -
       Icurr[gMC4R, AMC4R[t], V[t], EMC4R])/Cm,
   (* AgRP *)
   A'[t] == (-Ipc[A[t], aKA[t]] - Icurr[gGLU, aGLUV[t], A[t], EGLU] -
       Icurr[gOX, aOX[t], A[t], EOX] - Icurr[gGHSR, AGHSR[t], A[t], EGHSRagrp] -
       Icurr[gLEPR, ALEPR[t], A[t], ELEPRagrp] - Icurr[gINSR, AINSR[t], A[t], EINSRagrp])/Cm,
   (* POMC  (CORREGIDO: la corriente de leptina usa P[t]) *)
   P'[t] == (-Ipc[P[t], aKP[t]] - Icurr[gGABA, GABAA[t], P[t], EGABA] +
       Icurr[gOX, aOX[t], P[t], EOX] - Icurr[gINSR, AINSR[t], P[t], EINSRpomc] -
       Icurr[gLEPR, ALEPR[t], P[t], ELEPRpomc])/Cm,
   (* LHA *)
   L'[t] == (-Ipc[L[t], aKL[t]] - Icurr[gMC4R, AMC4R[t], L[t], EMC4R] + Icirc[t])/Cm,
   (* Potasio *)
   aKV'[t] == (sig[V[t], sK, v0K] - aKV[t])/tauK,
   aKA'[t] == (sig[A[t], sK, v0K] - aKA[t])/tauK,
   aKP'[t] == (sig[P[t], sK, v0K] - aKP[t])/tauK,
   aKL'[t] == (sig[L[t], sK, v0K] - aKL[t])/tauK,
   (* Sinapsis ionotropicas *)
   aGLUP'[t] == (sig[P[t], sSyn, vSpike] - aGLUP[t])/tauGLU,
   aGLUV'[t] == (sig[V[t], sSyn, vSpike] - aGLUV[t])/tauGLU,
   GABAA'[t] == (sig[A[t], sSyn, vSpike] - GABAA[t])/tauGABA,
   GABAL'[t] == (sig[L[t], sSyn, vSpike] - GABAL[t])/tauGABA,
   (* Orexina *)
   aOX'[t] == (sig[L[t], sSyn, vSpike] MOX[t] - aOX[t])/tauOX,
   MOX'[t] == (Mmax - MOX[t])/tauInc - sig[L[t], sSyn, vSpike] MOX[t]/tauDec,
   (* Neuropeptidos en la hendidura y MC4R *)
   aMSH'[t] == maxMSH sig[P[t], sSyn, vSpike]/tauAcoMSH - aMSH[t]/tauEliMSH,
   AgRP'[t] == maxAgRP sig[A[t], sSyn, vSpike]/tauAcoAgRP - AgRP[t]/tauEliAgRP,
   AMC4R'[t] == RaMC4R[aMSH[t], AgRP[t]]/tauActMC4R - AMC4R[t]/tauInacMC4R};

initialConditions = {
   V[0] == -40, A[0] == -60, P[0] == -60, L[0] == -60,
   aKV[0] == 0, aKA[0] == 0, aKP[0] == 0, aKL[0] == 0,
   aGLUP[0] == 0, aGLUV[0] == 0, GABAA[0] == 0, GABAL[0] == 0,
   aOX[0] == 0, MOX[0] == 1, aMSH[0] == 0, AgRP[0] == 0, AMC4R[0] == 1.38,
   GL1[0] == 0, GL2[0] == 0.01, GHR[0] == 670, INS[0] == 79, LEP[0] == 22,
   AINSR[0] == 0.11, AGHSR[0] == 0.6, ALEPR[0] == 0.6};

modelEqs[diet_List] := Join[neuronEqs, hormoneEqs[diet], receptorEqs, initialConditions];

(* ::Section:: *)
(* 5. Simulacion *)
(* Los spikes se detectan con WhenEvent (cruce de -20 mV hacia arriba). *)

popNames = {"AgRP", "POMC", "PVH", "LHA"};

simulate[dietName_String, hours_: 24] := Module[{sol, ev},
   {sol, ev} = Reap[
     NDSolve[Join[modelEqs[diets[dietName]],
       {WhenEvent[A[t] > -20, Sow[t, "AgRP"]], WhenEvent[P[t] > -20, Sow[t, "POMC"]],
        WhenEvent[V[t] > -20, Sow[t, "PVH"]], WhenEvent[L[t] > -20, Sow[t, "LHA"]]}],
      vars, {t, 0, hours}, MaxStepSize -> 1/60, MaxSteps -> Infinity],
     popNames];
   <|"diet" -> dietName, "hours" -> hours, "sol" -> First[sol],
     "spikes" -> AssociationThread[popNames, If[# === {}, {}, First[#]] & /@ ev]|>];

(* Solo hormonas (Fase 1): rapido, sin neuronas *)
simulateHormones[dietName_String, hours_: 72] := First@NDSolve[
    Join[hormoneEqs[diets[dietName]],
     {GL1[0] == 0, GL2[0] == 0.01, GHR[0] == 670, INS[0] == 79, LEP[0] == 22}],
    hormoneVars, {t, 0, hours}, MaxStepSize -> 1/60];

(* Spikes por bin de 30 min (mismo criterio que el paper) *)
spikeCounts[sim_, pop_String] := BinCounts[sim["spikes"][pop], {0, sim["hours"], 0.5}];

(* ::Section:: *)
(* 6. Graficas *)

plotSim[sim_] := Module[{sol = sim["sol"], T = sim["hours"], diet = diets[sim["diet"]], st, x},
   st = {Frame -> True, ImageSize -> Medium, PlotStyle -> Black, LabelStyle -> {10, Black, Bold},
     PlotRangePadding -> {0, Scaled[0.05]}};
   x = Range[0.25, T, 0.5];
   Column[{
     Plot[{food[diet][h], Icirc[h]/Ac}, {h, 0, T}, Filling -> Bottom, PlotRange -> {0, 2.5}, PlotPoints -> 400,
      PlotStyle -> {Gray, Darker[Yellow]}, FrameLabel -> {None, "Food\n(L/h)"}, PlotLabel -> sim["diet"],
      AspectRatio -> 1/5, Frame -> True, ImageSize -> Medium, LabelStyle -> {10, Black, Bold}],
     Plot[Evaluate[GL1[h] /. sol], {h, 0, T}, PlotRange -> {0, 0.8}, FrameLabel -> {None, "GI\:2081\n(L)"}, AspectRatio -> 1/6, Evaluate[Sequence @@ st]],
     Plot[Evaluate[GL2[h] /. sol], {h, 0, T}, PlotRange -> {0, 0.8}, FrameLabel -> {None, "GI\:2082\n(L)"}, AspectRatio -> 1/6, Evaluate[Sequence @@ st]],
     Plot[Evaluate[ghrToPM GHR[h] /. sol], {h, 0, T}, PlotRange -> {0, 300}, FrameLabel -> {None, "Ghrelin\n(pM)"}, AspectRatio -> 1/6, PlotStyle -> Blue, Evaluate[Sequence @@ st]],
     Plot[Evaluate[INS[h] /. sol], {h, 0, T}, PlotRange -> {0, 600}, FrameLabel -> {None, "Insulin\n(pM)"}, AspectRatio -> 1/6, PlotStyle -> Red, Evaluate[Sequence @@ st]],
     Plot[Evaluate[lepToPM LEP[h] /. sol], {h, 0, T}, PlotRange -> {0, 20}, FrameLabel -> {None, "Leptin\n(pM)"}, AspectRatio -> 1/6, PlotStyle -> Darker[Green], Evaluate[Sequence @@ st]],
     Sequence @@ Table[
       ListLinePlot[Transpose[{x, spikeCounts[sim, pop]}], PlotRange -> {0, 80}, Mesh -> Full,
        FrameLabel -> {If[pop == "LHA", "Time (h)", None], pop <> "\n(spikes/30 min)"},
        AspectRatio -> 1/6, Evaluate[Sequence @@ st]], {pop, popNames}]},
    Background -> White, Spacings -> 0.3]];

(* ::Section:: *)
(* 7. Validacion contra la referencia en Python (analysis/reference_model/ref/) *)
(* Criterio: hormonas con error maximo < 1e-3 (normalizado por el maximo de cada hormona);
   spikes: total por poblacion dentro de +-3 % y diferencia maxima por bin <= 5. *)

validate[sim_, refDir_String] := Module[{d = sim["diet"], hor, spk, sol = sim["sol"],
    model, hormErr, rows},
   hor = Rest@Import[FileNameJoin[{refDir, "hormonas_" <> d <> "_24h.csv"}], "CSV"];
   spk = Rest@Import[FileNameJoin[{refDir, "spikes_" <> d <> "_24h.csv"}], "CSV"];
   model = Table[{GL1[r[[1]]], GL2[r[[1]]], ghrToPM GHR[r[[1]]], INS[r[[1]]], lepToPM LEP[r[[1]]]} /. sol, {r, hor}];
   (* error normalizado por el maximo de cada hormona *)
   hormErr = Max[Abs[model - hor[[All, 2 ;;]]]/ConstantArray[Max /@ Transpose[Abs[hor[[All, 2 ;;]]]], Length[hor]]];
   rows = Table[With[{w = spikeCounts[sim, popNames[[k]]][[1 ;; 48]], r = spk[[All, k + 1]]},
      {popNames[[k]], Total[r], Total[w], N@Round[100 (Total[w] - Total[r])/Max[1, Total[r]], 0.1],
       Max[Abs[w - r]],
       If[Abs[Total[w] - Total[r]] <= 0.03 Total[r] && Max[Abs[w - r]] <= 5, "OK", "REVISAR"]}], {k, 4}];
   Column[{
     Row[{"Hormonas, error maximo normalizado: ", ScientificForm[hormErr, 2], "   ",
       If[hormErr < 10^-3, "OK", "REVISAR"]}],
     Grid[Prepend[rows, {"Poblacion", "Spikes Python", "Spikes Wolfram", "Dif %", "Max dif/bin", ""}],
      Frame -> All]}]];

(* ::Section:: *)
(* 8. Tabla de conversion de unidades (registro) *)

unitTable = Grid[Prepend[{
     {"Tiempo t", "s", "h", "t_h = t_s/3600"},
     {"Capacitancia / ecuaciones de voltaje", "C = 1 (t en s)", "Cm = 1/3600", "dV/dt_h = 3600 dV/dt_s"},
     {"Constantes de tiempo neuronales (tauK, tauGLU, tauGABA, tauOX, tauInc, tauDec, MC4R, receptores)", "s", "h", "tau_h = tau_s/3600"},
     {"Tasas hormonales (k12, k2x, khx, kix, beta, beta2, theta, rho, epsilon, R0)", "num[x] = x/3600 1/s", "x 1/h", "se usa x directamente"},
     {"Ingesta (TRE, CD, IF)", "num[2.1] por s", "2.1 L/h", "x directamente"},
     {"Duracion de comidas", "2340, 1676, 1400 s", "0.650, 0.466, 0.389 h", "/3600"},
     {"Frecuencia circadiana wc", "0.000073 rad/s", "0.2628 rad/h", "x 3600"},
     {"Conductancias, potenciales, EC50, concentraciones", "-", "-", "sin cambio"}},
    {"Magnitud", "Original (s)", "Nuevo (h)", "Conversion"}], Frame -> All, Alignment -> Left];

Print["HypothalamusModel.wl cargado. Dietas: ", Keys[diets],
  ". Usa simulate[\"TRE\", 24], plotSim[sim], validate[sim, dir], unitTable."];
