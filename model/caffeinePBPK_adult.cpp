//Caffeine PBPK model for women

$PROB caffeine PBPK women

$PARAM

//Partition coefficients specific to Caffeine calculated with PT
Kpad = 0.1131233  //adipose:plasma PT
Kpbo = 0.4533072  //bone:plasma PT
Kpbr = 0.7779381  //brain:plasma PT
Kpgu = 0.6984072 //gut:plasma PT
Kphe = 0.7512029  //heart:plasma PT
Kpki = 0.7426488  //kidney:plasma PT
Kpli = 0.714217 //liver:plasma PT
Kplu = 0.72736  //lungs:plasma PT
Kpmu = 0.7127541  //muscle:plasma PT
Kpsp = 0.7333799  //spleen:plasma PT
Kpsk = 0.6382312 //skin:plasma PT


//nonpregnant specific parameters
BW = 60 //Kg
Qc = 354 // L/hr

//organ volumes (L)
Vad = 22.5
Vbo = 7.8
Vbr = 1.3
Vgu = 1.03
Vhe = 0.250
Vki = 0.275
Vli = 1.4
Vlu = .42
Vmu = 17.5
Vsk = 2.3
Vsp = .13
Vbl = 3.9
Vpa = .12

//tissue blood flow (l/h)
Qad = 30.09
Qbo = 17.7
Qbr = 42.48
Qgu = 60.18
Qhe = 17.7
Qki = 60.18
Qli = 95.58
Qmu = 42.48
Qsk = 17.7
Qsp = 10.62
Qpa = 3.54


//other parameters specific to Caffeine
f = 0.93 //bioavailability
fg = 1 //gut availability
fa = 1 //fraction available for absorption from dosage form
Ka = 2.18 // First-order absorption rate constant (1/h)
fup = 0.681 //Fraction of drug unbound in plasma
BP = 0.98 //blood:plasma ratio
Cl_hep = 8 //hepatic clearance for single dose (l/h), equal to sum of CYP single dose clearances
Cl_r = 0.0383 //renal clearance for single dose (l/h)


$CMT
LUNG ADIPOSE BONE BRAIN HEART KIDNEY MUSCLE 
SKIN LIVER SPLEEN GUT ART VEN D PANCREAS REST

$MAIN
//additional derivations
double Kppa = (Kpbo + Kpbr + Kpgu + Kphe + Kpki + Kpli + Kplu + Kpmu + Kpsp + Kpsk)/10; //pancreas:plasma 
double Kpre =  (Kpbo + Kpbr + Kpgu + Kphe + Kpki + Kpli + Kplu + Kpmu + Kpsp + Kpsk)/10;//rest of body:plasma
double Vre = BW-(Vad+Vbo+Vbr+Vgu+Vhe+Vki+Vli+Vlu+Vmu+Vsk+Vsp+Vbl+Vpa);//rest of body
double Qha = Qli - Qsp - Qgu - Qpa; //hepatic artery blood flow
double Qre = Qc-(Qad+Qbo+Qbr+Qhe+Qki+Qli+Qmu+Qsk);//rest of the body
double Var = Vbl/3; //arterial
double Vve = Vbl*2/3; //venous
double Qlu = Qc; //lungs
double Qbl = Qc;
double Qar = Qc; //arterial blood
double Qve = Qc; //venous blood

$ODE
//Calculation of tissue drug concentrations (mg/L)
double Cadipose = ADIPOSE/Vad;
double Cbone = BONE/Vbo;
double Cbrain = BRAIN/Vbr; 
double Cgut = GUT/Vgu; 
double Cheart = HEART/Vhe; 
double Ckidney = KIDNEY/Vki;
double Cliver = LIVER/Vli; 
double Cskin = SKIN/Vsk;
double Clung = LUNG/Vlu; 
double Cmuscle = MUSCLE/Vmu;
double Cspleen = SPLEEN/Vsp;
double Carterial = ART/Var;
double Cvenous = VEN/Vve;
double Crest = REST/Vre;
double Cpancreas = PANCREAS/Vpa;

//Free Concentration Calculations
double Cliverfree = Cliver*fup; //mg/L
double Cplasma = Cvenous/BP; //defining plasma concentration
double Ckidneyfree = Ckidney*fup;

//Diferential equations
dxdt_ADIPOSE = Qad*(Carterial - Cadipose/(Kpad/BP)); 
dxdt_BRAIN = Qbr*(Carterial - Cbrain/(Kpbr/BP));
dxdt_HEART = Qhe*(Carterial - Cheart/(Kphe/BP));
dxdt_KIDNEY = Qki*(Carterial - Ckidney/(Kpki/BP))- Ckidney * (Qki * ((fup * Cl_r)/(Qki + fup * Cl_r)));
dxdt_GUT = Qgu*(Carterial - Cgut/(Kpgu/BP)) + Ka*D*fg*fa;
dxdt_LIVER = Qha*Carterial + Qsp*Cspleen/(Kpsp/BP) + Qgu*Cgut/(Kpgu/BP) + Qpa*Cpancreas/(Kppa/BP) - 
  Qli*Cliver/(Kpli/BP) - Cliver * (Qli * ((fup * Cl_hep)/(Qli + fup*Cl_hep))); //need to account for free drug in liver
dxdt_LUNG = Qlu*(Cvenous - Clung/(Kplu/BP));
dxdt_MUSCLE = Qmu*(Carterial - Cmuscle/(Kpmu/BP));
dxdt_SPLEEN = Qsp*(Carterial - Cspleen/(Kpsp/BP));
dxdt_BONE = Qbo*(Carterial - Cbone/(Kpbo/BP));
dxdt_SKIN = Qsk*(Carterial - Cskin/(Kpsk/BP));
dxdt_VEN = Qad*Cadipose/(Kpad/BP) +
  Qbr*Cbrain/(Kpbr/BP) +
  Qhe*Cheart/(Kphe/BP) +
  Qki*Ckidney/(Kpki/BP) +
  Qli*Cliver/(Kpli/BP) +
  Qmu*Cmuscle/(Kpmu/BP) +
  Qbo*Cbone/(Kpbo/BP) +
  Qre*Crest/(Kpre/BP) +
  Qsk*Cskin/(Kpsk/BP) - 
  Qlu*Cvenous;
dxdt_ART = Qlu * (Clung/(Kplu/BP) - Carterial);
dxdt_D = -Ka*fg*fa*D;
dxdt_PANCREAS = Qpa*(Carterial-Cpancreas/(Kppa/BP));
dxdt_REST = Qre*(Carterial-Crest/(Kpre/BP));



$CAPTURE Cplasma Cvenous
  