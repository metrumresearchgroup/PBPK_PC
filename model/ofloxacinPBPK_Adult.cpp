//Ofloxacin PBPK model for a typical adult male


$PROB ofloxacinPBPK_Adult


$PARAM 
//physiological parameters
WEIGHT = 73 //(kg)
C_OUTPUT = 6.5 //(l/min); from ICRP Publication 89

//drug-related parameters
Ka = 0.849 //absorption rate constant(/hr) ; ignored since we don't have the actual Ka
fup = 0.77 //fraction of unbound drug in plasma
BP = 0.92 //blood to plasma ratio; initial estimate

//partition coefficients estimated by Poulin and Theil method (2001)
Kpad = 0.1279002, Kpbo = 0.4461523, Kpbr = 0.7939131, Kpgu = 0.7136217, Kphe = 0.751919, Kpki = 0.765771, Kpli = 0.734187
Kplu = 0.7647848, Kpmu = 0.7453043, Kpsp = 0.7650025,  Kpsk = 0.6562872

//clearance; https://www.karger.com/Article/Abstract/238934
CLrenal = 8.92
CLhepatic = 2.76

$CMT 
D ADIPOSE BRAIN GUT HEART BONE 
KIDNEY LIVER LUNG MUSCLE SPLEEN REST 
ART VEN


$MAIN
//Kp for rest of body
double Kpre =  (Kpbo + Kpbr + Kpgu + Kphe + Kpki + Kpli + Kplu + Kpmu + Kpsp + Kpsk)/10;

//Tissue volumes (L); from ICRP Publication 89
double Vad = 18.2; //adipose
double Vbo = 10.5; //bone; whole skeleton
double Vbr = 1.45; //brain
double VguWall = 1.21; //gut wall
double VguLumen = 0.9; //gut lumen
double Vgu = VguWall + VguLumen; //total gut
double Vhe = 0.33; //heart
double Vki = 0.31; //kidneys
double Vli = 1.8; //liver
double Vlu = 0.5; //lungs
double Vmu = 29; //muscle
double Vsp = 0.15; //spleen
double Vbl = 5.6; //total blood
double Vve = 0.705*Vbl; //venous blood
double Var = 0.295*Vbl; //arterial blood
double Vre = WEIGHT - (Vli+Vki+Vsp+Vhe+Vlu+Vbo+Vbr+Vmu+Vad+Vgu+Vbl); //volume of rest of the body compartment

//fractions of tissue blood flows; from ICRP Publication 89
double FQad = 0.05;
double FQbo = 0.05;  
double FQbr = 0.12;
double FQgu = 0.16;
double FQhe = 0.04;
double FQki = 0.19;
double FQli = 0.255;
double FQmu = 0.17;
double FQsp = 0.03;

//computing the blood flows for each tissue
double CO = C_OUTPUT*60; //scaled cardiac output (L/hr)
double Qad = FQad*CO;
double Qbo = FQbo*CO;
double Qbr = FQbr*CO;
double Qgu = FQgu*CO;
double Qhe = FQhe*CO;
double Qki = FQki*CO;
double Qli = FQli*CO;
double Qmu = FQmu*CO;
double Qsp = FQsp*CO;
double Qha = Qli - (Qgu + Qsp); //hepatic artery 
double Qtot = Qli + Qki + Qbo + Qhe + Qmu + Qad + Qbr;
double Qre = CO - Qtot;
double Qlu = CO;


$ODE
//Calculation of tissue drug concentrations (mg/L)
double Cadipose = ADIPOSE/Vad;
double Cbone = BONE/Vbo;
double Cbrain = BRAIN/Vbr; 
double Cgut = GUT/VguWall; 
double Cheart = HEART/Vhe; 
double Ckidney = KIDNEY/Vki;
double Cliver = LIVER/Vli; 
double Clung = LUNG/Vlu; 
double Cmuscle = MUSCLE/Vmu;
double Cspleen = SPLEEN/Vsp;
double Crest = REST/Vre;
double Carterial = ART/Var;
double Cvenous = VEN/Vve;
double Cgut_D = D/Vgu;


//ODEs
dxdt_D = -Ka*D;
dxdt_ADIPOSE = Qad*(Carterial - Cadipose/(Kpad/BP)); 
dxdt_BRAIN = Qbr*(Carterial - Cbrain/(Kpbr/BP));
dxdt_HEART = Qhe*(Carterial - Cheart/(Kphe/BP));
dxdt_KIDNEY = Qki*(Carterial - Ckidney/(Kpki/BP)) - CLrenal*(Ckidney/(Kpki/BP));
dxdt_GUT = Ka*D + Qgu*(Carterial - Cgut/(Kpgu/BP)); 
dxdt_LIVER = Qgu*(Cgut/(Kpgu/BP)) + Qsp*(Cspleen/(Kpsp/BP)) + Qha*(Carterial) - Qli*(Cliver/(Kpli/BP)) - 
  CLhepatic*(Cliver/(Kpli/BP)); 
dxdt_LUNG = Qlu*(Cvenous - Clung/(Kplu/BP));
dxdt_MUSCLE = Qmu*(Carterial - Cmuscle/(Kpmu/BP));
dxdt_SPLEEN = Qsp*(Carterial - Cspleen/(Kpsp/BP));
dxdt_BONE = Qbo*(Carterial - Cbone/(Kpbo/BP));
dxdt_REST = Qre*(Carterial - Crest/(Kpre/BP));
dxdt_VEN = Qad*(Cadipose/(Kpad/BP)) + Qbr*(Cbrain/(Kpbr/BP)) +
  Qhe*(Cheart/(Kphe/BP)) + Qki*(Ckidney/(Kpki/BP)) + Qli*(Cliver/(Kpli/BP)) + 
  Qmu*(Cmuscle/(Kpmu/BP)) + Qbo*(Cbone/(Kpbo/BP)) + Qre*(Crest/(Kpre/BP)) - Qlu*Cvenous;
dxdt_ART = Qlu*(Clung/(Kplu/BP) - Carterial);

$TABLE
double Cplasma = Cvenous/BP;

$CAPTURE Cvenous Cplasma
  
  
  