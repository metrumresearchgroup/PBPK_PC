//Alfentanil PBPK model for representative adult male/female mix
//http://anesthesiology.pubs.asahq.org/article.aspx?articleid=1948182


$PROB alfentanilPBPK_Adult


$PARAM 
//physiological parameters
WEIGHT = 70 //(kg) average weight as reported in http://anesthesiology.pubs.asahq.org/article.aspx?articleid=1948182
C_OUTPUT = 6.5 //(l/min); from ICRP Publication 89

//drug-related parameters
Ka = 0.849 //absorption rate constant(/hr) 
fup = 0.11 //fraction of unbound drug in plasma
BP = 0.63 //blood to plasma ratio; initial estimate

//partition coefficients estimated by Poulin and Theil method (2001)
Kpad = 0.9381411, Kpbo = 4.48555, Kpbr = 3.901444, Kpgu = 3.38532, Kphe = 5.743885, Kpki = 2.666754, Kpli = 2.869263
Kplu = 0.6260866, Kpmu = 1.167515, Kpsp = 1.452498, Kpsk=2.535692

//clearance; https://www.karger.com/Article/Abstract/238934
CLrenal = 0
CLhepatic = 0.36*60*1.8  //https://bjanaesthesia.org/article/S0007-0912(17)42730-9/pdf

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

double fub = fup/BP;  //free fraction in blood

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
double Cplasma = (Cvenous/BP)*1000;
double Ca = Carterial/BP;

$CAPTURE Cvenous Cplasma Ca
  
  
  