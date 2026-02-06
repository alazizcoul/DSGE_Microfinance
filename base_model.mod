// Modèle Microfinance DSGE 2026 .
// Ce modèles est une modification de EA_GNSS10 model
// THEME --  RISQUE DE CREDIT ET SOLVABILITE EN MICROFINANCE : DETERMINATION D'UN SEUIL OPTIMAL DU PAR90 et AJUSTEMENT DE LA 
//                                                             NORME DE CAPITALISATION
// Auteur : AL AZIZ N'GOLO COULIBALY
//          Ingénieur Statisticien Economiste (ENSEA)

// Version: 02 01 2026
// BBS has all items dated (t) - Microfinance capital in real terms is defined as: K_m(t)/p(t)
// Banking profits are defined in the model code at time (t) - All equations in exp form
// This is the full version (monop. competitive banking sector & Microfinance capital à la Gerali with sticky rates)  


var 
c_p       // 1  PATIENT   HHs
d_p       // 3  PATIENT   HHs
l_p       // 4  PATIENT   HHs
lam_p     // 5  PATIENT   HHs
t_p       //    PATIENT   HHs
J_R       // 6  PATIENT   HHs
j_m       // 7  PATIENT   HHs
pie_wp    // 8  PATIENT   HHs  
c_i       // 9  IMPATIENT HHs
b_i       // 11 IMPATIENT HHs 
l_i       // 12 IMPATIENT HHs
lam_i     // 13 IMPATIENT HHs
//s_i       // 14 IMPATIENT HHs
pie_wi    // 15 IMPATIENT HHs
I         // 16 CAPITAL PRODUCERS
q_k       // 17 CAPITAL PRODUCERS
c_e       // 18 ENTREPRENEURS
k_e       // 19 ENTREPRENEURS
l_pd      // 20 ENTREPRENEURS
l_id      // 21 ENTREPRENEURS
b_ee      // 22 ENTREPRENEURS
y_e       // 23 ENTREPRENEURS
lam_e     // 24 ENTREPRENEURS
s_e       // 25 ENTREPRENEURS x
u         // 26 ENTREPRENEURS Capital utilization rate
d_b       // 27 IMFs
b_h       // 28 IMFs
b_e       // 29 IMFs
r_d       // 30 IMFs
r_bh      // 31 IMFs
r_be      // 32 IMFs
R_b       // 33 IMFs
K_m       // 34 IMFs
pie       // 35 RETAILERS
x         // 36 RETAILERS
C         // 37 AGGREGATION & EQUILIBRIUM
Y         // 38 AGGREGATION & EQUILIBRIUM
D         // 39 AGGREGATION & EQUILIBRIUM
BE        // 40 AGGREGATION & EQUILIBRIUM
BH        // 41 AGGREGATION & EQUILIBRIUM
B         // 42 AGGREGATION & EQUILIBRIUM
w_p       // 43 AGGREGATION & EQUILIBRIUM
w_i       // 44 AGGREGATION & EQUILIBRIUM
J_m       // 45 AGGREGATION & EQUILIBRIUM
K         // 47 AGGREGATION & EQUILIBRIUM
PIW       // 48 AGGREGATION & EQUILIBRIUM
r_ib      // 49 MONETARY POLICY
r_k       // 50 CAPITAL RENTAL RATE
ee_z      // 51 EXOGENOUS PROCESSES
A_e       // 52 EXOGENOUS PROCESSES
mk_d      // 54 EXOGENOUS PROCESSES
mk_be     // 55 EXOGENOUS PROCESSES
mk_bh     // 56 EXOGENOUS PROCESSES
ee_qk     // 57 EXOGENOUS PROCESSES
m_i       // 58 EXOGENOUS PROCESSES (IMPATIENT LTV)
m_e       // 59 EXOGENOUS PROCESSES (ENTREPRENEURS LTV)
eps_y     // 60 EXOGENOUS PROCESSES
eps_l     // 61 EXOGENOUS PROCESSES
eps_K_m   // 62 EXOGENOUS PROCESSES
Y1        // 63 output a prezzi di ss
rr_e      // 64 Entrep. Real Rate
RatioCap  // 65 auxiliary variable
bm        // 66 IMFs intermediation margins
spr_b     // 67 average bank spread (active-passive) (67-3 = 64)
zeta_e
zeta_i
vi;

//**************************************************************************
// Replication Variables                                                  //*
 //ROA interestPol interestH interestF inflation
 //loansH loansF output consumption investment deposits interestDep IMFcapital;
//**************************************************************************

varexo  e_zeta_e e_zeta_i  e_A_e e_eps_K_m e_l e_me e_mi e_mk_be e_mk_bh e_mk_d e_r_ib e_qk e_y e_z;   //14  
  
parameters  
            beta_p  phi beta_i m_i_ss beta_e m_e_ss alpha eksi_1 eksi_2     // HOUSEHOLDS & ENTREPRENEURS , j
            a_i a_p a_e gamma_p gamma_i gamma_e    ni                       // HOUSEHOLDS & ENTREPRENEURS , h
            eps_l_ss kappa_w                                                // HOUSEHOLDS (labor params)
            eps_d eps_bh eps_be                                             // IMFs 
            mk_d_ss mk_bh_ss mk_be_ss r_be_ss  r_bh_ss r_k_ss               // IMFs (SS)
            gamma_b  deltakm   kappa_km                                     // IMFs 
            eps_y_ss kappa_p ind_p ind_w                                    // RETAILERS
            kappa_i kappa_d kappa_be kappa_bh deltak                        // OTHERS
            ind_d ind_be ind_bh                                             // OTHERS
            rho_ib phi_pie phi_y                                            // POLICY
            piss  r_ib_ss                                                   // STEADY STATE
            rho_ee_z rho_A_e  rho_mi rho_me rho_eps_y                       // SHOCKS
            rho_mk_d rho_mk_be rho_mk_bh rho_ee_qk rho_eps_l rho_eps_K_m    // SHOCKS
            rho_vi vi_ss chi_nu_zeta 
            zeta_bar zeta_e_ss zeta_i_ss rho_zeta_e rho_zeta_i 
            chi_zeta_y chi_nu_y vi_bar  omega_m Z_i  Z_p
            ;
             
% *********************			
% CALIBRATED PARAMETERS
% *********************

beta_p       = 0.9892 ;                                                   % discount factor patient households
beta_i       = 0.965 ;                                                    % discount factor impatient households     
beta_e       = beta_i;                                                    % discount factor entrepreneurs
phi          = 2;                                                         % inverse Frisch elasticity of labor supply
m_i_ss       = 0.35  ;                                                    % loan-to-value ratio impatient households
m_e_ss       = 0.35 ;                                                     % loan-to-value ratio entrepreneurs
alpha        = 1/3 ;                                                      % capital share in the production function
eps_d        = -1.3; % -1.46025;                                          % elast. of subst. of deposits 
eps_bh       = 2.9; %2.932806;
eps_be       = 2.9; %2.932806; 
 
mk_d_ss      = eps_d   / (eps_d  - 1) ;                                    % steady state markdown on D (ok if eps_d<0; if eps_d>0 it should be eps_d/(eps_d+1) )
mk_bh_ss     = eps_bh  / (eps_bh - 1) ;                                    % steady state markup on loans to I
mk_be_ss     = eps_be  / (eps_be - 1) ;                                    % steady state markup on loans to E
eps_y_ss     = 4.5; %6;                                                          % 
eps_l_ss     = 4 ; %3.5; %4;                                                          % 
gamma_p      = 1;                                                          % shares of patient households
gamma_i      = 1; //1/3;                                                   % shares of impatient households
ni           = 0.4;                                                        % wage share of patient households
gamma_b      = 1; //0.10;												   % shares of bankers
gamma_e      = 1; //1 - gamma_p - gamma_i;								   % shares of entrepreneurs
deltak       = 0.035; %0.025                                                      % depreciation rate for physical capital
piss         = 1;                                                          % steady state gross inflation rate

zeta_e_ss    = 0.07 ;
zeta_i_ss    = 0.07 ;
r_ib_ss      = (1/beta_p - 1) /mk_d_ss ;                       % steady state gross nominal interest rate 
r_be_ss      = (r_ib_ss + zeta_e_ss)*eps_be/((eps_be-1)*(1 - zeta_e_ss)) ;								   % steady state interest rate on loans to E
r_bh_ss      = (r_ib_ss + zeta_i_ss)*eps_bh/((eps_bh-1)*(1 - zeta_i_ss)) ;								   % steady state interest rate on loans to H

rho_ib      =	0.513821 ;   
phi_pie     =   2;
phi_y       =	-0.16086462 ; 

% POLITIQUE MACRO PRUDENTIELLE
% =============Paremètres d'évolution du risque de crédit(RC)=========
rho_zeta_e    = 0.838502 ;   % persistence du RC
rho_zeta_i    = 0.838502 ;   % persistence du RC

chi_zeta_y  = -0.087765793; % effet de la croissance économique sur le RC

% =============Règle d'évolution de la norme de capitalisation========
chi_nu_y    = 0 ;  %0.002574401 ; %limite = 0.2; %est 0.002574401; % Sensibilité de l'Autorité Macroprudentielle par rapport niveau des prêts/PIB
rho_vi      = 0.96264 ;  %0;      % persistence de la norme de solvabilité (ici 0, pas de dynamisme dans la norme)
zeta_bar    = 0.03; %-0.07;     % Niveau cible du risque de crédit pour l'AM
vi_bar      = 0.15 ;         % Niveau cible de la norme de cap

chi_nu_zeta = 0 ;          % Sensibilité de l'Autorité Macroprudentielle par rapport au risque de crédit


vi_ss        = vi_bar + chi_nu_zeta*( 0.5 * zeta_e_ss + 0.5 * zeta_i_ss - zeta_bar) ;         %      Norme  de capitalisation  à l'état stationnaire
%=========================================================================================================================================================================================================


// Sous hypothèse zeta_e_ss = zeta_i_ss et eps_be = eps_bh
// Or on a: deltakm  /omega_m    = (1-zeta_i_ss)*(1 + r_bh_ss) b_i /K_m + (1-zeta_e_ss)*(1 + r_be_ss) b_e /K_m - r_d * D/K_m
%=========================================================================================================================================================================================================
omega_m      = 0.6 ;
deltakm    = omega_m * ( (1-zeta_i_ss)*(1 + r_bh_ss)/vi_ss - 1/beta_p *(1 - vi_ss)/vi_ss ) ;

%=========================================================================================================================================================================================================
r_k_ss       = 1/beta_e -(1-deltak)-m_e_ss*(1-deltak)/beta_e*(1/(1+r_be_ss)-beta_e * (1-zeta_e_ss));                       % steady state rental rate of capital
eksi_1       = r_k_ss;   
eksi_2       = 0.1*r_k_ss; 
%=========================================================================================================================================================================================================

Z_i = 1 + ( 1 - (1 - zeta_i_ss) * (1 + r_bh_ss) ) * m_i_ss /(1 + r_bh_ss) ;
Z_p =  1 + 1/ni * ( ( 1/beta_p - 1 ) * (1 - vi_ss) + ( 1/omega_m - 1 ) * vi_ss * deltakm )*( m_e_ss *(1-deltak)* alpha/( (1 + r_be_ss) * r_k_ss * (1 - alpha)) + m_i_ss * (1 - ni)/(1 + r_bh_ss)) 
          + 1/(ni * (1 - alpha)*gamma_e * (eps_y_ss - 1)); 


%Z_p = (eps_y_ss-1)/eps_y_ss * ( 1 + 1/ni * ( ( 1/beta_p - 1 ) * (1 - vi_ss) + ( 1/omega_m - 1 ) * vi_ss * deltakm )*( m_e_ss *(1-deltak)* alpha/( (1 + r_be_ss) * r_k_ss * (1 - alpha)) + m_i_ss * (1 - ni)/(1 + r_bh_ss)) 
%      + (1/eps_y_ss) * ( 1/(ni * (1 - alpha) + ( 1 - (1 - zeta_i_ss) * (1 + r_bh_ss) ) * m_i_ss * (1 - ni) /((1 + r_bh_ss) * ni) + ( 1 - (1 - zeta_e_ss) * (1 + r_be_ss) ) * m_e_ss * (1 - deltak)*alpha /((1 + r_be_ss) * r_k_ss * ni * (1-alpha))  ))) ;
%=========================================================================================================================================================================================================
%BBB =  m_e_ss * (1 - deltak)*alpha /((1 + r_be_ss) * r_k_ss * ni * (1-alpha))

ind_d        = 0;                   % indexation deposit rates
ind_be       = 0;                   % indexation rates on loans to firms
ind_bh       = 0;                   % indexation rates on loans to households
% *****************************************************************
% LOADING MEDIAN OF POSTERIOR: USES EXTRACT_MEDIAN_FROM_POSTERIOR.m (dummy way)
% *****************************************************************
load median_values.txt;
coeffs = median_values;

rho_ee_z	=	coeffs(1);  % 0.385953438168178	;
rho_A_e     =	coeffs(2);  % 0.93816527333294	;
rho_me      =	coeffs(4);  % 0.90129485520182	;
rho_mi      =	coeffs(5);  % 0.922378382753078	;
rho_mk_d	=	coeffs(6);  % 0.892731352899547	;
rho_mk_bh	=	coeffs(7);  % 0.851229673864555	;
rho_mk_be	=	coeffs(8);  % 0.873901213475799	;
rho_ee_qk	=	coeffs(9);  % 0.571692383714171	;
rho_eps_y	=	coeffs(10); % 0.294182239567384	;
rho_eps_l	=	coeffs(11); % 0.596186440884132	;
rho_eps_K_m	=	coeffs(12); % 0.813022758608552	;
kappa_p     =	coeffs(13); % 33.7705265016395	;
kappa_w     =	coeffs(14); % 107.352040072465	;
kappa_i     =	coeffs(15); % 10.0305562248008	;
kappa_d     =	coeffs(16); % 2.77537377104213	;
kappa_be	=	coeffs(17); % 7.98005959044637	;
kappa_bh	=	coeffs(18); % 9.04426718749482	;
kappa_km	=	coeffs(19); % 8.91481958034669	;
%===phi_pie     =   2;  % coeffs(20); % 2.00384780180824  ;
%===rho_ib      =	0.508944;   % 0.750481873084311	;
%===phi_y       =	-0.157796667; % 0.303247771697294	;
ind_p       =	coeffs(23); % 0.158112794106546	;
ind_w       =	coeffs(24); % 0.300197804017489	;
a_i	        =	coeffs(25); % 0.867003766306404	;
a_e         =   coeffs(25) ;% 0.0     ;   % degree of habit formation: entrepreneurs
a_p         =   coeffs(25); %0.0     ;   % degree of habit formation: patient households

//%------------------------------------------------------------
//% Model equations
//%------------------------------------------------------------

model;

//**************************************************************************
// Definition of Replication Variables in Terms of Original Model Variables //*
/*
interestPol   = 400*exp(r_ib); 
interestH     = 400*exp(r_bh);
interestF     = 400*exp(r_be);                                         
inflation     = pie*100;    
loansH        = BH*100;
loansF        = BE*100;
output        = Y1*100; 
consumption   = C*100;
investment    = I*100;   
deposits      = D*100;
interestDep   = 400*(exp(r_d));
IMFcapital    = 100*K_m;                 
ROA           = exp(J_m)/exp(B) ;       
*/
//**************************************************************************

// Model code:

////***********   1) PATIENT HHs ********************************************************6

(1-a_p)*exp(ee_z)*(exp(c_p) - a_p*exp(c_p(-1)))^(-1) = exp(lam_p); // CON rescaling  (1) where a_p = a_e = a_i

exp(lam_p)  = beta_p * exp(lam_p(+1)) * (1+exp(r_d)) / exp(pie(+1)); // (2)

(1 - exp(eps_l)) * exp(l_p) + exp(l_p) ^(1+phi) / exp(w_p) * exp(eps_l)/exp(lam_p) 
                                         - kappa_w *( exp(pie_wp)     - exp(pie(-1)) ^ ind_w * piss ^ (1-ind_w) ) * exp(pie_wp)
   +  beta_p * exp(lam_p(+1))/exp(lam_p) * kappa_w *( exp(pie_wp(+1)) - exp(pie)     ^ ind_w * piss ^ (1-ind_w) ) * exp(pie_wp(+1)) ^2 / exp(pie(+1)) = 0 ; // Unions, labor supply (3)

exp(pie_wp) = exp(w_p) / exp(w_p(-1)) * exp(pie); // definition of wage inflation (4)

exp(c_p)  + exp(d_p)  = exp(w_p) * exp(l_p)
   + (1+exp(r_d(-1)))*exp(d_p(-1))/exp(pie) + exp(t_p)  ;  //patient household budget constraint (5), exp(J_R)/gamma_p = t_p

exp(t_p) = exp(J_R)/gamma_p +  (1-omega_m)* exp(J_m) ; // (6)

////***********   2) IMPATIENT HHs ********************************************************7

(1-a_i)*exp(ee_z)*(exp(c_i) - a_i*exp(c_i(-1)))^(-1)  = exp(lam_i); // CON rescaling (7)

(1 - exp(eps_l)) * exp(l_i) + exp(l_i) ^(1+phi) / exp(w_i) * exp(eps_l)/exp(lam_i) 
                                         - kappa_w *( exp(pie_wi)     - exp(pie(-1))^ind_w * piss ^ (1-ind_w) ) * exp(pie_wi)
   +  beta_i * exp(lam_i(+1))/exp(lam_i) * kappa_w *( exp(pie_wi(+1)) - exp(pie)    ^ind_w * piss ^ (1-ind_w) ) * exp(pie_wi(+1)) ^2 / exp(pie(+1)) = 0; // (8)

exp(pie_wi) = exp(w_i) / exp(w_i(-1)) * exp(pie); // (9)

exp(c_i)  + (1-exp(zeta_i))*(1+exp(r_bh(-1)))*exp(b_i(-1))/exp(pie) =  exp(w_i) * exp(l_i) + exp(b_i)  ;  // (10) budget constraint

(1+exp(r_bh)) * exp(b_i) = exp(m_i) * exp(w_i(+1)) *exp(l_i) * exp(pie_wi(+1));     // (11) borrowing constraint impatient household 

////***********  3) CAPITAL PRODUCERS *****************************************************

exp(K) = (1-deltak) * exp(K(-1)) + ( 1 - kappa_i/2 * (exp(I)*exp(ee_qk)/exp(I(-1)) - 1)^2 ) * exp(I) ; // (12)  

1 = exp(q_k) * ( 1 -  kappa_i/2 * (exp(I)*exp(ee_qk)/exp(I(-1)) - 1)^2  - kappa_i * (exp(I)*exp(ee_qk)/exp(I(-1)) - 1) * exp(I)*exp(ee_qk)/exp(I(-1)) ) 
  + beta_e * exp(lam_e(+1)) / exp(lam_e) * exp(q_k(+1)) *   kappa_i * (exp(I(+1))*exp(ee_qk(+1))/exp(I) - 1) * exp(ee_qk(+1)) * (exp(I(+1))/exp(I))^2 ; // real price of capital (13)


////************  4) ENTREPRENEURS *********************************************************

(1-a_e)*(exp(c_e) - a_e*exp(c_e(-1)))^(-1) = exp(lam_e);         // CON rescaling  (14) FOC consumption
//        (exp(c_e) - a_e*exp(c_e(-1)))^(-1) = exp(lam_e);         // SENZA rescaling 

exp(s_e)  * exp(m_e) * exp(q_k(+1)) * exp(pie(+1)) * (1-deltak) 
     + beta_e * exp(lam_e(+1)) * ( exp(q_k(+1))*(1-deltak) + exp(r_k(+1))*exp(u(+1))
     - ( eksi_1*(exp(u(+1))-1)+eksi_2/2*( (exp(u(+1))-1)^2 ) ) )   = exp(lam_e) * exp(q_k) ;  // (15) FOC capital

exp(w_p) =    ni  * (1-alpha) * exp(y_e) / ( exp(l_pd) * exp(x) ); // (16) FOC labor patient households
exp(w_i) = (1-ni) * (1-alpha) * exp(y_e) / ( exp(l_id) * exp(x) ); // (17) FOC labor impatient households

exp(lam_e) - exp(s_e)  * (1+exp(r_be)) = beta_e * exp(lam_e(+1)) * (1+exp(r_be))*(1-exp(zeta_e(+1))) / exp(pie(+1));  // (20) FOC credit demand
exp(r_k)  = eksi_1 + eksi_2 * (exp(u)-1); // (18) FOC utilization rate

exp(c_e) + (1-exp(zeta_e))*((1+exp(r_be(-1))) * exp(b_ee(-1)) / exp(pie) ) +  (exp(w_p)*exp(l_pd) + exp(w_i)*exp(l_id)) + exp(q_k) * exp(k_e) 
   + ( eksi_1*(exp(u)-1)+eksi_2/2*(exp(u)-1)^2 ) * exp(k_e(-1)) = 
    exp(y_e) / exp(x) + exp(b_ee) + exp(q_k) * (1-deltak) * exp(k_e(-1))  ;   // budget constraint entrepreneurs (19)

exp(y_e) = exp(A_e) * (exp(u)*exp(k_e(-1)))^(alpha) * ( exp(l_pd)^ni * exp(l_id)^(1-ni) ) ^ (1-alpha); // production technology (20)

(1+exp(r_be)) * exp(b_ee) = exp(m_e) * exp(q_k(+1))  *exp(pie(+1)) * exp(k_e) * (1-deltak); // borrowing constraint entrepreneurs (21)

exp(r_k) = alpha * exp(A_e) * exp(u)^(alpha-1) * exp(k_e(-1))^(alpha-1) * ( exp(l_pd)^ni * exp(l_id)^(1-ni) ) ^ (1-alpha) /exp(x);  // definition (22)

////*************  5)IMFs ****************************************************************

exp(R_b) = exp(r_ib) - kappa_km * ( exp(K_m) / exp(B) - exp(vi) ) * (exp(K_m)/exp(B))^2   ; // (23) FOC wholesale branch, where r_ib is policy rate

exp(K_m) * exp(pie) = (1-deltakm) * exp(K_m(-1)) / exp(eps_K_m) + omega_m* exp(j_m(-1)) ; // (24) accumulation of bank capital

gamma_b * exp(d_b)  = gamma_p * exp(d_p) ;    //(25)                                  
gamma_b * exp(b_h)  = gamma_i* exp(b_i) ;   //(26)
gamma_b * exp(b_e)  = gamma_e * exp(b_ee);  //(27)

exp(b_h) + exp(b_e)  =  exp(d_b) + exp(K_m) ; //(28)

/// PRICING in terms of MK ///

// (28) FOC deposit branch:
- 1 + exp(mk_d)/(exp(mk_d)-1)  - exp(mk_d)/(exp(mk_d)-1)  * exp(r_ib)/exp(r_d)  - kappa_d  * ( exp(r_d)/exp(r_d(-1)) - ( exp(r_d(-1)) / exp(r_d(-2)) )^ind_d  )  * exp(r_d)/exp(r_d(-1)) 
  + beta_p * ( exp(lam_p(+1))/exp(lam_p) ) * kappa_d  * ( exp(r_d(+1))/exp(r_d) - ( exp(r_d)/exp(r_d(-1)))^ind_d )   * ( (exp(r_d(+1))/exp(r_d))^2 )   * (exp(d_b(+1))/exp(d_b)) = 0;//(29) 
  
// (29) FOC loan branch, entrepreneurs:
(+ 1 - exp(mk_be)/(exp(mk_be)-1) )*(1 - exp(zeta_e(+1)))  +  exp(mk_be)/(exp(mk_be)-1)  * (exp(R_b) + exp(zeta_e(+1)))/exp(r_be) - kappa_be * (exp(r_be)/exp(r_be(-1)) - ( exp(r_be(-1)) / exp(r_be(-2)) )^ind_be ) * exp(r_be)/exp(r_be(-1)) 
  + beta_p * ( exp(lam_p(+1))/exp(lam_p) ) * kappa_be * ( exp(r_be(+1))/exp(r_be) - ( exp(r_be)/exp(r_be(-1)))^ind_be ) * ( (exp(r_be(+1))/exp(r_be))^2 ) * (exp(b_e(+1))/exp(b_e)) = 0;//   (30)
  
// (30) FOC loan brach, impatient households
(+ 1 - exp(mk_bh)/(exp(mk_bh)-1) )*(1 - exp(zeta_i(+1)))  +  exp(mk_bh)/(exp(mk_bh)-1)  * (exp(R_b) + exp(zeta_i(+1)))/exp(r_bh) - kappa_bh * (exp(r_bh)/exp(r_bh(-1)) - ( exp(r_bh(-1)) / exp(r_bh(-2)))^ind_bh ) * exp(r_bh)/exp(r_bh(-1)) 
  + beta_p * ( exp(lam_p(+1))/exp(lam_p) ) * kappa_bh * ( exp(r_bh(+1))/exp(r_bh) - ( exp(r_bh)/exp(r_bh(-1)))^ind_bh ) * ( (exp(r_bh(+1))/exp(r_bh))^2 ) * (exp(b_h(+1))/exp(b_h)) = 0;// (31)

// (31) overall microfinance profits:
exp(j_m) = +((1 -exp(zeta_i(+1)))* exp(r_bh) - exp(zeta_i(+1)))*  exp(b_h)
           +((1 -exp(zeta_e(+1)))* exp(r_be) -exp(zeta_e(+1))) *  exp(b_e) 
           - exp(r_d)   *  exp(d_b)           
           - kappa_d/2  * ( (exp(r_d)/exp(r_d(-1))-1)^2)   * exp(r_d) *exp(d_b) 
           - kappa_be/2 * ( (exp(r_be)/exp(r_be(-1))-1)^2) * exp(r_be)*exp(b_e) 
           - kappa_bh/2 * ( (exp(r_bh)/exp(r_bh(-1))-1)^2) * exp(r_bh)*exp(b_h)
           - kappa_km/2 * ( (exp(K_m) / exp(B)  - exp(vi)) ^2) * exp(K_m); //  32
%%========================================================================================================================================================================================================================
 
%% Calcul de vi selon les normes dynamique ou statique
     % Norme Statique
   % exp(vi) = vi_ss;  % Situation statique, vi est fixé à sa valeur stationnaire
    
     % Norme Dynamique
     exp(vi) = exp(vi(-1)) * rho_vi + vi_ss * (1 - rho_vi) +  (0.5 * exp(zeta_e) + 0.5 * exp(zeta_i) ) * (1 - rho_vi) * chi_nu_zeta +  (exp(B) /exp(Y)  -exp(steady_state(B))/ exp(steady_state(Y)) ) * (1 - rho_vi) * chi_nu_y;
 %========================================================================================================================================================================================================================

////***********  6)RETAILERS **************************************************************

exp(J_R)  = exp(Y)*(1 - (1/exp(x))  - (kappa_p/2) * (exp(pie) - ( exp(pie(-1)) ^ ind_p * piss ^ (1-ind_p) ))^2 ) ; // (32) aggregate retail profits

// (33) New Keynesian Phillips curve:
1 - exp(eps_y) + exp(eps_y) / exp(x) -    kappa_p * (exp(pie)     - ( exp(pie(-1)) ^ ind_p * piss ^ (1-ind_p) )) * exp(pie) 
    + beta_p*(exp(lam_p(+1))/exp(lam_p))* kappa_p * (exp(pie(+1)) - ( exp(pie)     ^ ind_p * piss ^ (1-ind_p) )) * exp(pie(+1)) * (exp(Y(+1))/exp(Y)) = 0;  // 34

////************  7) AGGREGATION & EQUILIBRIUM  ************************************************

exp(C)              = gamma_p * exp(c_p) + gamma_i * exp(c_i) + gamma_e * exp(c_e);
exp(BH)             = gamma_b * exp(b_h);
exp(BE)             = gamma_b * exp(b_e);
exp(B)              = (exp(BH) + exp(BE)); //
exp(D)              = gamma_p * exp(d_p) ; // oppure: (gamma_b * exp(d_b))
exp(J_m)            = gamma_b * exp(j_m);  // 
gamma_e * exp(l_pd) = gamma_p * exp(l_p);
gamma_e * exp(l_id) = gamma_i * exp(l_i);
exp(K)              = gamma_e * exp(k_e); //
% exp(Y1)             = exp(C) +    1     * (exp(K)-(1-deltak)*exp(K(-1))) + deltakm * exp(K_m(-1)); //   47
exp(Y1)             = exp(C) +    1     * (exp(K)-(1-deltak)*exp(K(-1))) ; //   47
exp(Y)              = gamma_e * exp(y_e); //
//exp(Y)            = exp(C) + exp(q_k) * (exp(K)-(1-deltak)*exp(K(-1))) + deltakm * exp(K_m(-1))
//                      + (eksi_1*(exp(u)-1) + eksi_2/2*((exp(u)-1)^2))
//                      + kappa_p/2  * (  exp(pie) - ( exp(pie(-1)) ^ ind_p * piss ^ (1-ind_p) ))^2 * exp(Y)
//                      + kappa_d/2  * ( (exp(r_d(-1))/exp(r_d(-2))-1)^2)   * exp(r_d(-1)) *exp(d_b(-1)) 
//                      + kappa_be/2 * ( (exp(r_be(-1))/exp(r_be(-2))-1)^2) * exp(r_be(-1))*exp(b_e(-1)) 
//                      + kappa_bh/2 * ( (exp(r_bh(-1))/exp(r_bh(-2))-1)^2) * exp(r_bh(-1))*exp(b_h(-1))
//                         ;  //
exp(PIW)            = ( exp(w_p) + exp(w_i) )  / ( exp(w_p(-1)) + exp(w_i(-1)) ) * exp(pie);


////***********  8) TAYLOR RULE & PROFITS CB *****************************************************                                                  

(1+exp(r_ib)) = (1+r_ib_ss)^(1 - rho_ib ) * (1+exp(r_ib(-1)))^rho_ib * (( exp(pie) / piss ) ^phi_pie *  
             (exp(Y1)/exp(Y1(-1)))^phi_y  ) ^ ( 1 - rho_ib ) * (1+e_r_ib) ;// (34) monetary policy

 
////***********  9) EXOGENOUS PROCESSES ****************************************************12

exp(ee_z)     = (1 - rho_ee_z) *    1          + rho_ee_z   * exp(ee_z(-1))    + e_z;
exp(A_e)      = (1 - rho_A_e)  *    1          + rho_A_e    * exp(A_e(-1))     + e_A_e;
exp(m_i)      = (1-rho_mi)     *  m_i_ss       + rho_mi     * exp(m_i(-1))     + e_mi;
exp(m_e)      = (1-rho_me)     *  m_e_ss       + rho_me     * exp(m_e(-1))     + e_me;
exp(mk_d)     = (1-rho_mk_d)   * mk_d_ss       + rho_mk_d   * exp(mk_d(-1))    + e_mk_d;
exp(mk_be)    = (1-rho_mk_be)  * mk_be_ss      + rho_mk_be  * exp(mk_be(-1))   + e_mk_be;
exp(mk_bh)    = (1-rho_mk_bh)  * mk_bh_ss      + rho_mk_bh  * exp(mk_bh(-1))   + e_mk_bh;
exp(ee_qk)    =  1-rho_ee_qk   *    1          + rho_ee_qk  * exp(ee_qk(-1))   + e_qk;
exp(eps_y)    = (1-rho_eps_y)  * eps_y_ss      + rho_eps_y  * exp(eps_y(-1))   + e_y;
exp(eps_l)    = (1-rho_eps_l)  * eps_l_ss      + rho_eps_l  * exp(eps_l(-1))   + e_l;
exp(eps_K_m)  = (1-rho_eps_K_m)*    1          + rho_eps_K_m* exp(eps_K_m(-1)) + e_eps_K_m;

%exp(zeta_e) =   rho_zeta_e*  (exp(zeta_e(-1))) + (1-rho_zeta_e)*zeta_e_ss + (1-rho_zeta_e)*chi_zeta_y*(exp(Y)/exp(Y(-1))-1) + e_zeta_e; %
%exp(zeta_i) =   rho_zeta_i*  (exp(zeta_i(-1))) + (1-rho_zeta_i)*zeta_i_ss + (1-rho_zeta_i)*chi_zeta_y*(exp(Y)/exp(Y(-1))-1) + e_zeta_i; %

exp(zeta_e) =   rho_zeta_e*  (exp(zeta_e(-1))) + (1-rho_zeta_e)*zeta_e_ss  + e_zeta_e; 
exp(zeta_i) =   rho_zeta_i*  (exp(zeta_i(-1))) + (1-rho_zeta_i)*zeta_i_ss  + e_zeta_i; 


////***********  10) AUXILIARY VARIABLES *****************************************************4

rr_e         = exp(lam_e) - beta_e*exp(lam_e(+1))*(1+exp(r_be))/exp(pie(+1));
RatioCap     = exp(K_m)/exp(B);
exp(bm)      = (exp(b_h(-1))/(exp(b_h(-1))+exp(b_e(-1))) * exp(r_bh(-1)) + exp(b_e(-1))/(exp(b_h(-1))+exp(b_e(-1))) * exp(r_be(-1))) - exp(r_d(-1));
exp(spr_b)   =  0.5*exp(r_bh) + 0.5*exp(r_be) - exp(r_d);                                 //64

end;





steady_state_model ; 
ee_z  =  log(1);
A_e = log(1);
ee_qk  =  log(1);
eps_K_m  =  log(1);
pie_wp  = log (1) ; 
pie_wi = log(1);
pie = log(piss) ;
u = log(1) ;
q_k = log(1);

mk_d   = log(mk_d_ss);
mk_be  = log(mk_be_ss);
mk_bh  = log (mk_bh_ss);
m_i    = log (m_i_ss);
m_e    = log( m_e_ss);
eps_y  = log(eps_y_ss);
eps_l  = log(eps_l_ss);
zeta_i   = log(zeta_i_ss);
zeta_e   = log(zeta_e_ss);
vi   =   log(vi_ss);
r_k  = log(r_k_ss) ;

r_ib = log(r_ib_ss) ;
r_bh = log(r_bh_ss) ;
r_be = log(r_be_ss) ;

x   = log( exp(eps_y)/(exp(eps_y)-1)) ;
r_d  = log ( exp(mk_d) * exp(r_ib)) ;
R_b  = r_ib ;
l_i = log (( (exp(eps_l) - 1)/(exp(eps_l) * Z_i))^(1/(1 + phi))) ;
l_id = l_i ;
l_p = log( ( (exp(eps_l) - 1)/(exp(eps_l) * Z_p))^(1 /(1 + phi)) );
l_pd = l_p ;

y_e = log ( ( ( exp(A_e) * (alpha/(exp(x)*exp(r_k)))^alpha)^(1/(1-alpha)) ) * (exp(l_p)^ni * exp(l_i)^(1-ni)) ) ;
w_i = log ((1-ni)*(1-alpha)*exp(y_e)/(exp(x)*exp(l_i)) );
w_p = log (ni*(1-alpha)*exp(y_e)/(exp(x)*exp(l_p)) );

k_e = log ( alpha * exp(y_e)/(exp(x) * exp(r_k)) ) ;
K = k_e ;
I = log(deltak * exp(K));

b_i =log ( exp(m_i) * exp(w_i)* exp(l_i) /(1 + exp(r_bh)) ) ; 
b_e = log (exp(m_e) * (1 - deltak) * exp(k_e)/(1 + exp(r_be))) ;
B = log (exp(b_i) + exp(b_e)) ;
d_p = log ( (1 - exp(vi)) *exp(B) ) ;
D = d_p ;
K_m = log( exp(vi) * exp(B) );
J_m = log ( deltakm * exp(K_m)/omega_m );
J_R = log ((1 - 1/exp(x))* exp(y_e)) ;
t_p = log (exp(J_R) + (1 - omega_m) * exp(J_m)) ;

c_p = log (exp(w_p) * exp(l_p) + exp(r_d) * exp(d_p) + exp(t_p)) ;
c_i = log (exp(w_i) * exp(l_i) + ( 1 - (1 - exp(zeta_i)) * (1 + exp(r_bh)) ) * exp(b_i)) ;
c_e = log ( alpha * exp(y_e) /exp(x) - exp(I) + ( 1 - (1 - exp(zeta_e)) * (1 + exp(r_be)) ) * exp(b_e));

lam_i = log( 1/exp(c_i)) ; 
lam_p = log( 1/exp(c_p)) ; 
lam_e = log( 1/exp(c_e)) ; 

//s_i = log (exp(lam_i) / (1 + exp(r_bh)) * ( 1 - beta_i * (1 - exp(zeta_i) )  )) ;
s_e = log (exp(lam_e) / (1 + exp(r_be)) * ( 1 - beta_e * (1 - exp(zeta_e) ) )) ;

b_ee = b_e ;
d_b = d_p ;
b_h = b_i; 

C     = log (exp(c_e)  + exp(c_p) + exp(c_i));
Y     = log (exp(C) + exp(I)) ;
BE    = b_e;
BH    = b_h ;
j_m   = J_m;
PIW  = log(1) ;
Y1     = y_e ; // log (exp(C) + exp(I));
rr_e   = exp(lam_e) - beta_e*exp(lam_e)*(1+exp(r_be));
RatioCap  = exp(K_m)/exp(B);
bm      = log (exp(b_h)/( exp(b_h)+exp(b_e)) * exp(r_bh) + exp(b_e)/( exp(b_h)+exp(b_e)) * exp(r_be) - exp(r_d));
spr_b   = log( 0.5*exp(r_bh) + 0.5*exp(r_be) - exp(r_d)); 

/*
interestPol   = 400*exp(r_ib); 
interestH = 400*exp(r_bh);
interestF = 400*exp(r_be);                                         
inflation  = pie*100;    
loansH = BH*100;
loansF = BE*100;
output     = Y1*100; 
consumption = C*100;
investment = I*100;   
deposits = D*100;
interestDep = 400*(exp(r_d));
IMFcapital = 100*K_m;                 
ROA        = exp(J_m)/exp(B) ;
*/
//**************************************************************************

end;


//check ;
//model_diagnostics ;

steady;

shocks;
// Variances des chocs (écarts-types au carré)
var e_z         = 0.0144^2;    // Choc technologique
var e_A_e       = 0.0062^2;    // Choc de productivité entrepreneurs
var e_me        = 0.0034^2;    // Choc mark-up entrepreneurs
var e_mi        = 0.0023^2;    // Choc mark-up intermédiaires
var e_mk_d      = 0.0488^2;    // Choc mark-up dépôts
var e_mk_bh     = 0.0051^2;    // Choc mark-up prêts ménages
var e_mk_be     = 0.1454^2;    // Choc mark-up prêts entreprises
var e_qk        = 0.0128^2;    // Choc prix du capital
var e_r_ib      = 0.0018^2;    // Choc taux interbancaire
var e_y         = 1.0099^2;    // Choc production
var e_l         = 0.3721^2;    // Choc travail
var e_eps_K_m   = 0.637^2; %0.05  // Choc fonds propres des microfinances
var e_zeta_e      = 0.050^2;     // Choc risque de crédit chez les entrepreneurs
var e_zeta_i      = 0.050^2;     // Choc risque de crédit chez les ménages

end;



