$Title Storage Investment Model

$OnText
Developed by

   Diego Alejandro Tejada Arango
   Instituto de Investigación Tecnológica
   Escuela Técnica Superior de Ingeniería - ICAI
   UNIVERSIDAD PONTIFICIA COMILLAS
   Alberto Aguilera 23
   28015 Madrid, Spain
   dtejada@comillas.edu

$OffText

*-------------------------------------------------------------------------------
*                              MODEL OPTIONS
*-------------------------------------------------------------------------------

* To use a algebraic formulation (formulation without defining the sets)
$OnEmpty OnMulti OffListing

* definition of symbol for comments at the end of the line
$EOLCOM //

* optimizer definition
OPTION   lp = cplex ;
OPTION  mip = cplex ;
OPTION rmip = cplex ;

* general options
OPTION optcr    =      0 ;   // tolerance to solve MIP until IntGap < OptcR
OPTION reslim   = 129600 ;   // maximum run time [sec]
OPTION threads  =     -1 ;   // number of cores
OPTION solprint =     on ;   // print the final solution in the .lst file
OPTION limrow   =   1000 ;   // maximum number of equations in the .lst file
OPTION limcol   =   1000 ;   // maximum number of variables in the .lst file

* profile options
OPTION profile=1, profileTol = 0.01 ;

* Statement of sets and indices
SETS
p             Time Periods
g             Generators
n             Electrical nodes
tech          Set of types of generator technologies
c             Circuit ID
s             System states
ch            State changes
rp            Representative periods (e.g. days weeks)
* Statement of dynamic sets: they are subsets of the previous sets.
t       (g     ) Thermal Generators
h       (g     ) Hydraulic (or storage) Generators
hr      (g     ) Hydraulic (or storage) Generators with        reservoir
hf      (g     ) Hydraulic (or storage) Generators with "fast" reservoir (i.e. battery)
hs      (g     ) Hydraulic (or storage) Generators with "slow" reservoir (i.e. hydro  )
pa      (p     ) Active periods
ps      (p     ) Active periods for "slow" storage in RP-TM model
ns      (n     ) Nodes without slack
gn      (g,n   ) Generators connected to node n
gtech   (g,tech) generators of technology type tech
hindex  (p,s   ) Hour index indicating which hour belongs to which state
hindexRP(p,p   ) Hour index indicating which hour belongs to which hour of representative period
Omega2  (n,n   ) Set with nodes with at least one line connection
Omega   (n,n,c ) Set with connection between nodes
rpp     (rp ,p ) relation among representative periods and periods

ALIAS(n,m),(s,ss), (p,pp,ppp), (rp,rrpp) ;

* Statement of parameters
PARAMETERS
* General Parameters
pISF     (n,m,c,n) Injection Shift Factors (for node i to line i->ii circuitID c)
pXnm     (n,m,c)   Reactance between nodes [GW]
pTCmax   (n,m,c)   Maximum transmission capacity for a line[GW]
pYBUS    (n,m)     Susceptance matrix
pYBUSInv (n,m)     Susceptance matrix inverse Z = B^-1
pInflows (p,g)     Inflows to the reservoir of generator g in period p [GWh]
pDemand  (p,n)     Hourly demand per node [GW]
pWind    (p,n)     Available wind  production data per period and node [GW]
pSolar   (p,n)     Available solar production data per period and node [GW]
pUo      (g)       Initial status at the beginning of the first hour {1 0}
pAlfa    (g)       Variable fuel consumption of generator g [MTh per GWh]
pBeta    (g)       Fixed fuel conxumption fo generator g [MTh per h]
pGamma   (g)       Fuel consumption of generator g at startup  [MTh]
pTheta   (g)       Fuel consumption of generator g at shutdown [MTh]
pQmax    (g)       Maximum gross power of generator g [GW]
pQmin    (g)       Minimum gross power of generator g [GW]
pBmax    (g)       Maximum gross pumping power of generatior g [GW]
pWmax    (g)       Maximum reserve level of catchment basin of generator g [GWh]
pWmin    (g)       Minimum reserve level of the catchment basin of generator g [GWh]
pEta     (g)       Efficiency of the pumping-turbine cycle of generator g [p.u.]
pCostfuel(g)       Cost of fuel used by generator g [k€ per MTh]
pCostOM  (g)       Variable cost of operation and maintenaince of the generator g [kEuros per GWh]
pW0      (g)       Initial reserve level of catchment basin for generator g [GWh]
pWfin    (g)       Allocation of final reserve level of the catchment basin for generator g [GWh]
pSRR     (g)       Ten minute ramp of g (for spinning reserve) [GW per 10 min]
pK       (g)       Gross to net power conversion factor for generator g [p.u.]
pXres              Operating reserve as fraction of total demand [p.u.]
pCostENS           Cost of energy non-supplied  [k€ per GWh]
pTransNet          Parameter to control if the model considers or not the transmission network {0 1}
pSbase             Base power for p.u. calculations
pSlack             Slack node
pModel             Indicator of model (1-5) we solve {1-hourly 2-SystemStates 3-SS_ReducFmatrix 4-RepreDays 5-RP+TM}

* Parameters for storage investment
pStorageInvest     Parameter to control if the model considers or not storage investment {0 1}
pInveCost      (g) Storage investment cost       per storage technology [k€ per GW]
pMaxStorI      (g) Maximum storage investment    per storage technology [       GW]
pEPRmax         (g) Maximum Energy to Power Ratio per storage technology [h        ]
pEPRmin        (g) Minimum Energy to Power Ratio per storage technology [h        ]

* Parameters for system states formulation
pFrequency  (s,ss,ch) Frequency with which the state change from s to ss appears at state change ch
pFreqMatRed (s,ss,ch) Frequency  matrix reduced with which the state change from s to ss appears at state change ch
pTransMatrix(s,ss)    Transition matrix representing number of jumps from s to ss
pDemand_s   (s,n)     Demand per system state per node [GW]
pWind_s     (s,n)     Wind  production data per system state [GW]
pSolar_s    (s,n)     Solar production data per system state [GW]
pInflows_s  (s,g)     Inflows to the reservoir of generator g in system state s [GWh]
pDuration_s (s)       Duration of system state s [h]
pHourChanges(ch)      Hour after which a state changes
pHourChRFM  (ch)      Hour after which a state changes for reduced frequency matrix

* Parameters for representatives periods formulation
pWeight_rp      (rp   ) Representatives periods weight [h]
pNumPer_rp      (rp   ) Number of periods at each representative period [h]
pFirstP_rp      (rp,p ) First period of each representative period
pTransMatrix_rp (rp,rp) Transition matrix representing number of jumps from rp to rrpp
pNumPerMinus1_rp(rp   ) Number of periods minus 1 at each representative period [h]
pMaxTM_rp               Maximum value of representative period transition matrix
pMinTransition          Minimum transition to be considered in RP_TM model
pStorMovWindow          Storage moving window for representative periods [h]
pLastStorMovW           parameter to indicate the last hour on ps(p) set
pScalar                 To activate periods to calculate "slow" storage level
;

* Statement of free variables
VARIABLES
vObjectiveFunction             Value of objective function
vCircuitPowerFlow   (p,n,n,c)  Electrical power flow between node n and m per period [GW]
* Variables for system states formulation
vCircuitPowerFlow_s (s,n,n,c)  Electrical power flow between node n and m per system state s [GW]
vDeltaStoReserve    (s,ss,g)   Difference of storage level from the beginning of an hour in s to the beginning of the hour in ss [GWh]

* Statement of positive variables
POSITIVE VARIABLES
vProduction     (p,g)    Power production by generator g per period [GW]
vRampProduction (p,g)    Power production by generator g above minimum stable load [GW]
vConsumption    (p,g)    Gross power consumed by the generator g when acting as a pump [GW]
vSpinningReserve(p,g)    Spinning reserve provided by generator g [GW]
vStorageReserve (p,g)    Energy storage in the reservoir of the generator g at the end of p [GWh]
vSpillage       (p,g)    Spillages of the reservoir of generator g during period p [GW]
vWind           (p,n)    Wind  production per period and node [GW]
vSolar          (p,n)    Solar production per period and node [GW]
vPowerNotSupply (p,n)    Power not supply per period and node [GW]
vStorageInvest  (  g)    Storage investment [GW]

* Variables for system states levels formulation
vProduction_s     (s,g)  Power production by generator g in system state s [GW]
vRampProduction_s (s,g)  Power dispatched by generator g above minimum stable load [GW]
vConsumption_s    (s,g)  Gross power consumed by the generator g when acting as a pump [GW]
vSpinningReserve_s(s,g)  Spinning reserve provided by generator g per load level[GW]
vStorageReserve_s (s,g)  Energy storage in the reservoir of the generator g at the end of system state s [GWh]
vSpillage_s       (s,g)  Spillages of the reservoir of generator g during system state s [GW]
vWind_s           (s,n)  Wind  production per system state and node [GW]
vSolar_s          (s,n)  Solar production per system state and node [GW]
vPowerNotSupply_s (s,n)  Power not supply per system state and node [GW]

* Statement of binary variables
BINARY VARIABLES
vConnected   (p,   g)    Binary variable indicating whether unit g is connected (1) or disconnected (0) in period
vStartUp     (p,   g)    Start-up decision for unit g in period p [binary]

* Variables for system states formulation
vConnected_s (s,   g)    Binary variable indicating whether unit g is connected (1) or disconnected (0) in system state
vStartUp_s   (s,ss,g)    Start-up decision for unit g of state ss that started in s and ended in ss

* Statement of Equations
EQUATIONS
E_FOBJ           Objective Function [1000 M€]
E_DMND (p  )     Meet the demand of the system (only if one node model is executed)
E_DMNDN(p,n)     Meet the demand per node
E_RROD1(p,g)     Spinning reserve for a thermal unit
E_RROD2(p,g)     Spinning reserve (10 minute reserve)
E_RROD3(p  )     Spinning reserve (overhead)
E_ACOP (p,g)     Startup decision
E_INI  (  g)     Initial condition for thermal group
E_QMAXT(p,g)     Maximum power of thermal group
E_QMINT(p,g)     Minimum power of thermal group
E_RSRVH(p,g)     Evolution of reserves g
E_PF   (p,n,n,c) Power flow between nodes
E_QMAXS(p,g)     Maximum production  of storage technology considering investment
E_BMAXS(p,g)     Maximum consumption of storage technology considering investment
E_WMAX (p,g)     Maximum storage level considering investment
E_WMIN (p,g)     Minimum storage level considering investment

* Equations for system states formulation
E_FOBJ_S              Objective Function with system states [1000 M€]
E_DMND_S    (s    )   Meet the demand per system state
E_DMNDN_S   (s,  n)   Meet the demand per node per system state
E_RROD1_S   (s,  g)   Spinning reserve for a thermal unit
E_RROD2_S   (s,  g)   Spinning reserve (10 minute reserve)
E_RROD3_S   (s    )   Spinning reserve (overhead)
E_ACOP_S    (s,s,g)   Logic of start ups and shut downs of thermal group between system states
E_INI_S     (    g)   Initial condition for thermal group
E_QMAXT_S   (s,  g)   Maximum power of thermal group
E_QMINT_S   (s,  g)   Minimum power of thermal group
E_DEF_DELTAW(s,s,g)   Definition of variable deltaw
E_UB_DELTAW (    g)   Upper bound for overall time horizon
E_LB_DELTAW (    g)   Lower bound for overall time horizon
E_LB_APPROX (ch, g)   Lower bound for overall storage state change
E_UB_APPROX (ch, g)   Upper bound for overall storage state change
E_QMAXS_S   (s,  g)   Maximum production  of storage technology considering investment
E_BMAXS_S   (s,  g)   Maximum consumption of storage technology considering investment

* System states model reduced F matrix
E_LB_APPROX2(ch, g)   Lower bound for overall storage state change for reduced F matrix
E_UB_APPROX2(ch, g)   Upper bound for overall storage state change for reduced F matrix
E_LB_APPROX3(ch, g)   Lower bound for overall storage state change for reduced F matrix
E_UB_APPROX3(ch, g)   Upper bound for overall storage state change for reduced F matrix
E_PF_S      (s,n,n,c) Power flow between nodes per system state

* Equations for representative periods formulation
E_FOBJ_RP             Objective Function for representative periods [1000 M€]
E_LB_STOR_RP(rp,p,g)  Lower bound for storage units for representative periods

* Equation for representative periods formulation with transition matrix
E_ACOP_RP      (rp,p,rp,g) Logic of start ups and shut downs of thermal group between representatives
E_LB_STOR_RPTM (rp,p,rp,g) Lower bound for storage units for representative periods using Transition Matrix
E_RSRVH_RP1    (   p   ,g) Evolution of reserves for storage units with a moving window
E_RSRVH_RP2    (   p   ,g) Evolution of reserves for storage units with a moving window for last period
E_WMAX_RP1      (   p   ,g) Maximum storage level considering investment
E_WMIN_RP1      (   p   ,g) Minimum storage level considering investment
E_WMAX_RP2      (   p   ,g) Maximum storage level considering investment
E_WMIN_RP2      (   p   ,g) Minimum storage level considering investment

;

*-------------------------------------------------------------------------------
*                                EQUATIONS
*-------------------------------------------------------------------------------

* ------------------------------ hourly model ----------------------------------
E_FOBJ ..
  vObjectiveFunction =E= SUM[pa(p),
*                            penalize the power not supply
                             pCostENS * SUM[n, vPowerNotSupply(p,n)]
*                            penalize spillages
                             + SUM[hr, vSpillage(p,hr) * 0.0001 ]
*                            Thermal production cost
                             + SUM[t,
                                  pCostfuel(t) * [ pBeta(t) * vConnected(p,t)
                                         + (pGamma(t)+pTheta(t)) * vStartUp(p,t)
                                         + pAlfa(t) * vProduction(p,t) / pK(t)
                                         ]
                                  + pCostOM(t) * vProduction(p,t) / pK(t)
                                  ]
                            ] * 1e-6                                         //1e-6 to convert k€ -> 1000 M€
                        + SUM[h, vStorageInvest(h) * pInveCost(h) ] * 1e-6 ; //1e-6 to convert k€ -> 1000 M€

E_DMND(pa(p))    $ [NOT pTransNet] ..
*              Thermal production in the node
                 SUM[t, vProduction(p,t)]
*              Hydro production and charging/pumping consumption in the node
               + SUM[h, vProduction(p,h)-vConsumption(p,h)]
*              Wind  production
               + SUM[n, vWind (p,n)]
*              Solar production
               + SUM[n, vSolar(p,n)]
*              Power not supply in the node
               + SUM[n, vPowerNotSupply(p,n)]
*              Demand in the node
               =E= SUM[n, pDemand(p,n)] ;

E_DMNDN(pa(p),n) $ pTransNet ..
*              Thermal production in the node
                 SUM[gn(t,n), vProduction(p,t)]
*              Hydro production and charging/pumping consumption in the node
               + SUM[gn(h,n), vProduction(p,h)-vConsumption(p,h)]
*              Wind  production
               + vWind (p,n)
*              Solar production
               + vSolar(p,n)
*              Power flow entering through the lines connected to the node
               + SUM[Omega(m,n,c), vCircuitPowerFlow(p,m,n,c)]
*              Power flow outgoing through the lines connected to the node
               - SUM[Omega(n,m,c), vCircuitPowerFlow(p,n,m,c)]
*              Power not supply in the node
               + vPowerNotSupply(p,n)
*              Demand in the node
               =E= pDemand(p,n) ;

E_RROD1(pa(p),t)..
               vSpinningReserve(p,t) + vProduction(p,t) =L= vConnected(p,t) * pQmax(t) * pK(t) ;

E_RROD2(pa(p),t)..
               vSpinningReserve(p,t) =L= pSRR(t) ;

E_RROD3(pa(p))  ..
               SUM[t, vSpinningReserve(p,t)] =G= pXres * SUM[n, pDemand(p,n) - vPowerNotSupply(p,n)] ;

E_ACOP(pa(p),t) $[ORD(p)>1] ..
               vStartUp(p,t) =G= vConnected(p,t) - vConnected(p-1,t) ;

E_INI(t) ..
               SUM{p $[ORD(p)=1],vConnected(p,t)} =E= SUM{p $[ORD(p)=CARD(p)], vConnected(p,t)} ;

E_QMAXT(pa(p),t)..
               vProduction(p,t) =E= vConnected(p,t) * pQmin(t) * pK(t) + vRampProduction(p,t) ;

E_QMINT(pa(p),t)..
               vRampProduction(p,t)=L= (pQmax(t)-pQmin(t)) * pK(t) * vConnected(p,t) ;

E_RSRVH(pa(p),hr(h))..
               vStorageReserve(p,h) =E=
*              reserve value at p-1
               vStorageReserve(p-1,h)$[ORD(p) > 1] +  pW0(h)$[ORD(p)=1]
*              Hydro production
               - vProduction(p,h)
*              Spillage
               - vSpillage(p,h)
*              charging/pumping consumption
               + pEta(h) * vConsumption(p,h)
*              inflows
               + pInflows(p,h) ;

E_PF(pa(p),n,m,c) $[Omega(n,m,c)] ..
             vCircuitPowerFlow(p,n,m,c) =E= SUM[ns, pISF(n,m,c,ns)            *
                                            [SUM[gn(t,ns), vProduction (p,t)] +
                                             SUM[gn(h,ns), vProduction (p,h)] -
                                             SUM[gn(h,ns), vConsumption(p,h)] +
                                             vWind  (p,ns)                    +
                                             vSolar (p,ns)                    -
                                             pDemand(p,ns)]] ;

E_QMAXS  (pa(p),h)..
           vProduction    (p,h) =L= pk   (h) * [pQmax(h)   + vStorageInvest(h)];

E_BMAXS(pa(p),h)..
           vConsumption   (p,h) =L= pBmax(h) +  pEta (h)   * vStorageInvest(h) ;

E_WMAX (pa(p),h)..
           vStorageReserve(p,h) =L= pWmax(h) +  pEPRmax(h) * vStorageInvest(h) ;

E_WMIN (pa(p),h)..
           vStorageReserve(p,h) =G= pWmin(h) +  pEPRmin(h) * vStorageInvest(h) ;

* ------------------------- system states model --------------------------------
E_FOBJ_S ..
  vObjectiveFunction =E= SUM[s,
*                            penalize the power not supply
                             pCostENS * SUM[n, vPowerNotSupply_s(s,n)] * pDuration_s(s)
*                            penalize spillages
                             + SUM[hr, vSpillage_s(s,hr) * 0.0001 * pDuration_s(s)]
*                            Thermal production cost
                             + SUM[t,
                                  pCostfuel(t) * [ pBeta(t) * vConnected_s(s,t) * pDuration_s(s)
                                         + SUM[ss$(ORD(ss) <> ORD(s)),(pGamma(t)+pTheta(t)) * vStartUp_s(ss,s,t) * pTransMatrix(ss,s)]
                                         + pAlfa(t) * vProduction_s(s,t) / pK(t) * pDuration_s(s)
* preguntar a Sonia, porque pk(t) si es una eficiencia, porque divide si el dato de entrada es menor a cero.
                                                 ]
                                  + pCostOM(t) * vProduction_s(s,t) / pK(t) * pDuration_s(s)
                                  ]
                            ] * 1e-6                                         //1e-6 to convert k€ -> 1000 M€
                        + SUM[h, vStorageInvest(h) * pInveCost(h) ] * 1e-6 ; //1e-6 to convert k€ -> 1000 M€
E_DMND_S(s) $ [NOT pTransNet] ..
*              Thermal production in the node
                 SUM[t, vProduction_s(s,t)]
*              Hydro production and charging/pumping consumption in the node
               + SUM[h, vProduction_s(s,h)-vConsumption_s(s,h)]
*              Wind  production
               + SUM[n, vWind_s (s,n)]
*              Solar production
               + SUM[n, vSolar_s(s,n)]
*              Power not supply in the node
               + SUM[n, vPowerNotSupply_s(s,n)]
*              Demand in the node
               =E= SUM[n, pDemand_s(s,n)] ;

E_DMNDN_S(s,n) $ pTransNet ..
*              Thermal production in the node
                 SUM[gn(t,n), vProduction_s(s,t)]
*              Hydro production and charging/pumping consumption in the node
               + SUM[gn(h,n), vProduction_s(s,h)-vConsumption_s(s,h)]
*              Wind  production
               + vWind_s (s,n)
*              Solar production
               + vSolar_s(s,n)
*              Power flow entering through the lines connected to the node
               + SUM[Omega(m,n,c), vCircuitPowerFlow_s(s,m,n,c)]
*              Power flow outgoing through the lines connected to the node
               - SUM[Omega(n,m,c), vCircuitPowerFlow_s(s,n,m,c)]
*              Power not supply in the node
               + vPowerNotSupply_s(s,n)
*              Demand in the node
               =E= pDemand_s(s,n) ;

E_RROD1_S(s,t)..
               vSpinningReserve_s(s,t) + vProduction_s(s,t) =L= vConnected_s(s,t) * pQmax(t) * pK(t) ;

E_RROD2_S(s,t)..
               vSpinningReserve_s(s,t) =L= pSRR(t) ;

E_RROD3_S(s)  ..
               SUM[t, vSpinningReserve_s(s,t)] =G= pXres * SUM[n, pDemand_s(s,n) - vPowerNotSupply_s(s,n)] ;

E_ACOP_S (s,ss,t) $[ORD(ss) <> ORD(s) AND pTransMatrix(ss,s)>0]..
               vStartUp_s(ss,s,t) =G= vConnected_s(s,t) - vConnected_s(ss,t) ;

E_INI_S(t) ..
               SUM{s $[ORD(s)=1],vConnected_s(s,t)} =E= SUM{s $[ORD(s)=CARD(s)], vConnected_s(s,t)} ;

E_QMAXT_S(s,t)..
               vProduction_s(s,t) =E= vConnected_s(s,t) * pQmin(t) * pK(t) + vRampProduction_s(s,t) ;

E_QMINT_S(s,t)..
               vRampProduction_s(s,t) =L= (pQmax(t)-pQmin(t)) * pK(t) * vConnected_s(s,t) ;

E_DEF_DELTAW(s,ss,hr) $[pTransMatrix(s,ss)>0] ..

  vDeltaStoReserve(s,ss,hr) =E=  0.5 * [- vProduction_s(s, hr) - vSpillage_s(s, hr) + pEta(hr) * vConsumption_s(s, hr) + pInflows_s(s, hr)]
                               + 0.5 * [- vProduction_s(ss,hr) - vSpillage_s(ss,hr) + pEta(hr) * vConsumption_s(ss,hr) + pInflows_s(ss,hr)] ;

E_UB_DELTAW (hr)..
               SUM[(s,ss), vDeltaStoReserve(s,ss,hr) * pTransMatrix(s,ss)]
* finite differences at center of each state
               + 0.5 * SUM{s$[ORD(s)=CARD(s) OR ORD(s)=1],- vProduction_s(s, hr) - vSpillage_s(s, hr) + pEta(hr) * vConsumption_s(s, hr) + pInflows_s(s, hr)}
               =L= pWmax(hr) - pW0(hr) + pEPRmax(hr) * vStorageInvest(hr) ;

E_LB_DELTAW (hr)..
               SUM[(s,ss), vDeltaStoReserve(s,ss,hr) * pTransMatrix(s,ss)]
* finite differences at center of each state: for first and last half hour we add:
               + 0.5 * SUM{s$[ORD(s)=CARD(s) OR ORD(s)=1],- vProduction_s(s, hr) - vSpillage_s(s, hr) + pEta(hr) * vConsumption_s(s, hr) + pInflows_s(s, hr)}
               =G= pWfin(hr) - pW0(hr) + pEPRmin(hr) * vStorageInvest(hr) ;

E_LB_APPROX (ch,hr) ..
               SUM{(s,ss)$[pFrequency(s,ss,ch)>0], vDeltaStoReserve(s,ss,hr) * pFrequency(s,ss,ch)}
** finite differences at center of each state: for first and last half hour we add:
               + 0.5 * SUM{s$[ORD(s)=1],- vProduction_s(s, hr) - vSpillage_s(s, hr) + pEta(hr) * vConsumption_s(s, hr) + pInflows_s(s, hr)}
               =G= pWmin(hr) - pW0(hr) + pEPRmin(hr) * vStorageInvest(hr) ;

E_UB_APPROX (ch,hr) ..
               SUM{(s,ss)$[pFrequency(s,ss,ch)>0], vDeltaStoReserve(s,ss,hr) * pFrequency(s,ss,ch)}
** finite differences at center of each state: for first and last half hour we add:
               + 0.5 * SUM{s$[ORD(s)=1],- vProduction_s(s, hr) - vSpillage_s(s, hr) + pEta(hr) * vConsumption_s(s, hr) + pInflows_s(s, hr)}
               =L= pWmax(hr) - pW0(hr) + pEPRmax(hr) * vStorageInvest(hr) ;

E_PF_S(s,n,m,c) $[Omega(n,m,c)] ..
           vCircuitPowerFlow_s(s,n,m,c) =E= SUM[ns, pISF(n,m,c,ns)              *
                                            [SUM[gn(t,ns), vProduction_s (s,t)] +
                                             SUM[gn(h,ns), vProduction_s (s,h)] -
                                             SUM[gn(h,ns), vConsumption_s(s,h)] +
                                             vWind_s  (s,ns)                    +
                                             vSolar_s (s,ns)                    -
                                             pDemand_s(s,ns)]] ;

E_QMAXS_S(s,h)..
           vProduction_s  (s,h) =L= pk   (h) * [pQmax(h)   + vStorageInvest(h)];

E_BMAXS_S(s,h)..
           vConsumption_s (s,h) =L= pBmax(h) +  pEta (h)   * vStorageInvest(h) ;

* ---------------- system states model reduced F matrix-------------------------

E_LB_APPROX2 (ch,hs) $ pHourChRFM(ch)..
               SUM{(s,ss)$[pFrequency(s,ss,ch)>0], vDeltaStoReserve(s,ss,hs) * pFrequency(s,ss,ch)}
** finite differences at center of each state: for first and last half hour we add:
               + 0.5 * SUM{s$[ORD(s)=1],- vProduction_s(s, hs) - vSpillage_s(s, hs) + pEta(hs) * vConsumption_s(s, hs) + pInflows_s(s, hs)}
               =G= pWmin(hs) - pW0(hs) + pEPRmin(hs) * vStorageInvest(hs) ;

E_UB_APPROX2 (ch,hs) $ pHourChRFM(ch)..
               SUM{(s,ss)$[pFrequency(s,ss,ch)>0], vDeltaStoReserve(s,ss,hs) * pFrequency(s,ss,ch)}
** finite differences at center of each state: for first and last half hour we add:
               + 0.5 * SUM{s$[ORD(s)=1],- vProduction_s(s, hs) - vSpillage_s(s, hs) + pEta(hs) * vConsumption_s(s, hs) + pInflows_s(s, hs)}
               =L= pWmax(hs) - pW0(hs) + pEPRmax(hs) * vStorageInvest(hs) ;

E_LB_APPROX3(ch,hf) ..
               SUM{(s,ss)$[pFreqMatRed(s,ss,ch)>0],vDeltaStoReserve(s,ss,hf) * pFreqMatRed(s,ss,ch)}
               =G= pWmin(hf) - pW0(hf) + pEPRmin(hf) * vStorageInvest(hf) ;

E_UB_APPROX3(ch,hf) ..
               SUM{(s,ss)$[pFreqMatRed(s,ss,ch)>0],vDeltaStoReserve(s,ss,hf) * pFreqMatRed(s,ss,ch)}
               =L= pWmax(hf) - pW0(hf) + pEPRmax(hf) * vStorageInvest(hf) ;

* --------------------- representative periods model ---------------------------
E_FOBJ_RP ..
  vObjectiveFunction =E= SUM[rpp(rp,pa(p)), pWeight_rp(rp) *
                                [
*                               penalize the power not supply
                                pCostENS * SUM[n, vPowerNotSupply(p,n)]
*                               penalize spillages
                                + SUM[h, vSpillage(p,h) * 0.0001 ]
*                               Thermal production cost
                                + SUM[t,
                                     + pCostfuel(t) *
                                         [ pBeta(t) * vConnected(p,t)
                                         + (pGamma(t)+pTheta(t)) * vStartUp(p,t)
                                         + pAlfa(t) * vProduction(p,t) / pK(t)
                                         ]
                                     + pCostOM(t) * vProduction(p,t) / pK(t)
                                     ]
                                ]
                            ] * 1e-6                                         //1e-6 to convert k€ -> 1000 M€
                        + SUM[h, vStorageInvest(h) * pInveCost(h) ] * 1e-6 ; //1e-6 to convert k€ -> 1000 M€

E_LB_STOR_RP(rp,p,hr(h)) $[pFirstP_rp(rp,p) AND pa(p)] ..
         SUM[pp $[ORD(pp)= pFirstP_rp(rp,p)+pNumPer_rp(rp)-1],
             vStorageReserve(pp,h)] =G= vStorageReserve(p,h) ;

* ------------ representative periods model with Transition Matrix -------------

E_ACOP_RP      (rp,p,rrpp,t    )
  $[pa(p) AND pFirstP_rp     (rp,p   )
          AND pTransMatrix_rp(rp,rrpp) > pMinTransition]..

  vConnected    (p,t) =E=
                    SUM[pp $[pFirstP_rp(rrpp,pp)],
                        vConnected     (pp+pNumPerMinus1_rp(rrpp),t)] ;

E_RSRVH_RP1(ps(p),h) $[ORD(p) < CARD(p)]..
               vStorageReserve(p,h) =E=
               // Reserve value at p-pStorMovWindow (Moving Window)
               vStorageReserve(p-pStorMovWindow,h)$[ORD(p) >= pStorMovWindow]
               // Reserve value for first period
            +  pW0(h)$[ORD(p)=1]
               // Sum of production/consumption during the moving window
            +  SUM[pp$[ORD(pp) >= ORD(p)+1-pStorMovWindow
                   AND ORD(pp) <= ORD(p)],
               // Sum taking into acount the RP relation among periods
               +  SUM[hindexRP(pp,ppp),
               //     Hydro production
                    - vProduction(ppp,h)
               //     Spillage
                    - vSpillage  (ppp,h)
               //     charging/pumping consumption
                    + pEta(h) * vConsumption(ppp,h)
               //     inflows
                    + pInflows(ppp,h)]
                  ] ;

E_RSRVH_RP2(ps(p),h) $[ORD(p) = CARD(p)]..
               vStorageReserve(p,h) =E=
               // Reserve value at last period os set ps(p)
               vStorageReserve(p-pLastStorMovW,h)
               // Sum of production/consumption during the moving window
            +  SUM[pp$[ORD(pp) >= ORD(p)+1-pLastStorMovW
                   AND ORD(pp) <= ORD(p)],
               // Sum taking into acount the RP relation among periods
               +  SUM[hindexRP(pp,ppp),
               //     Hydro production
                    - vProduction(ppp,h)
               //     Spillage
                    - vSpillage  (ppp,h)
               //     charging/pumping consumption
                    + pEta(h) * vConsumption(ppp,h)
               //     inflows
                    + pInflows(ppp,h)]
                  ] ;

E_WMAX_RP1 (ps(p),hs(h))..
           vStorageReserve(p,h) =L= pWmax(h) +  pEPRmax(h) * vStorageInvest(h) ;

E_WMIN_RP1 (ps(p),hs(h))..
           vStorageReserve(p,h) =G= pWmin(h) +  pEPRmin(h) * vStorageInvest(h) ;

E_WMAX_RP2 (pa(p),hf(h))..
           vStorageReserve(p,h) =L= pWmax(h) +  pEPRmax(h) * vStorageInvest(h) ;

E_WMIN_RP2 (pa(p),hf(h))..
           vStorageReserve(p,h) =G= pWmin(h) +  pEPRmin(h) * vStorageInvest(h) ;



*-------------------------------------------------------------------------------
*                                MODELS
*-------------------------------------------------------------------------------

* Specification of equations to be used in the model
MODEL MSEM
/
E_FOBJ
E_DMND
E_DMNDN
E_RROD1
E_RROD2
E_RROD3
E_ACOP
E_INI
E_QMAXT
E_QMINT
E_RSRVH
E_PF
E_QMAXS
E_BMAXS
E_WMAX
E_WMIN
/;
MSEM.holdfixed = 1 ; MSEM.optfile = 1 ;

* Model with system states
MODEL MSEM_S
/
E_FOBJ_S
E_DMND_S
E_DMNDN_S
E_RROD1_S
E_RROD2_S
E_RROD3_S
E_ACOP_S
E_INI_S
E_QMAXT_S
E_QMINT_S
E_DEF_DELTAW
E_UB_DELTAW
E_LB_DELTAW
E_LB_APPROX
E_UB_APPROX
E_PF_S
E_QMAXS_S
E_BMAXS_S
/;
MSEM_S.holdfixed = 1 ; MSEM_S.optfile = 1 ;

* Model with system states
MODEL MSEM_SS_RFM
/
  MSEM_S
- E_LB_APPROX
- E_UB_APPROX
+ E_LB_APPROX2
+ E_UB_APPROX2
+ E_LB_APPROX3
+ E_UB_APPROX3
/;
MSEM_SS_RFM.holdfixed = 1 ; MSEM_SS_RFM.optfile = 1 ;

* Model with representative periods
MODEL MSEM_RP
/
  MSEM
- E_FOBJ
+ E_FOBJ_RP
+ E_LB_STOR_RP
/;
MSEM_RP.holdfixed = 1 ; MSEM_RP.optfile = 1 ;

* Model with representative periods and transition matrix
MODEL MSEM_RP_TM
/
  MSEM_RP
+ E_ACOP_RP
- E_LB_STOR_RP
+ E_RSRVH_RP1
+ E_RSRVH_RP2
- E_WMAX
- E_WMIN
+ E_WMAX_RP1
+ E_WMIN_RP1
+ E_WMAX_RP2
+ E_WMIN_RP2
/;

MSEM_RP_TM.holdfixed = 1 ; MSEM_RP_TM.optfile = 1 ;

*-------------------------------------------------------------------------------
*                   OPTIONS DEFINITION FOR SOLVERS
*-------------------------------------------------------------------------------
FILE     GOPT / gurobi.opt /
PUT      GOPT 'Method   2' / 'IntFeasTol 1E-6' / 'OptimalityTol 1E-6' /
PUT      GOPT 'FeasibilityTol 1E-9' / 'IIS 1' / 'Crossover 0'    / 'Names  1' /
PUTCLOSE GOPT

FILE     COPT / cplex.opt /
PUT      COPT 'LpMethod 4' / 'Eprhs 1E-6'      / 'Epopt 1E-6'         /
PUT      COPT 'Scaind 0'            / 'IIS yes'/ 'BarCrossalg -1' / 'Names YES' /
PUTCLOSE COPT

*-------------------------------------------------------------------------------
*                             DATA FROM EXCEL FILE
*-------------------------------------------------------------------------------
FILE TMP / tmp.txt /
$onecho > tmp.txt
   i="%gams.user1%.xlsm"
   r1=Demand
   o1=Demand
   r2=Demand_s
   o2=Demand_s
   r3=Duration
   o3=Duration
   r4=Frequency
   o4=Frequency
   r5=Hindex
   o5=Hindex
   r6=HourChanges
   o6=HourChanges
   r7=Indices
   o7=Indices
   r8=Inflows
   o8=Inflows
   r9=Inflows_s
   o9=Inflows_s
   r10=Network
   o10=Network
   r11=Param
   o11=Param
   r12=Set_gn
   o12=Set_gn
   r13=StorageUnits
   o13=StorageUnits
   r14=ThermalUnits
   o14=ThermalUnits
   r15=TransitionMatrix
   o15=TransitionMatrix
   r16=Wind
   o16=Wind
   r17=Wind_s
   o17=Wind_s
   r18=Solar
   o18=Solar
   r19=Solar_s
   o19=Solar_s
   r20=ActivePeriods
   o20=ActivePeriods
   r21=ReprePeriods
   o21=ReprePeriods
   r22=Weight
   o22=Weight
   r23=Hindex_rp
   o23=Hindex_rp
   r24=Set_gtech
   o24=Set_gtech
   r25=HourChangesRFM
   o25=HourChangesRFM
   r26=TransitionMatrix_rp
   o26=TransitionMatrix_rp
$offecho
$call xls2gms m @"tmp.txt"

* Sets used to define the structure and size of the problem
SETS
$INCLUDE Indices
$INCLUDE Set_gn
$INCLUDE Set_gtech
$INCLUDE Hindex
$INCLUDE ReprePeriods
$INCLUDE Hindex_rp
;

PARAMETERS
$INCLUDE Duration
$INCLUDE Frequency
$INCLUDE HourChanges
$INCLUDE HourChangesRFM
$INCLUDE Weight
;

* Parameters (penalizations and flags)
$INCLUDE Param

* Information from tables
TABLE    pDemand(p,n)
$INCLUDE Demand
TABLE    pDemand_s(s,n)
$INCLUDE Demand_s
TABLE    pWind(p,n)
$INCLUDE Wind
TABLE    pWind_s(s,n)
$INCLUDE Wind_s
TABLE    pSolar(p,n)
$INCLUDE Solar
TABLE    pSolar_s(s,n)
$INCLUDE Solar_s
TABLE    pInflows(p,g)
$INCLUDE Inflows
TABLE    pInflows_s(s,g)
$INCLUDE Inflows_s
TABLE    pThermalUnits(g,*)
$INCLUDE ThermalUnits
TABLE    pStorageUnits(g,*)
$INCLUDE StorageUnits
TABLE    pNetwork(n,n,c,*)
$INCLUDE Network
TABLE    pTransMatrix(s,ss)
$INCLUDE TransitionMatrix
TABLE    pActivePeriods(p,*)
$INCLUDE ActivePeriods
TABLE    pTransMatrix_rp(rp,rrpp)
$INCLUDE TransitionMatrix_rp
;

* Delete the loaded ranges from memory
EXECUTE 'del tmp.txt Demand Demand_s Duration Frequency Hindex HourChanges Indices' ;
EXECUTE 'del tmp.txt Inflows Inflows_s Network Param Set_gn StorageUnits ThermalUnits' ;
EXECUTE 'del tmp.txt TransitionMatrix Wind Wind_s Solar Solar_s ActivePeriods' ;
EXECUTE 'del tmp.txt ReprePeriods Weight Hindex_rp Set_gtech HourChangesRFM ' ;
EXECUTE 'del tmp.txt TransitionMatrix_rp ' ;

*-------------------------------------------------------------------------------
*                        SET VALUE VARIABLES AND DYNAMIC SETS
*-------------------------------------------------------------------------------
* Specify dynamic sets
t (g)$[    pThermalUnits (g,'FuelCost'  )] = YES ;
h (g)$[NOT pThermalUnits (g,'FuelCost'  )] = YES ;
hf(h)$[    pStorageUnits (h,'type'      )] = YES ; // activating "fast" storage technology
hs(h)$[NOT pStorageUnits (h,'type'      )] = YES ; // activating "slow" storage technology
pa(p)$[    pActivePeriods(p,'Active'    )] = YES ;

* Activating set hr depending on model execution
IF(pModel = 1 OR pModel = 2 OR pModel = 3 OR pModel = 4,
   hr(h)$[hf(h) OR hs(h)] = YES ; // storage units with reservoir are both "slow" and "fast"
);

IF(pModel = 5,
   hr(h)$[hf(h)         ] = YES ; // storage units with reservoir are only "fast", "slow" are limited by a new constraint in the model
);

* Activating set of periods to calculate storage level for "slow" units in RP-TM model

FOR(pScalar = 1 to CARD(p) by pStorMovWindow,
   ps(p) $[ORD(p) = pScalar] = yes;
);
   ps(p) $[ORD(p) = CARD(p)] = YES ; // also including the last period

* Parameter to indicate the last hour on ps(p) set
IF(Mod[CARD(p),pStorMovWindow]=0,
   pLastStorMovW = pStorMovWindow - 1 ;
ELSE
   pLastStorMovW = CARD(p) - Floor[CARD(p)/pStorMovWindow] * pStorMovWindow - 1 ;
);

* Loading information for thermal units
pAlfa    (t)   = pThermalUnits(t,'alfa'    ) ;
pBeta    (t)   = pThermalUnits(t,'beta'    ) ;
pGamma   (t)   = pThermalUnits(t,'gamma'   ) ;
pTheta   (t)   = pThermalUnits(t,'theta'   ) ;
pCostfuel(t)   = pThermalUnits(t,'FuelCost') ;
pCostOM  (t)   = pThermalUnits(t,'OMCost'  ) ;
pQmax    (t)   = pThermalUnits(t,'Qmax'    ) ;
pQmin    (t)   = pThermalUnits(t,'Qmin'    ) ;
pUo      (t)   = pThermalUnits(t,'Uo'      ) ;
pK       (t)   = pThermalUnits(t,'k'       ) ;
pSRR     (t)   = pThermalUnits(t,'SRR'     ) ;

* Loading information for storage units
pWmax    (h)   = pStorageUnits(h,'MaxReserve'      ) ;
pWmin    (h)   = pStorageUnits(h,'MinReserve'      ) ;
pW0      (h)   = pStorageUnits(h,'IniReserve'      ) ;
pWfin    (h)   = pStorageUnits(h,'FinReserve'      ) ;
pQmax    (h)   = pStorageUnits(h,'Qmax'            ) ;
pQmin    (h)   = pStorageUnits(h,'Qmin'            ) ;
pK       (h)   = pStorageUnits(h,'k'               ) ;
pBmax    (h)   = pStorageUnits(h,'Bmax'            ) ;
pEta     (h)   = pStorageUnits(h,'Efficiency'      ) ;
pInveCost(h)   = pStorageUnits(h,'AnnualInvestCost') ;
pMaxStorI(h)   = pStorageUnits(h,'MaxStorInvest'   ) ;
pEPRmax  (h)   = pStorageUnits(h,'EPRmax'          ) ;
pEPRmin  (h)   = pStorageUnits(h,'EPRmin'          ) ;

* Reduced Frequency Matrix calculation
pFreqMatRed (s,ss,ch) $[ORD(ch) = 1]
               = pFrequency(s,ss,ch)                         ;
pFreqMatRed (s,ss,ch) $[ORD(ch) > 1]
               = pFrequency(s,ss,ch) - pFrequency(s,ss,ch-1) ;

* Scaling Transition Matrix for Representative Periods
pMaxTM_rp = SMAX[(rp,rrpp), pTransMatrix_rp(rp,rrpp)]           ;
pTransMatrix_rp  (rp,rrpp)= pTransMatrix_rp(rp,rrpp) / pMaxTM_rp;

* Loading information for network
pTCmax   (n,m,c) = pNetwork(n,m,c,'TCmax') * 1e-3            ;
pXnm     (n,m,c)  $pNetwork(n,m,c,'X') = pNetwork(n,m,c,'X') ;

* Line conecttions
ns    (n    ) $ [ORD(n) <> pSlack]  = YES ;
Omega (n,m,c) $ pNetwork(n,m,c,'X') = YES ;
Omega (n,m,c) $ [NOT pTransNet]     = NO  ;
Omega2(n,m  )       = SUM[c,Omega(n,m,c)] ;

* Representative periods information
pNumPer_rp(rp)                       = SUM[p $[rpp(rp,p)],1]           ; // Periods per Representative
pWeight_rp(rp) $[    pNumPer_rp(rp)] = pWeight_rp(rp) / pNumPer_rp(rp) ; // Weight per period
pWeight_rp(rp) $[NOT pNumPer_rp(rp)] = 0                               ;

* Procedure to find the first period at each representative day
pFirstP_rp(rp,p) $[NOT rpp(rp,p)] = 0      ;                           // 1. Set to zero periods not within representatives periods
pFirstP_rp(rp,p) $[    rpp(rp,p)] = ORD(p) ;                           // 2. Assign value to periods within representatives periods
pFirstP_rp(rp,p)                  = + pFirstP_rp(rp,p  )               // 3. Rule to identify first period
                                    - pFirstP_rp(rp,p-1)$[ORD(p) > 1]
                                    + 1                 $[ORD(p) = 1];
pFirstP_rp(rp,p) $[NOT rpp(rp,p)] = 0 ;                                // 4. Set to zero periods not within representatives periods
pFirstP_rp(rp,p) $[    rpp(rp,p) AND pFirstP_rp(rp,p) = 1] = 0 ;       // 5. Only the first period of each rp is stored
pFirstP_rp(rp,p) $[    rpp(rp,p) AND           ORD(p) = 1] = 1 ;       // 6. Expecific rule for first period

* Definition of parameter indicating the number of periods in the RP minus 1
pNumPerMinus1_rp(rp) = pNumPer_rp(rp) - 1 ;                            // To model the startup constraint when RP use the Transition Matrix

*-------------------------------------------------------------------------------
*                              Data Validation
*-------------------------------------------------------------------------------
FILE WARN /warnings.out/ ;
PUT  WARN ;
PUT  'Warning file Storage Investment Model'/;

* Checking that Qmin production vs Inflows on each representative period
*-----------------------------------------------------------------------
LOOP((rp,h),
   IF(SUM[p $rpp(rp,p), pk(h) * pQmin(h)] > SUM[p $rpp(rp,p), pInflows(p,h)],
      PUT WARN ' Minimum production for ' h.TL:0 ' over the representative period ' rp.TL:0 ' is higher than inflows.' / ;
      PUT WARN ' Inflows are increased considering the minimum production needed' / ;
      LOOP(p $[rpp(rp,p) AND pNumPer_rp(rp)],
         pInflows(p,h) = + pInflows(p,h)
                         + [+ SUM[pp $rpp(rp,pp), pk(h)*pQmin(h)]
                            - SUM[pp $rpp(rp,pp), pInflows(pp,h)]
                           ] / pNumPer_rp(rp);
      );
   );
);
PUTCLOSE WARN ;
*-------------------------------------------------------------------------------
*                      Injection Sensitivity Factors
*-------------------------------------------------------------------------------
IF(pTransNet,
* Obtaining susceptance matrix with multiple circuits defined as n->m
    pYBUS(Omega2(n,m)) = -SUM[c $ (Omega(n,m,c) AND pXnm(n,m,c)),1/pXnm(n,m,c)]-SUM[c $ (Omega(n,m,c) AND pXnm (m,n,c)),1/pXnm(m,n,c)] ;

* Creating the symmetric matrix
    pYBUS(m,n) $Omega2(n,m) = pYBUS(n,m) ;

* Creating the diagonal for the B matriz
    pYBUS(m,m) = SUM(n,-pYBUS(n,m)) ;

*Eliminaing the Slack/reference node
    pYBUS(n,m) $[NOT(ns(n))] = 0 ;
    pYBUS(n,m) $[NOT(ns(m))] = 0 ;

* Creating gdx files for inverse calculation
    EXECUTE_UNLOAD 'gdxfrominverse.gdx' ns ;
    EXECUTE_UNLOAD 'gdxforinverse.gdx ' ns ;

*Obtaining the inverse of pYBUS and saving it in pYBUSInv
    EXECUTE_UNLOAD 'gdxforinverse.gdx' ns,pYBUS                            ;
    EXECUTE 'invert gdxforinverse.gdx ns pYBUS gdxfrominverse.gdx pYBUSInv';
    EXECUTE_LOAD 'gdxfrominverse.gdx' , pYBUSInv                           ;

* Obtaining the Injection Shift Factors
    pISF(n,m,c,ns) $ [pXnm(n,m,c) AND Omega(n,m,c)] = (pYBUSInv(n,ns) - pYBUSInv(m,ns))/pXnm(n,m,c) ;
* eliminating the pISF lower than 1E-9 (contribution lower than 1W)
    pISF(n,m,c,ns) $ (ABS(pISF(n,m,c,ns)) < 1E-9) = 0 ;
);
*-------------------------------------------------------------------------------
*                        SET BOUNDS FOR VARIABLES
*-------------------------------------------------------------------------------

* Storage Investment bounds
vStorageInvest.FX (h) $[NOT pStorageInvest] = 0            ; //Investment mode disabled
vStorageInvest.UP (h) $[    pStorageInvest] = pMaxStorI(h) ;

* Set bound for the power flow
vCircuitPowerFlow.UP   (pa,n,m,c) =  pTCmax(n,m,c) ;
vCircuitPowerFlow.LO   (pa,n,m,c) = -pTCmax(n,m,c) ;
vCircuitPowerFlow_s.UP (s ,n,m,c) =  pTCmax(n,m,c) ;
vCircuitPowerFlow_s.LO (s ,n,m,c) = -pTCmax(n,m,c) ;

* Set upper bound for wind  generation
vWind.UP   (pa,n) = pWind   (pa,n) ;
vWind_s.UP (s ,n) = pWind_s (s ,n) ;

* Set upper bound for solar generation
vSolar.UP  (pa,n) = pSolar  (pa,n) ;
vSolar_s.UP(s ,n) = pSolar_s(s ,n) ;

* Set limits of the hydro variables
vProduction.UP      (pa,h) $[NOT pStorageInvest] = pk   (h) *  pQmax  (h)                 ;
vProduction.LO      (pa,h) $[NOT pStorageInvest] = pk   (h) *  pQmin  (h)                 ;
vProduction_s.UP    (s ,h) $[NOT pStorageInvest] = pk   (h) *  pQmax  (h)                 ;
vProduction_s.LO    (s ,h) $[NOT pStorageInvest] = pk   (h) *  pQmin  (h)                 ;
vConsumption.UP     (p ,h) $[NOT pStorageInvest] = pBmax(h)                               ;
vConsumption_s.UP   (s ,h) $[NOT pStorageInvest] = pBmax(h)                               ;
vStorageReserve.UP  (p ,h) $[NOT pStorageInvest] = pWmax(h)                               ;
vStorageReserve.LO  (p ,h) $[NOT pStorageInvest] = pWmin(h)                               ;
vProduction.UP      (pa,h) $[    pStorageInvest] = pk   (h) * [pQmax  (h) + pMaxStorI(h)] ;
vProduction.LO      (pa,h) $[    pStorageInvest] = pk   (h) *  pQmin  (h)                 ;
vProduction_s.UP    (s ,h) $[    pStorageInvest] = pk   (h) * [pQmax  (h) + pMaxStorI(h)] ;
vProduction_s.LO    (s ,h) $[    pStorageInvest] = pk   (h) *  pQmin  (h)                 ;
vConsumption.UP     (p ,h) $[    pStorageInvest] = pBmax(h) +  pEta   (h) * pMaxStorI(h)  ;
vConsumption_s.UP   (s ,h) $[    pStorageInvest] = pBmax(h) +  pEta   (h) * pMaxStorI(h)  ;
vStorageReserve.UP  (p ,h) $[    pStorageInvest] = pWmax(h) +  pEPRmax(h) * pMaxStorI(h)  ;
vStorageReserve.LO  (p ,h) $[    pStorageInvest] = pWmin(h) +  pEPRmin(h) * pMaxStorI(h)  ;

vStorageReserve.LO  (pa(p),h)$[ORD(p) = CARD(p)]=pWfin(h) ;

* Set the difference in storage level to zero for thermal generators and if no transition exists
vDeltaStoReserve.FX(s,ss,g)$[NOT h(g)] = 0             ;
vDeltaStoReserve.FX(s,ss,g)$[pTransMatrix(s,ss)=0] = 0 ;

* Upper bounds for Start-up (hourly and system states)
vStartUp.UP  (pa  ,t) = 1 ;
vStartUp_s.UP(s,ss,t) = 1 ;

* Set start-ups to zero in first period. Only possible in hourly model
vStartUp.FX(pa(p),t)$[ORD(p) = 1] = 0 ;

* If there is no transition between two states, there can be no start-ups between them either.
vStartUp_s.FX(s,ss,t)$[pTransMatrix(s,ss)=0] = 0 ;

* Set power not supply upper bound equal to the demand
vPowerNotSupply.UP  (pa,n) = pDemand  (pa,n)   ;
vPowerNotSupply_s.UP(s ,n) = pDemand_s(s ,n)   ;
vPowerNotSupply.FX  (pa,n) $[pENSflag = 0]= 0 ;
vPowerNotSupply_s.FX(s ,n) $[pENSflag = 0]= 0 ;

*-------------------------------------------------------------------------------
*                                MODEL SOLVES
*-------------------------------------------------------------------------------

* ------------------------------ hourly model ----------------------------------
IF(pModel = 1,
    SOLVE MSEM
          USING MIP MINIMIZING vObjectiveFunction;
);

* ------------------------- system states model --------------------------------
IF(pModel = 2,
    SOLVE MSEM_S
          USING MIP MINIMIZING vObjectiveFunction;
);

* ---------------- system states model reduced F matrix-------------------------
IF(pModel = 3,
    SOLVE MSEM_SS_RFM
          USING MIP MINIMIZING vObjectiveFunction;
);

* --------------------- representative periods model ---------------------------
IF(pModel = 4,



    SOLVE MSEM_RP
          USING MIP MINIMIZING vObjectiveFunction;
);

* ---------- representative periods model with Transition Matrix ---------------
IF(pModel = 5,
    pa(p) = NO ;                      // 1. reset active periods set
    pa(p) $[SUM[rpp(rp,p),1]] = YES ; // 2. set active periods equal to representatives
    SOLVE MSEM_RP_TM
          USING MIP MINIMIZING vObjectiveFunction;
);

*-------------------------------------------------------------------------------
*                       CALCULATE EX POST PARAMETERS
*-------------------------------------------------------------------------------
* Parameters for output in hourly model
PARAMETERS
OF_Cost    (*  )    Total cost model                   [k€]
GenCPUTime (*  )    Generation CPU time                [ s]
SolCPUTime (*  )    Solve      CPU time                [ s]
NumVar     (*  )    Number of variables
NumDVar    (*  )    Number of discrete variables
NumEqu     (*  )    Number of equations
NumNZ      (*  )    Number of nonzero entries in the model coefficient matrix
BestSol    (*  )    The estimate of the best possible solution for a MIP
StorInvest (g,*)    Storage Investment                 [MW ]
Qhourly    (p,g)    Hourly production                  [GWh]
Bhourly    (p,g)    Hourly pumping                     [GWh]
Shourly    (p,g)    Hourly spillage                    [GWh]
Uhourly    (p,g)    Hourly conection                   [0-1]
Whourly    (p,g)    Hourly storage level               [GWh]
Windhourly (p,n)    Hourly wind                        [GWh]
Solarhourly(p,n)    Hourly solar                       [GWh]
WindSpill  (p,n)    Hourly wind  spillage              [GWh]
SolarSpill (p,n)    Hourly solar spillage              [GWh]
PNShourly  (p,n)    Hourly pns                         [GWh]
PRICEhourly(p,*,*)  Hourly price                 [€ per MWh]
Fhourly  (p,n,n,c)  Hourly power flow                  [GWh]
Tot_tech (p,tech,*) Total  production per technology   [GWh]
TotTechSU(  tech,*) Total startups per technology      [ # ]
MovWindow1 (p  )    Moving Window auxiliar parameter
MovWindow2 (p  )    Moving Window auxiliar parameter
;

IF(pModel = 1,
    OF_Cost    ('Obj Func  Model      [1000 M€]') = vObjectiveFunction.L + EPS ;

    GenCPUTime ('CPU Time  Model generation [s]') = MSEM.resGen  ;
    SolCPUTime ('CPU Time  Model solution   [s]') = MSEM.resUsd  ;
    NumVar     ('Number of variables           ') = MSEM.numVar  ;
    NumDVar    ('Number of discrete variables  ') = MSEM.numDVar ;
    NumEqu     ('Number of equations           ') = MSEM.numEqu  ;
    NumNZ      ('Number of nonzero elements    ') = MSEM.numNZ   ;
    BestSol    ('Best possible solution for MIP') = MSEM.objest  ;

    Qhourly    (pa,g) = vProduction.L    (pa,g)  +EPS ;
    Uhourly    (pa,t) = vConnected.L     (pa,t)  +EPS ;
    Bhourly    (pa,h) = vConsumption.L   (pa,h)  +EPS ;
    Shourly    (pa,h) = vSpillage.L      (pa,h)  +EPS ;
    Whourly    (pa,h) = vStorageReserve.L(pa,h)  +EPS ;
    Windhourly (pa,n) = vWind.L          (pa,n)  +EPS ;
    Solarhourly(pa,n) = vSolar.L         (pa,n)  +EPS ;
    PNShourly  (pa,n) = vPowerNotSupply.L(pa,n)  +EPS ;
    Fhourly  (pa,n,m,c) $Omega(n,m,c) = vCircuitPowerFlow.L(pa,n,m,c)         + EPS ;

    WindSpill  (pa,n) = pWind (pa,n) - vWind.L (pa,n)  +EPS ;
    SolarSpill (pa,n) = pSolar(pa,n) - vSolar.L(pa,n)  +EPS ;

    Tot_tech(pa,tech,'[GWh]') = SUM[g$gtech(g,tech), vProduction.L(pa,g)]     + EPS ;
    TotTechSU(  tech,'[# startups]') = SUM[(pa,t)$gtech(t,tech),
                                             vStartUp.L(pa,t  )]              + EPS ;

    PRICEhourly(pa,n     ,'[€/MWh]') $[    pTransNet] = 1e6 * E_DMNDN.M(pa,n) + EPS ;
    PRICEhourly(pa,'BUS1','[€/MWh]') $[NOT pTransNet] = 1e6 * E_DMND.M (pa  ) + EPS ;

    StorInvest (h,'[MW]') =  vStorageInvest.L (h) * 1e3 + EPS ;
);

IF(pModel = 2,
    OF_Cost    ('Obj Func  Model      [1000 M€]') = vObjectiveFunction.L + EPS ;

    GenCPUTime ('CPU Time  Model generation [s]') = MSEM_S.resGen  ;
    SolCPUTime ('CPU Time  Model solution   [s]') = MSEM_S.resUsd  ;
    NumVar     ('Number of variables           ') = MSEM_S.numVar  ;
    NumDVar    ('Number of discrete variables  ') = MSEM_S.numDVar ;
    NumEqu     ('Number of equations           ') = MSEM_S.numEqu  ;
    numNZ      ('Number of nonzero elements    ') = MSEM_S.numNZ   ;
    BestSol    ('Best possible solution for MIP') = MSEM_S.objest  ;

    Qhourly    (pa,g) = SUM[hindex(pa,s),vProduction_s.L    (s,g)]  + EPS ;
    Uhourly    (pa,t) = SUM[hindex(pa,s),vConnected_s.L     (s,t)]  + EPS ;
    Bhourly    (pa,h) = SUM[hindex(pa,s),vConsumption_s.L   (s,h)]  + EPS ;
    Shourly    (pa,h) = SUM[hindex(pa,s),vSpillage_s.L      (s,h)]  + EPS ;
    Windhourly (pa,n) = SUM[hindex(pa,s),vWind_s.L          (s,n)]  + EPS ;
    Solarhourly(pa,n) = SUM[hindex(pa,s),vSolar_s.L         (s,n)]  + EPS ;
    PNShourly  (pa,n) = SUM[hindex(pa,s),vPowerNotSupply_s.L(s,n)]  + EPS ;

    WindSpill  (pa,n) = SUM[hindex(pa,s),pWind_s (s,n) - vWind_s.L (s,n)]  + EPS ;
    SolarSpill (pa,n) = SUM[hindex(pa,s),pSolar_s(s,n) - vSolar_s.L(s,n)]  + EPS ;

    Whourly    (pa(p),h)  = pW0(h) + SUM[pp$[ORD(pp)<=ORD(p)],
                                        -           Qhourly (pp,h)
                                        -           Shourly (pp,h)
                                        + pEta(h) * Bhourly (pp,h)
                                        +           pInflows(pp,h)] + EPS ;

    Fhourly    (pa,n,m,c) $Omega(n,m,c) = SUM[hindex(pa,s),
                                              vCircuitPowerFlow_s.L(s,n,m,c)]     + EPS ;

    Tot_tech  (pa,tech,'[GWh]') = SUM[(s,g)$hindex(pa,s),
                                              vProduction_s.L(s,g)$gtech(g,tech)] + EPS ;
    TotTechSU (   tech,'[# startups]') = SUM[(s,t), SUM[ss$(ORD(ss) <> ORD(s)),
                                              vStartUp_s.L(ss,s,t)$gtech(t,tech) * pTransMatrix(ss,s)]]         + EPS ;

    PRICEhourly(pa,n     ,'[€/MWh]') $[    pTransNet] = SUM[hindex(pa,s),1e6 * E_DMNDN_S.M(s,n)/pDuration_s(s)] + EPS ;
    PRICEhourly(pa,'BUS1','[€/MWh]') $[NOT pTransNet] = SUM[hindex(pa,s),1e6 * E_DMND_S.M (s  )/pDuration_s(s)] + EPS ;

    StorInvest (h,'[MW]') =  vStorageInvest.L (h) * 1e3 + EPS ;
);

IF(pModel = 3,
    OF_Cost    ('Obj Func  Model      [1000 M€]') = vObjectiveFunction.L + EPS ;

    GenCPUTime ('CPU Time  Model generation [s]') = MSEM_SS_RFM.resGen  ;
    SolCPUTime ('CPU Time  Model solution   [s]') = MSEM_SS_RFM.resUsd  ;
    NumVar     ('Number of variables           ') = MSEM_SS_RFM.numVar  ;
    NumDVar    ('Number of discrete variables  ') = MSEM_SS_RFM.numDVar ;
    NumEqu     ('Number of equations           ') = MSEM_SS_RFM.numEqu  ;
    NumNZ      ('Number of nonzero elements    ') = MSEM_SS_RFM.numNZ   ;
    BestSol    ('Best possible solution for MIP') = MSEM_SS_RFM.objest  ;

    Qhourly    (pa,g) = SUM[hindex(pa,s),vProduction_s.L    (s,g)] + EPS ;
    Uhourly    (pa,t) = SUM[hindex(pa,s),vConnected_s.L     (s,t)] + EPS ;
    Bhourly    (pa,h) = SUM[hindex(pa,s),vConsumption_s.L   (s,h)] + EPS ;
    Shourly    (pa,h) = SUM[hindex(pa,s),vSpillage_s.L      (s,h)] + EPS ;
    Windhourly (pa,n) = SUM[hindex(pa,s),vWind_s.L          (s,n)] + EPS ;
    Solarhourly(pa,n) = SUM[hindex(pa,s),vSolar_s.L         (s,n)] + EPS ;
    PNShourly  (pa,n) = SUM[hindex(pa,s),vPowerNotSupply_s.L(s,n)] + EPS ;

    WindSpill  (pa,n) = SUM[hindex(pa,s),pWind_s (s,n) - vWind_s.L (s,n)]  + EPS ;
    SolarSpill (pa,n) = SUM[hindex(pa,s),pSolar_s(s,n) - vSolar_s.L(s,n)]  + EPS ;

    // Calculation of Moving Window for System States using a the Reduced Frequency Matrix
    LOOP(ch,
        MovWindow1(p) $[ORD(p)= 1                             ] = 1                                    ; // 1. The first hour has to be equal to 1
        MovWindow1(p) $[ORD(p)= pHourChanges(ch) AND ORD(ch)=1] = pHourChanges(ch)                     ; // 2. The first change keep the same value of hour
        MovWindow1(p) $[ORD(p)= pHourChanges(ch)              ] = pHourChanges(ch) - pHourChanges(ch-1); // 3. We calcule the difference in order to obtain the number of hours between the changes or hours in which the bounds will be imposed
    );

    MovWindow2(p) = 1; // 4.We need a second parameter with all elements equal to 1 at the beginning
    LOOP(pp,
        MovWindow2(p) $[MovWindow1(p) = 0 AND ORD(pp)<ORD(p)]  // 5. this loop and if condition ensure that always the count restart at the hours that have been selected for bounds
                      = MovWindow2(p-1) + 1 ;
    );
    // end of Moving Window calculation

    Whourly (pa(p),hf(h))  = pW0(h) + SUM[pp$[ORD(pp)>=ORD(p)+1-MovWindow2(p)
                                          AND ORD(pp)<=ORD(p)],
                                        -           Qhourly (pp,h)
                                        -           Shourly (pp,h)
                                        + pEta(h) * Bhourly (pp,h)
                                        +           pInflows(pp,h)] + EPS ;

    Whourly (pa(p),hs(h))  = pW0(h) + SUM[pp$[ORD(pp)<=ORD(p)],
                                        -           Qhourly (pp,h)
                                        -           Shourly (pp,h)
                                        + pEta(h) * Bhourly (pp,h)
                                        +           pInflows(pp,h)] + EPS ;

    Fhourly    (pa,n,m,c) $Omega(n,m,c) = SUM[hindex(pa,s),
                                              vCircuitPowerFlow_s.L(s,n,m,c)]     + EPS ;

    Tot_tech  (pa,tech,'[GWh]') = SUM[(s,g)$hindex(pa,s),
                                              vProduction_s.L(s,g)$gtech(g,tech)] + EPS ;

    TotTechSU (   tech,'[# startups]') = SUM[(s,t), SUM[ss$(ORD(ss) <> ORD(s)),
                                              vStartUp_s.L(ss,s,t)$gtech(t,tech) * pTransMatrix(ss,s)]]         + EPS ;

    PRICEhourly(pa,n     ,'[€/MWh]') $[    pTransNet] = SUM[hindex(pa,s),1e6 * E_DMNDN_S.M(s,n)/pDuration_s(s)] + EPS ;
    PRICEhourly(pa,'BUS1','[€/MWh]') $[NOT pTransNet] = SUM[hindex(pa,s),1e6 * E_DMND_S.M (s  )/pDuration_s(s)] + EPS ;

    StorInvest (h,'[MW]') =  vStorageInvest.L (h) * 1e3 + EPS ;
);

IF(pModel = 4,
    OF_Cost    ('Obj Func  Model      [1000 M€]') = vObjectiveFunction.L + EPS ;

    GenCPUTime ('CPU Time  Model generation [s]') = MSEM_RP.resGen  ;
    SolCPUTime ('CPU Time  Model solution   [s]') = MSEM_RP.resUsd  ;
    NumVar     ('Number of variables           ') = MSEM_RP.numVar  ;
    NumDVar    ('Number of discrete variables  ') = MSEM_RP.numDVar ;
    NumEqu     ('Number of equations           ') = MSEM_RP.numEqu  ;
    NumNZ      ('Number of nonzero elements    ') = MSEM_RP.numNZ   ;
    BestSol    ('Best possible solution for MIP') = MSEM_RP.objest  ;

    Qhourly    (p,g) = SUM[hindexRP(p,pp),vProduction.L    (pp,g)] + EPS ;
    Uhourly    (p,t) = SUM[hindexRP(p,pp),vConnected.L     (pp,t)] + EPS ;
    Bhourly    (p,h) = SUM[hindexRP(p,pp),vConsumption.L   (pp,h)] + EPS ;
    Shourly    (p,h) = SUM[hindexRP(p,pp),vSpillage.L      (pp,h)] + EPS ;
    Windhourly (p,n) = SUM[hindexRP(p,pp),vWind.L          (pp,n)] + EPS ;
    Solarhourly(p,n) = SUM[hindexRP(p,pp),vSolar.L         (pp,n)] + EPS ;
    PNShourly  (p,n) = SUM[hindexRP(p,pp),vPowerNotSupply.L(pp,n)] + EPS ;
    Whourly    (p,h) = SUM[hindexRP(p,pp),vStorageReserve.L(pp,h)] + EPS ;

    WindSpill  (p,n) = SUM[hindexRP(p,pp),pWind (pp,n) - vWind.L (pp,n)]  + EPS ;
    SolarSpill (p,n) = SUM[hindexRP(p,pp),pSolar(pp,n) - vSolar.L(pp,n)]  + EPS ;

    Fhourly(p,n,m,c) $Omega(n,m,c) = SUM[hindexRP(p,pp),
                                         vCircuitPowerFlow.L(pp,n,m,c)]     + EPS ;

    Tot_tech(p,tech,'[GWh]') = SUM[(pp,g)$hindexRP(p,pp),
                                         vProduction.L(pp,g)$gtech(g,tech)] + EPS ;
    TotTechSU (   tech,'[# startups]') = SUM[rpp(rp,pa(p)),
                                         pWeight_rp(rp) * SUM[t, vStartUp.L(p,t)$gtech(t,tech)]] + EPS ;

    PRICEhourly(p,n     ,'[€/MWh]') $[    pTransNet] = SUM[hindexRP(p,pp),1e6 * E_DMNDN.M(pp,n) / SUM[rp $rpp(rp,pp), pWeight_rp(rp)]] + EPS ;
    PRICEhourly(p,'BUS1','[€/MWh]') $[NOT pTransNet] = SUM[hindexRP(p,pp),1e6 * E_DMND.M (pp  ) / SUM[rp $rpp(rp,pp), pWeight_rp(rp)]] + EPS ;

    StorInvest (h,'[MW]') =  vStorageInvest.L (h) * 1e3 + EPS ;
);

IF(pModel = 5,
    OF_Cost    ('Obj Func  Model      [1000 M€]') = vObjectiveFunction.L + EPS ;

    GenCPUTime ('CPU Time  Model generation [s]') = MSEM_RP_TM.resGen  ;
    SolCPUTime ('CPU Time  Model solution   [s]') = MSEM_RP_TM.resUsd  ;
    NumVar     ('Number of variables           ') = MSEM_RP_TM.numVar  ;
    NumDVar    ('Number of discrete variables  ') = MSEM_RP_TM.numDVar ;
    NumEqu     ('Number of equations           ') = MSEM_RP_TM.numEqu  ;
    NumNZ      ('Number of nonzero elements    ') = MSEM_RP_TM.numNZ   ;
    BestSol    ('Best possible solution for MIP') = MSEM_RP_TM.objest  ;

    Qhourly    (p,g) = SUM[hindexRP(p,pp),vProduction.L    (pp,g)] + EPS ;
    Uhourly    (p,t) = SUM[hindexRP(p,pp),vConnected.L     (pp,t)] + EPS ;
    Bhourly    (p,h) = SUM[hindexRP(p,pp),vConsumption.L   (pp,h)] + EPS ;
    Shourly    (p,h) = SUM[hindexRP(p,pp),vSpillage.L      (pp,h)] + EPS ;
    Windhourly (p,n) = SUM[hindexRP(p,pp),vWind.L          (pp,n)] + EPS ;
    Solarhourly(p,n) = SUM[hindexRP(p,pp),vSolar.L         (pp,n)] + EPS ;
    PNShourly  (p,n) = SUM[hindexRP(p,pp),vPowerNotSupply.L(pp,n)] + EPS ;

    WindSpill  (p,n) = SUM[hindexRP(p,pp),pWind (pp,n) - vWind.L (pp,n)]  + EPS ;
    SolarSpill (p,n) = SUM[hindexRP(p,pp),pSolar(pp,n) - vSolar.L(pp,n)]  + EPS ;

    Whourly    (p,hf(h)) = SUM[hindexRP(p,pp),vStorageReserve.L(pp,h)] + EPS ;
    Whourly    (p,hs(h)) =                    vStorageReserve.L(p ,h)  + EPS ;

    Fhourly(p,n,m,c) $Omega(n,m,c) = SUM[hindexRP(p,pp),
                                         vCircuitPowerFlow.L(pp,n,m,c)]     + EPS ;

    Tot_tech(p,tech,'[GWh]') = SUM[(pp,g)$hindexRP(p,pp),
                                         vProduction.L(pp,g)$gtech(g,tech)] + EPS ;
    TotTechSU (   tech,'[# startups]') = SUM[rpp(rp,pa(p)),
                                         pWeight_rp(rp) * SUM[t, vStartUp.L(p,t)$gtech(t,tech)]] + EPS ;

    PRICEhourly(p,n     ,'[€/MWh]') $[    pTransNet] = SUM[hindexRP(p,pp),1e6 * E_DMNDN.M(pp,n) / SUM[rp $rpp(rp,pp), pWeight_rp(rp)]] + EPS ;
    PRICEhourly(p,'BUS1','[€/MWh]') $[NOT pTransNet] = SUM[hindexRP(p,pp),1e6 * E_DMND.M (pp  ) / SUM[rp $rpp(rp,pp), pWeight_rp(rp)]] + EPS ;

    StorInvest (h,'[MW]') =  vStorageInvest.L (h) * 1e3 + EPS ;
);

* Data output to xls file
EXECUTE 'del tmp.xlsx' // Deleting tmp file

IF(pModel = 1,
    PUT TMP PUT 'par=Whourly rdim=1 rng=Whourly!a1' / 'par=Qhourly rdim=1 rng=Qhourly!a1' / 'par=Bhourly rdim=1 rng=Bhourly!a1' / 'par=Shourly rdim=1 rng=Shourly!a1' /
    PUT TMP PUT 'par=Uhourly rdim=1 rng=Uhourly!a1' / 'par=Tot_tech rdim=1 rng=Thermal!a2' / 'par=Fhourly rdim=1 rng=Fhourly!a1' / 'par=OF_Cost rdim=1 rng=OperCost_hourly!a1' /
    PUT TMP PUT 'par=PNShourly rdim=1 rng=PNShourly!a1' / 'par=PRICEhourly rdim=1 rng=PRICEhourly!a1' / 'par=Windhourly rdim=1 rng=Windhourly!a1' /
    PUT TMP PUT 'par=GenCPUTime rdim=1 rng=OperCost_hourly!a2' / 'par=SolCPUTime rdim=1 rng=OperCost_hourly!a3' / 'par=NumVar rdim=1 rng=OperCost_hourly!a4' /
    PUT TMP PUT 'par=NumDVar rdim=1 rng=OperCost_hourly!a5' / 'par=NumEqu rdim=1 rng=OperCost_hourly!a6' / 'par=NumNZ rdim=1 rng=OperCost_hourly!a7' /
    PUT TMP PUT 'par=Solarhourly rdim=1 rng=Solarhourly!a1' / 'par=BestSol rdim=1 rng=OperCost_hourly!a8'/ 'par=TotTechSU rdim=1 rng=TotTechSU_hourly!a1' /
    PUT TMP PUT 'par=WindSpill rdim=1 rng=WindSpill_hourly!a1' / 'par=SolarSpill rdim=1 rng=SolarSpill_hourly!a1' / 'par=StorInvest rdim=1 rng=StorInvest_hourly!a1'
    PUTCLOSE
    EXECUTE_UNLOAD   'tmp.gdx' Whourly Qhourly Bhourly Shourly Uhourly PNShourly Fhourly Tot_tech OF_Cost PRICEhourly Windhourly GenCPUTime SolCPUTime NumVar NumDVar NumEqu NumNZ Solarhourly BestSol TotTechSU WindSpill SolarSpill StorInvest
    EXECUTE          'gdxxrw.exe tmp.gdx SQ=n EpsOut=0 O="tmp.xlsx" @tmp.txt'
    EXECUTE          'del        tmp.gdx                             tmp.txt'
);

IF(pModel = 2,
    PUT TMP PUT 'par=Fhourly rdim=1 rng=Fs_hourly!a1' / 'par=Qhourly rdim=1 rng=Qs_hourly!a1' / 'par=Bhourly rdim=1 rng=Bs_hourly!a1' / 'par=Shourly rdim=1 rng=Ss_hourly!a1' /
    PUT TMP PUT 'par=Uhourly rdim=1 rng=Us_hourly!a1' / 'par=Whourly rdim=1 rng=Ws_hourly!a1' / 'par=Tot_tech rdim=1 rng=Thermal_s!a2' / 'par=OF_Cost rdim=1 rng=OperCost_SS!a1' /
    PUT TMP PUT 'par=PNShourly rdim=1 rng=PNSs_hourly!a1' / 'par=PRICEhourly rdim=1 rng=PRICEs_hourly!a1' / 'par=Windhourly rdim=1 rng=Winds_hourly!a1' /
    PUT TMP PUT 'par=GenCPUTime rdim=1 rng=OperCost_SS!a2' / 'par=SolCPUTime rdim=1 rng=OperCost_SS!a3' / 'par=NumVar rdim=1 rng=OperCost_SS!a4' /
    PUT TMP PUT 'par=NumDVar rdim=1 rng=OperCost_SS!a5' / 'par=NumEqu rdim=1 rng=OperCost_SS!a6' / 'par=NumNZ rdim=1 rng=OperCost_SS!a7' /
    PUT TMP PUT 'par=Solarhourly rdim=1 rng=Solars_hourly!a1' / 'par=BestSol rdim=1 rng=OperCost_SS!a8' / 'par=TotTechSU rdim=1 rng=TotTechSU_SS!a1' /
    PUT TMP PUT 'par=WindSpill rdim=1 rng=WindSpills_hourly!a1' / 'par=SolarSpill rdim=1 rng=SolarSpills_hourly!a1' / 'par=StorInvest rdim=1 rng=StorInvest_SS!a1'
    PUTCLOSE
    EXECUTE_UNLOAD   'tmp.gdx' Whourly Qhourly Bhourly Shourly Uhourly PNShourly Fhourly Tot_tech OF_Cost PRICEhourly Windhourly GenCPUTime SolCPUTime NumVar NumDVar NumEqu NumNZ Solarhourly BestSol TotTechSU WindSpill SolarSpill StorInvest
    EXECUTE          'gdxxrw.exe tmp.gdx SQ=n EpsOut=0 O="tmp.xlsx" @tmp.txt'
    EXECUTE          'del        tmp.gdx                             tmp.txt'
);

IF(pModel = 3,
    PUT TMP PUT 'par=Fhourly rdim=1 rng=Fs_RFM!a1' / 'par=Qhourly rdim=1 rng=Qs_RFM!a1' / 'par=Bhourly rdim=1 rng=Bs_RFM!a1' / 'par=Shourly rdim=1 rng=Ss_RFM!a1' /
    PUT TMP PUT 'par=Uhourly rdim=1 rng=Us_RFM!a1' / 'par=Whourly rdim=1 rng=Ws_RFM!a1' / 'par=Tot_tech rdim=1 rng=Thermal_s_RFM!a2' / 'par=OF_Cost rdim=1 rng=OperCost_SS-RFM!a1' /
    PUT TMP PUT 'par=PNShourly rdim=1 rng=PNSs_RFM!a1' / 'par=PRICEhourly rdim=1 rng=PRICEs_RFM!a1' / 'par=Windhourly rdim=1 rng=Winds_RFM!a1' /
    PUT TMP PUT 'par=GenCPUTime rdim=1 rng=OperCost_SS-RFM!a2' / 'par=SolCPUTime rdim=1 rng=OperCost_SS-RFM!a3' / 'par=NumVar rdim=1 rng=OperCost_SS-RFM!a4' /
    PUT TMP PUT 'par=NumDVar rdim=1 rng=OperCost_SS-RFM!a5' / 'par=NumEqu rdim=1 rng=OperCost_SS-RFM!a6' / 'par=NumNZ rdim=1 rng=OperCost_SS-RFM!a7' /
    PUT TMP PUT 'par=Solarhourly rdim=1 rng=Solars_RFM!a1' / 'par=BestSol rdim=1 rng=OperCost_SS-RFM!a8' / 'par=TotTechSU rdim=1 rng=TotTechSU_SS-RFM!a1' /
    PUT TMP PUT 'par=WindSpill rdim=1 rng=WindSpill_RFM!a1' / 'par=SolarSpill rdim=1 rng=SolarSpill_RFM!a1' / 'par=StorInvest rdim=1 rng=StorInvest_SS-RFM!a1'
    PUTCLOSE
    EXECUTE_UNLOAD   'tmp.gdx' Whourly Qhourly Bhourly Shourly Uhourly PNShourly Fhourly Tot_tech OF_Cost PRICEhourly Windhourly GenCPUTime SolCPUTime NumVar NumDVar NumEqu NumNZ Solarhourly BestSol TotTechSU WindSpill SolarSpill StorInvest
    EXECUTE          'gdxxrw.exe tmp.gdx SQ=n EpsOut=0 O="tmp.xlsx" @tmp.txt'
    EXECUTE          'del        tmp.gdx                             tmp.txt'
);

IF(pModel = 4,
    PUT TMP PUT 'par=Whourly rdim=1 rng=W_RP!a1' / 'par=Qhourly rdim=1 rng=Q_RP!a1' / 'par=Bhourly rdim=1 rng=B_RP!a1' / 'par=Shourly rdim=1 rng=S_RP!a1' /
    PUT TMP PUT 'par=Uhourly rdim=1 rng=U_RP!a1' / 'par=Tot_tech rdim=1 rng=Thermal_RP!a2' / 'par=Fhourly rdim=1 rng=F_RP!a1' / 'par=OF_Cost rdim=1 rng=OperCost_RP!a1' /
    PUT TMP PUT 'par=PNShourly rdim=1 rng=PNS_RP!a1' / 'par=PRICEhourly rdim=1 rng=PRICE_RP!a1' / 'par=Windhourly rdim=1 rng=Wind_RP!a1' /
    PUT TMP PUT 'par=GenCPUTime rdim=1 rng=OperCost_RP!a2' / 'par=SolCPUTime rdim=1 rng=OperCost_RP!a3' / 'par=NumVar rdim=1 rng=OperCost_RP!a4' /
    PUT TMP PUT 'par=NumDVar rdim=1 rng=OperCost_RP!a5' / 'par=NumEqu rdim=1 rng=OperCost_RP!a6' / 'par=NumNZ rdim=1 rng=OperCost_RP!a7' /
    PUT TMP PUT 'par=Solarhourly rdim=1 rng=Solar_RP!a1' / 'par=BestSol rdim=1 rng=OperCost_RP!a8' / 'par=TotTechSU rdim=1 rng=TotTechSU_RP!a1' /
    PUT TMP PUT 'par=WindSpill rdim=1 rng=WindSpill_RP!a1' / 'par=SolarSpill rdim=1 rng=SolarSpill_RP!a1' / 'par=StorInvest rdim=1 rng=StorInvest_RP!a1'
    PUTCLOSE
    EXECUTE_UNLOAD   'tmp.gdx' Whourly Qhourly Bhourly Shourly Uhourly PNShourly Fhourly Tot_tech OF_Cost PRICEhourly Windhourly GenCPUTime SolCPUTime NumVar NumDVar NumEqu NumNZ Solarhourly BestSol TotTechSU WindSpill SolarSpill StorInvest
    EXECUTE          'gdxxrw.exe tmp.gdx SQ=n EpsOut=0 O="tmp.xlsx" @tmp.txt'
    EXECUTE          'del        tmp.gdx                             tmp.txt'
);

IF(pModel = 5,
    PUT TMP PUT 'par=Whourly rdim=1 rng=W_RP-TM!a1' / 'par=Qhourly rdim=1 rng=Q_RP-TM!a1' / 'par=Bhourly rdim=1 rng=B_RP-TM!a1' / 'par=Shourly rdim=1 rng=S_RP-TM!a1' /
    PUT TMP PUT 'par=Uhourly rdim=1 rng=U_RP-TM!a1' / 'par=Tot_tech rdim=1 rng=Thermal_RP-TM!a2' / 'par=Fhourly rdim=1 rng=F_RP-TM!a1' / 'par=OF_Cost rdim=1 rng=OperCost_RP-TM!a1' /
    PUT TMP PUT 'par=PNShourly rdim=1 rng=PNS_RP-TM!a1' / 'par=PRICEhourly rdim=1 rng=PRICE_RP-TM!a1' / 'par=Windhourly rdim=1 rng=Wind_RP-TM!a1' /
    PUT TMP PUT 'par=GenCPUTime rdim=1 rng=OperCost_RP-TM!a2' / 'par=SolCPUTime rdim=1 rng=OperCost_RP-TM!a3' / 'par=NumVar rdim=1 rng=OperCost_RP-TM!a4' /
    PUT TMP PUT 'par=NumDVar rdim=1 rng=OperCost_RP-TM!a5' / 'par=NumEqu rdim=1 rng=OperCost_RP-TM!a6' / 'par=NumNZ rdim=1 rng=OperCost_RP-TM!a7' /
    PUT TMP PUT 'par=Solarhourly rdim=1 rng=Solar_RP-TM!a1' / 'par=BestSol rdim=1 rng=OperCost_RP-TM!a8' / 'par=TotTechSU rdim=1 rng=TotTechSU_RP-TM!a1' /
    PUT TMP PUT 'par=WindSpill rdim=1 rng=WindSpill_RP-TM!a1' / 'par=SolarSpill rdim=1 rng=SolarSpill_RP-TM!a1' / 'par=StorInvest rdim=1 rng=StorInvest_RP-TM!a1'
    PUTCLOSE
    EXECUTE_UNLOAD   'tmp.gdx' Whourly Qhourly Bhourly Shourly Uhourly PNShourly Fhourly Tot_tech OF_Cost PRICEhourly Windhourly GenCPUTime SolCPUTime NumVar NumDVar NumEqu NumNZ Solarhourly BestSol TotTechSU WindSpill SolarSpill StorInvest
    EXECUTE          'gdxxrw.exe tmp.gdx SQ=n EpsOut=0 O="tmp.xlsx" @tmp.txt'
    EXECUTE          'del        tmp.gdx                             tmp.txt'
);

display pNumPer_rp, pa ;
