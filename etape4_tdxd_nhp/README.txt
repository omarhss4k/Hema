étape 4 — T-DXd NHP (singe cynomolgus)
=======================================

Cette étape contient la simulation PK/PD de T-DXd chez le singe cynomolgus,
analogue à l'étape 3 (rat) et l'étape 5 (humain).

À implémenter :
  - parameters_nhp.R        : baselines PD NHP (Fornari, allométrie depuis rat)
  - parameters_tdxd_nhp.R   : PK T-DXd NHP (TMDD, allométrie depuis humain)
  - pkpd_tdxd_nhp.R         : ODE fusionné T-DXd + Fornari PD
  - run_pkpd_tdxd_nhp.R     : simulations scénarios FDA BLA Table 7
  - plots_tdxd_nhp.R        : figures publication-ready

Note : les données FGFR2 NHP se trouvent dans etape6_fgfr2/
