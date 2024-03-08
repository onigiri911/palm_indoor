{salsa+simple.spc
Current revision
----------------
 20180903 Added SALSA variables               monakurppa
Former revisions
----------------
 $Id: smog.spc 2459 2017-09-13 14:10:33Z forkel $
}
#include atoms

 #DEFVAR
    O       = O ;      	    {oxygen atomic ground state (3P)}
    O3	    = 3O ;          {ozone}  
    NO	    = N + O ;       {nitric oxide}  
    NO2	    = N + 2O ;      {nitrogen dioxide} 
    NO3	    = N + 3O ;      {nitrogen trioxide} 
    N2O5    = 2N + 5O ;     {dinitrogen pentoxide}           
    HNO3    = H + N + 3O ;  { nitric acid }
    HNO4    = H + N + 4O ;  {HO2NO2 pernitric acid}
    H       = H ;           {hydrogen atomic ground state (2S)}
    OH	    = O + H ;       {hydroxyl radical}  
    HO2	    = H + 2O ;      {perhydroxyl radical}                     
    H2O2    = 2H + 2O ;     {hydrogen peroxide} 
    CH3     = C + 3H ;      {methyl radical}
    CH3O    = C + 3H + O ;  {methoxy radical}
    CH3O2   = C + 3H + 2O ; {methylperoxy radical}
    CH3OOH  = C + 4H + 2O ; {CH4O2 methylperoxy alcohol}
    HCO     = H + C + O ;   {CHO formyl radical}
    CH2O    = C + 2H + O ;  {formalydehyde}
    CO      = C + O ;       {carbon monoxide}   
    
    RH      = ignore ;      {alkanes} 
    RO2     = ignore ;      {alkyl peroxy radical}
    RCHO    = ignore ;      {carbonyl}
    RCOO2   = ignore ;      
    RCOO2NO2= ignore ;

    H2SO4   = 2H + S +4O ;  {sulfuric acid}
    NH3     = 3H + N ;      {ammonia}
    OCNV    = ignore ;      {non-volatile OC}
    OCSV    = ignore ;      {semi-volatile OC}

#DEFFIX
    H2O     = H + 2O ;      {water}
    H2      = 2H ;          {molecular hydrogen}
    O2      = 2O ;          {molecular oxygen}
    N2      = 2N ;          {molecular nitrogen}
    CH4     = C + 4H ;      {methane}
    CO2     = C + 2O ;      {carbon dioxide}

