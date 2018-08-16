c  FLTIMER = whether to do timing
c  FLAVER  = whether to do axial average
c  FLTHRM  = whether to include temperature
c  FLSTLE  = whether tt has settling velocity
c  FLGVARY = whether gravity vector changes with y 
c  FLPGVAR = whether pressure gradient changes with time 
      module flags_m
      LOGICAL FLTIMER,FLAVER,FLTHRM,FLSTLE,
     *        FLSV_YES,FLSV_NO,FLGVARY,FLPGVAR
      COMMON /flags/ FLTIMER,FLAVER,FLTHRM,
     *               FLSTLE,
     *               FLSV_YES,FLSV_NO,
     *               FLGVARY,FLPGVAR  
      end module
