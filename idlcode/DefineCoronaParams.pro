function DefineCoronaParams, Tbase, nbase, Q0, a, b, $
                             force_isothermal=force_isothermal, interpolB=interpolB, AddTR=AddTR
 mode=0L
 if keyword_set(force_isothermal) then mode=mode or 1L
 if keyword_set(interpolB) then mode=mode or 2L                            
 if keyword_set(AddTR) then mode=mode or 4L
                             
 coronaparms={Tbase: double(Tbase), $
              nbase: double(nbase), $
              Q0: double(Q0), $
              a: double(a), $
              b: double(b), $
              mode: mode}
              
 return, coronaparms
end