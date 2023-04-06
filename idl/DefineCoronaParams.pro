function DefineCoronaParams, Tbase, nbase, Q0, a, b, force_isothermal=force_isothermal
 coronaparms={Tbase: double(Tbase), $
              nbase: double(nbase), $
              Q0: double(Q0), $
              a: double(a), $
              b: double(b), $
              force_isothermal: keyword_set(force_isothermal) ? 1L : 0L}
              
 return, coronaparms
end