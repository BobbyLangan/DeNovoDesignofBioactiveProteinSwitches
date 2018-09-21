import sys
import numpy as np
import sympy

def residual( CLc, CLo, K, T, CLK, CLT, CLKT, CLtot, Ktot, Ttot, K_o, K_ck, K_lt ):

    if( Ktot > 0 and Ttot > 0 ):
        one = CLc + CLo + CLK + CLKT + CLT - CLtot
        two = K + CLK + CLKT - Ktot
        three = T + CLKT + CLT - Ttot
        four =  CLo / CLc - K_o
        five = ( (K*CLo) / CLK ) - K_ck
        six =  ( (CLK*T) / CLKT ) - K_lt
        seven = ( (CLo * T) / CLT ) - K_lt
        eight = ( (CLT * K) / CLKT ) - K_ck
        return [one, two, three, four, five, six, seven, eight]

    elif( Ktot > 0 ):
        one = CLc + CLo + CLK - CLtot
        two = K + CLK - Ktot
        three =  CLo / CLc - K_o
        four = ( (K*CLo) / CLK ) - K_ck
        return [one, two, three, four]

    elif( Ttot > 0 ):
        one = CLc + CLo + CLT - CLtot
        two = T + CLKT + CLT - Ttot
        three =  CLo / CLc - K_o
        four = ( (CLo * T) / CLT ) - K_lt
        return [one, two, three, four]

    else:
        return 1000

#Returns concentration of all species in solution and residual
def model( CLtot, Ktot, Ttot, K_o, K_ck, K_lt ):
    #Unknowns
    CLc = sympy.Symbol('CLc')
    CLo = sympy.Symbol('CLo')
    K = sympy.Symbol('K')
    T = sympy.Symbol('T')
    CLK = sympy.Symbol('CLK')
    CLT = sympy.Symbol('CLT')
    CLKT = sympy.Symbol('CLKT')
    
    #To determine
    fracT_bound = 0
    fracK_bound = 0
    fracCL_bound_to_key = 0
    fracCL_bound_to_T = 0
    r = 0
    
    CLc_conc = 0
    CLo_conc = 0
    K_conc = 0
    T_conc = 0
    CLK_conc = 0
    CLT_conc = 0
    CLKT_conc = 0

    #If all parts present
    if( Ktot > 0 and Ttot > 0 ):
        #print( "Fitting eight equations" )
        one = CLc + CLo + CLK + CLKT + CLT - CLtot
        two = K + CLK + CLKT - Ktot
        three = T + CLKT + CLT - Ttot
        four =  CLo / CLc - K_o
        five = ( (K*CLo) / CLK ) - K_ck
        six =  ( (CLK*T) / CLKT ) - K_lt
        seven = ( (CLo * T) / CLT ) - K_lt
        eight = ( (CLT * K) / CLKT ) - K_ck
        
        equations = (one, two, three, four, five, six, seven, eight)
        symbols = (CLc, CLo, K, T, CLK, CLT, CLKT)#, CLKT_alt)
        guess = [CLtot, CLtot * K_o, Ktot, Ttot, .1, 0.01, Ttot/2]
        conc_bounds = [ (0,0,0,0,0,0,0), (CLtot, CLtot, Ktot, Ttot, min(CLtot, Ktot), min(CLtot, Ttot), min(CLtot, Ktot, Ttot)) ]
    
        try:
            solution = sympy.nsolve( equations, symbols, guess, bounds=conc_bounds, solver='bisect' )
        except:
            solution = sympy.nsolve( [func.as_numer_denom()[0] for func in equations], symbols, guess, verify=False, bounds=conc_bounds, solver='bisect' )

        CLc_conc = sympy.N( solution[0] )
        CLo_conc = sympy.N( solution[1] )
        K_conc = sympy.N( solution[2] )
        T_conc = sympy.N( solution[3] )
        CLK_conc = sympy.N( solution[4] )
        CLT_conc = sympy.N( solution[5] )
        CLKT_conc = sympy.N( solution[6] )

        r = sum( residual(CLc_conc, CLo_conc, K_conc, T_conc, CLK_conc, CLT_conc, CLKT_conc, CLtot, Ktot, Ttot, K_o, K_ck, K_lt) )
        fracT_bound = (CLKT_conc + CLT_conc) / Ttot
        fracK_bound = (CLK_conc + CLKT_conc) / Ktot
        fracCL_bound_to_key = CLK_conc / CLtot
        fracCL_bound_to_T = (CLKT_conc + CLT_conc) / CLtot

    #If only key is present
    elif( Ktot > 0 ):
        #print( "No target present, fitting to cage/key equillibrium" )
        one = CLc + CLo + CLK - CLtot
        two = K + CLK - Ktot
        three =  CLo / CLc - K_o
        four = ( (K*CLo) / CLK ) - K_ck
        
        equations = (one, two, three, four)
        symbols = (CLc, CLo, K, CLK)
        guess = [CLtot, CLtot * K_o, Ktot, 0.1]

        try:
            solution = sympy.nsolve( equations, symbols, guess )
        except:
            solution = sympy.nsolve( [func.as_numer_denom()[0] for func in equations], symbols, guess, verify=False )
            
        CLc_conc = sympy.N( solution[0] )
        CLo_conc = sympy.N( solution[1] )
        K_conc = sympy.N( solution[2] )
        CLK_conc = sympy.N( solution[3] )
        
        r = sum( residual(CLc_conc, CLo_conc, K_conc, T_conc, CLK_conc, CLT_conc, CLKT_conc, CLtot, Ktot, Ttot, K_o, K_ck, K_lt) )
        fracT_bound = 0
        fracK_bound = (CLK_conc + CLKT_conc) / Ktot
        fracCL_bound_to_key = CLK_conc / CLtot
        fracCL_bound_to_T = 0

    #If only target is present
    elif( Ttot > 0 ):
        #print( "No key present, fitting to cage/target equillibrium" )
        one = CLc + CLo + CLT - CLtot
        two = T + CLT - Ttot
        three =  CLo / CLc - K_o
        four = ( (CLo * T) / CLT ) - K_lt

        equations = (one, two, three, four)
        symbols = (CLc, CLo, T, CLT)
        guess = [CLtot, CLtot * K_o, Ttot, 0.001]

        try:
            solution = sympy.nsolve( equations, symbols, guess )
        except:
            solution = sympy.nsolve( [func.as_numer_denom()[0] for func in equations], symbols, guess, verify=False )
            
        CLc_conc = sympy.N( solution[0] )
        CLo_conc = sympy.N( solution[1] )
        T_conc = sympy.N( solution[2] )
        CLT_conc = sympy.N( solution[3] )

        r = sum( residual(CLc_conc, CLo_conc, K_conc, T_conc, CLK_conc, CLT_conc, CLKT_conc, CLtot, Ktot, Ttot, K_o, K_ck, K_lt) )
        fracT_bound = (CLKT_conc + CLT_conc) / Ttot
        fracK_bound = 0
        fracCL_bound_to_key = 0
        fracCL_bound_to_T = (CLKT_conc + CLT_conc) / CLtot
    
    #If its just the switch
    else:
        print ("Why are you doing this, there is no key nor target present")

    return (CLc_conc, CLo_conc, K_conc, T_conc, CLK_conc, CLT_conc, CLKT_conc, r)