1

   Abaqus 2020                                  Date 05-Dec-2022   Time 13:10:22
   For use by UNIVERSITY OF CAMBRIDGE under license from Dassault Systemes or its subsidiary.

                                                                                 
                                                                                 
 STEP    1     INCREMENT     1     STEP TIME    0.00    


                        S T E P       1     S T A T I C   A N A L Y S I S


                                                                                          

     AUTOMATIC TIME CONTROL WITH -
          A SUGGESTED INITIAL TIME INCREMENT OF                 1.00    
          AND A TOTAL TIME PERIOD OF                            1.00    
          THE MINIMUM TIME INCREMENT ALLOWED IS                1.000E-05
          THE MAXIMUM TIME INCREMENT ALLOWED IS                 1.00    

     LINEAR EQUATION SOLVER TYPE         DIRECT SPARSE

 CONVERGENCE TOLERANCE PARAMETERS FOR FORCE    
     CRITERION FOR RESIDUAL FORCE     FOR A NONLINEAR PROBLEM          5.000E-03
     CRITERION FOR DISP.    CORRECTION IN A NONLINEAR PROBLEM          1.000E-02
     INITIAL VALUE OF TIME AVERAGE FORCE                               1.000E-02
     AVERAGE FORCE     IS TIME AVERAGE FORCE    
     ALTERNATE CRIT. FOR RESIDUAL FORCE     FOR A NONLINEAR PROBLEM    2.000E-02
     CRITERION FOR ZERO FORCE     RELATIVE TO TIME AVRG. FORCE         1.000E-05
     CRITERION FOR RESIDUAL FORCE     WHEN THERE IS ZERO FLUX          1.000E-05
     CRITERION FOR DISP.    CORRECTION WHEN THERE IS ZERO FLUX         1.000E-03
     CRITERION FOR RESIDUAL FORCE     FOR A LINEAR INCREMENT           1.000E-08
     FIELD CONVERSION RATIO                                             1.00    
     CRITERION FOR ZERO FORCE     REL. TO TIME AVRG. MAX. FORCE        1.000E-05
     CRITERION FOR ZERO DISP.    RELATIVE TO CHARACTERISTIC LENGTH     1.000E-08

     VOLUMETRIC STRAIN COMPATIBILITY TOLERANCE FOR HYBRID SOLIDS       1.000E-05
     AXIAL STRAIN COMPATIBILITY TOLERANCE FOR HYBRID BEAMS             1.000E-05
     TRANS. SHEAR STRAIN COMPATIBILITY TOLERANCE FOR HYBRID BEAMS      1.000E-05
     SOFT CONTACT CONSTRAINT COMPATIBILITY TOLERANCE FOR P>P0          5.000E-03
     SOFT CONTACT CONSTRAINT COMPATIBILITY TOLERANCE FOR P=0.0         0.100    
     CONTACT FORCE ERROR TOLERANCE FOR CONVERT SDI=YES                 1.00    
     DISPLACEMENT COMPATIBILITY TOLERANCE FOR DCOUP ELEMENTS           1.000E-05
     ROTATION COMPATIBILITY TOLERANCE FOR DCOUP ELEMENTS               1.000E-05

 EQUILIBRIUM WILL BE CHECKED FOR SEVERE DISCONTINUITY ITERATIONS

 TIME INCREMENTATION CONTROL PARAMETERS:
     FIRST EQUILIBRIUM ITERATION FOR CONSECUTIVE DIVERGENCE CHECK              4
     EQUILIBRIUM ITERATION AT WHICH LOG. CONVERGENCE RATE CHECK BEGINS         8
     EQUILIBRIUM ITERATION AFTER WHICH ALTERNATE RESIDUAL IS USED              9
     MAXIMUM EQUILIBRIUM ITERATIONS ALLOWED                                   16
     EQUILIBRIUM ITERATION COUNT FOR CUT-BACK IN NEXT INCREMENT               10
     MAXIMUM EQUILIB. ITERS IN TWO INCREMENTS FOR TIME INCREMENT INCREASE      4
     MAXIMUM ITERATIONS FOR SEVERE DISCONTINUITIES                            50
     MAXIMUM ATTEMPTS ALLOWED IN AN INCREMENT                                  5
     MAXIMUM DISCON. ITERS IN TWO INCREMENTS FOR TIME INCREMENT INCREASE      50
     MAXIMUM CONTACT AUGMENTATIONS FOR *SURFACE BEHAVIOR,AUGMENTED LAGRANGE   50
     CUT-BACK FACTOR AFTER DIVERGENCE                                 0.2500    
     CUT-BACK FACTOR FOR TOO SLOW CONVERGENCE                         0.5000    
     CUT-BACK FACTOR AFTER TOO MANY EQUILIBRIUM ITERATIONS            0.7500    
     CUT-BACK FACTOR AFTER TOO MANY SEVERE DISCONTINUITY ITERATIONS   0.2500    
     CUT-BACK FACTOR AFTER PROBLEMS IN ELEMENT ASSEMBLY               0.2500    
     INCREASE FACTOR AFTER TWO INCREMENTS THAT CONVERGE QUICKLY        1.500    
     MAX. TIME INCREMENT INCREASE FACTOR ALLOWED                       1.500    
     MAX. TIME INCREMENT INCREASE FACTOR ALLOWED (DYNAMICS)            1.250    
     MAX. TIME INCREMENT INCREASE FACTOR ALLOWED (DIFFUSION)           2.000    
     MINIMUM TIME INCREMENT RATIO FOR EXTRAPOLATION TO OCCUR          0.1000    
     MAX. RATIO OF TIME INCREMENT TO STABILITY LIMIT                   1.000    
     FRACTION OF STABILITY LIMIT FOR NEW TIME INCREMENT               0.9500    
     TIME INCREMENT INCREASE FACTOR BEFORE A TIME POINT                1.000    
     GLOBAL STABILIZATION CONTROL IS NOT USED

          PRINT OF INCREMENT NUMBER, TIME, ETC., EVERY    1  INCREMENTS

     THE MAXIMUM NUMBER OF INCREMENTS IN THIS STEP IS                     100

          LARGE DISPLACEMENT THEORY WILL BE USED

     LINEAR EXTRAPOLATION WILL BE USED

     CHARACTERISTIC ELEMENT LENGTH     5.000E-02

     DETAILS REGARDING ACTUAL SOLUTION WAVEFRONT REQUESTED

     DETAILED OUTPUT OF DIAGNOSTICS TO DATABASE REQUESTED

     PRINT OF INCREMENT NUMBER, TIME, ETC., TO THE MESSAGE FILE EVERY     1  INCREMENTS

     COLLECTING MODEL CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS

     COLLECTING STEP CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS


  INCREMENT     1 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT   1.00    

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.10000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  0.10000    
       ELAPSED WALLCLOCK TIME (SEC) =          1

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      5.880E-02   TIME AVG. FORCE       5.880E-02
 LARGEST RESIDUAL FORCE             -15.2       AT NODE         72   DOF  3
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         0.300       AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.        0.300       AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 32 NEGATIVE EIGENVALUES.

     NEGATIVE EIGENVALUES MEAN THAT THE SYSTEM MATRIX IS NOT 
     POSITIVE DEFINITE:
     FOR EXAMPLE, A BIFURCATION (BUCKLING) LOAD MAY HAVE BEEN EXCEEDED.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      0.144       TIME AVG. FORCE       0.144    
 LARGEST RESIDUAL FORCE             -57.0       AT NODE        718   DOF  3
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         0.300       AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.        0.142       AT NODE        727   DOF  1
   INSTANCE: PART-1                                                                          
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 54 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.
   ELEMENT PART-1.6 IS DISTORTING EXCESSIVELY.
   ELEMENT PART-1.8 IS DISTORTING EXCESSIVELY.
   ELEMENT PART-1.16 IS DISTORTING EXCESSIVELY.
   ELEMENT PART-1.56 IS DISTORTING EXCESSIVELY.
   ELEMENT PART-1.62 IS DISTORTING EXCESSIVELY.
   ELEMENT PART-1.64 IS DISTORTING EXCESSIVELY.
   ELEMENT PART-1.71 IS DISTORTING EXCESSIVELY.
   ELEMENT PART-1.72 IS DISTORTING EXCESSIVELY.
   ELEMENT PART-1.80 IS DISTORTING EXCESSIVELY.
   ELEMENT PART-1.120 IS DISTORTING EXCESSIVELY.
   ELEMENT PART-1.127 IS DISTORTING EXCESSIVELY.
   ELEMENT PART-1.128 IS DISTORTING EXCESSIVELY.

 ***WARNING: EXCESSIVE DISTORTION AT A TOTAL OF 12 INTEGRATION POINTS IN SOLID 
             (CONTINUUM) ELEMENTS
 

 ***NOTE: ELEMENTS ARE DISTORTING EXCESSIVELY. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT     1 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  0.250    

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.730E-03   TIME AVG. FORCE       4.730E-03
 LARGEST RESIDUAL FORCE            -0.454       AT NODE        720   DOF  3
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         7.500E-02   AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.        7.500E-02   AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.496E-03   TIME AVG. FORCE       4.496E-03
 LARGEST RESIDUAL FORCE             0.198       AT NODE        720   DOF  3
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         7.500E-02   AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.       -1.921E-02   AT NODE        720   DOF  3
   INSTANCE: PART-1                                                                          
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      4.463E-03   TIME AVG. FORCE       4.463E-03
 LARGEST RESIDUAL FORCE             5.966E-02   AT NODE        720   DOF  3
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         7.500E-02   AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.       -7.588E-03   AT NODE        729   DOF  3
   INSTANCE: PART-1                                                                          
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.10000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  0.10000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      4.474E-03   TIME AVG. FORCE       4.474E-03
 LARGEST RESIDUAL FORCE             3.294E-03   AT NODE        720   DOF  3
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         7.500E-02   AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.        3.546E-03   AT NODE        720   DOF  3
   INSTANCE: PART-1                                                                          
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      4.477E-03   TIME AVG. FORCE       4.477E-03
 LARGEST RESIDUAL FORCE             3.483E-05   AT NODE        720   DOF  3
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         7.500E-02   AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.       -2.029E-04   AT NODE        729   DOF  3
   INSTANCE: PART-1                                                                          
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      4.477E-03   TIME AVG. FORCE       4.477E-03
 LARGEST RESIDUAL FORCE             1.011E-09   AT NODE        720   DOF  3
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         7.500E-02   AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.        1.974E-06   AT NODE        720   DOF  3
   INSTANCE: PART-1                                                                          
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.250    ,  FRACTION OF STEP COMPLETED  0.250    
 STEP TIME COMPLETED       0.250    ,  TOTAL TIME COMPLETED        0.250    


  INCREMENT     2 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.250    

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 8 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      1.449E-02   TIME AVG. FORCE       9.486E-03
 LARGEST RESIDUAL FORCE             0.252       AT NODE        153   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         7.500E-02   AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.        3.344E-02   AT NODE         90   DOF  2
   INSTANCE: PART-1                                                                          
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      1.442E-02   TIME AVG. FORCE       9.448E-03
 LARGEST RESIDUAL FORCE            -5.726E-02   AT NODE        585   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         7.500E-02   AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.       -1.467E-02   AT NODE        576   DOF  2
   INSTANCE: PART-1                                                                          
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      1.419E-02   TIME AVG. FORCE       9.335E-03
 LARGEST RESIDUAL FORCE            -1.637E-02   AT NODE        666   DOF  3
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         7.500E-02   AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.       -1.203E-02   AT NODE        720   DOF  2
   INSTANCE: PART-1                                                                          
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      1.413E-02   TIME AVG. FORCE       9.304E-03
 LARGEST RESIDUAL FORCE            -9.014E-04   AT NODE        638   DOF  3
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         7.500E-02   AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.       -2.405E-03   AT NODE        576   DOF  2
   INSTANCE: PART-1                                                                          
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      1.414E-02   TIME AVG. FORCE       9.308E-03
 LARGEST RESIDUAL FORCE             3.876E-06   AT NODE        720   DOF  3
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         7.500E-02   AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.       -1.381E-04   AT NODE        648   DOF  1
   INSTANCE: PART-1                                                                          
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   5 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  5 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.250    ,  FRACTION OF STEP COMPLETED  0.500    
 STEP TIME COMPLETED       0.500    ,  TOTAL TIME COMPLETED        0.500    


  INCREMENT     3 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.250    

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 10 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      2.923E-02   TIME AVG. FORCE       1.595E-02
 LARGEST RESIDUAL FORCE             0.298       AT NODE        584   DOF  1
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         7.500E-02   AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.        5.634E-02   AT NODE        576   DOF  2
   INSTANCE: PART-1                                                                          
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      2.881E-02   TIME AVG. FORCE       1.581E-02
 LARGEST RESIDUAL FORCE             0.160       AT NODE        638   DOF  1
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         7.500E-02   AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.       -3.473E-02   AT NODE        648   DOF  2
   INSTANCE: PART-1                                                                          
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      2.814E-02   TIME AVG. FORCE       1.558E-02
 LARGEST RESIDUAL FORCE            -0.201       AT NODE        666   DOF  3
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         7.500E-02   AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.       -3.921E-02   AT NODE        576   DOF  2
   INSTANCE: PART-1                                                                          
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      2.780E-02   TIME AVG. FORCE       1.547E-02
 LARGEST RESIDUAL FORCE             1.720E-02   AT NODE        584   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         7.500E-02   AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.       -9.620E-03   AT NODE         90   DOF  2
   INSTANCE: PART-1                                                                          
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      2.777E-02   TIME AVG. FORCE       1.546E-02
 LARGEST RESIDUAL FORCE             3.119E-03   AT NODE        584   DOF  3
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         7.500E-02   AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.       -5.220E-03   AT NODE        666   DOF  2
   INSTANCE: PART-1                                                                          
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      2.778E-02   TIME AVG. FORCE       1.546E-02
 LARGEST RESIDUAL FORCE             7.012E-06   AT NODE        584   DOF  3
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         7.500E-02   AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.       -1.467E-04   AT NODE        720   DOF  2
   INSTANCE: PART-1                                                                          
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   6 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  6 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  0.250    ,  FRACTION OF STEP COMPLETED  0.750    
 STEP TIME COMPLETED       0.750    ,  TOTAL TIME COMPLETED        0.750    


  INCREMENT     4 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  0.250    

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 2 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      6.535E-02   TIME AVG. FORCE       2.794E-02
 LARGEST RESIDUAL FORCE             -12.7       AT NODE        160   DOF  1
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         0.114       AT NODE         79   DOF  1
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.        9.351E-02   AT NODE         79   DOF  1
   INSTANCE: PART-1                                                                          
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 20 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.
   ELEMENT PART-1.5 IS DISTORTING EXCESSIVELY.
   ELEMENT PART-1.7 IS DISTORTING EXCESSIVELY.
   ELEMENT PART-1.61 IS DISTORTING EXCESSIVELY.
   ELEMENT PART-1.63 IS DISTORTING EXCESSIVELY.

 ***WARNING: EXCESSIVE DISTORTION AT A TOTAL OF 4 INTEGRATION POINTS IN SOLID 
             (CONTINUUM) ELEMENTS
 

 ***NOTE: ELEMENTS ARE DISTORTING EXCESSIVELY. CONVERGENCE IS JUDGED UNLIKELY.


  INCREMENT     4 STARTS. ATTEMPT NUMBER  2, TIME INCREMENT  6.250E-02

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      3.267E-02   TIME AVG. FORCE       1.976E-02
 LARGEST RESIDUAL FORCE            -7.783E-02   AT NODE        638   DOF  3
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         1.875E-02   AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.        1.331E-02   AT NODE        576   DOF  2
   INSTANCE: PART-1                                                                          
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.274E-02   TIME AVG. FORCE       1.978E-02
 LARGEST RESIDUAL FORCE             1.229E-02   AT NODE        648   DOF  3
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         1.875E-02   AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.        9.134E-03   AT NODE        720   DOF  2
   INSTANCE: PART-1                                                                          
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.277E-02   TIME AVG. FORCE       1.979E-02
 LARGEST RESIDUAL FORCE            -5.454E-04   AT NODE        638   DOF  3
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         1.875E-02   AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.        1.708E-03   AT NODE        720   DOF  2
   INSTANCE: PART-1                                                                          
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      3.277E-02   TIME AVG. FORCE       1.979E-02
 LARGEST RESIDUAL FORCE             1.815E-06   AT NODE        720   DOF  3
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         1.875E-02   AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.        8.817E-05   AT NODE        666   DOF  2
   INSTANCE: PART-1                                                                          
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   4 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  4 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  6.250E-02,  FRACTION OF STEP COMPLETED  0.812    
 STEP TIME COMPLETED       0.812    ,  TOTAL TIME COMPLETED        0.812    


  INCREMENT     5 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  6.250E-02

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0
 

 ***WARNING: THE SYSTEM MATRIX HAS 6 NEGATIVE EIGENVALUES.
             EXPLANATIONS ARE SUGGESTED AFTER THE FIRST OCCURRENCE OF THIS MESSAGE.

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.009E-02   TIME AVG. FORCE       2.385E-02
 LARGEST RESIDUAL FORCE             -2.04       AT NODE        152   DOF  1
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.        -7.238E-02   AT NODE        162   DOF  1
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.       -7.520E-02   AT NODE        162   DOF  1
   INSTANCE: PART-1                                                                          
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =  0.10000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =  0.10000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      3.762E-02   TIME AVG. FORCE       2.336E-02
 LARGEST RESIDUAL FORCE            -0.251       AT NODE        152   DOF  1
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.        -5.438E-02   AT NODE        162   DOF  1
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.       -2.367E-02   AT NODE        576   DOF  2
   INSTANCE: PART-1                                                                          
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      3.768E-02   TIME AVG. FORCE       2.337E-02
 LARGEST RESIDUAL FORCE            -0.361       AT NODE        638   DOF  3
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         4.362E-02   AT NODE         72   DOF  3
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.       -4.404E-02   AT NODE        720   DOF  2
   INSTANCE: PART-1                                                                          
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =  0.10000    
       ELAPSED TOTAL CPU TIME (SEC) =  0.10000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     4


 AVERAGE FORCE                      3.775E-02   TIME AVG. FORCE       2.338E-02
 LARGEST RESIDUAL FORCE             9.827E-02   AT NODE        162   DOF  1
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.        -2.984E-02   AT NODE        162   DOF  1
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.       -2.315E-02   AT NODE         72   DOF  3
   INSTANCE: PART-1                                                                          
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     5


 AVERAGE FORCE                      3.777E-02   TIME AVG. FORCE       2.339E-02
 LARGEST RESIDUAL FORCE             2.454E-02   AT NODE         90   DOF  1
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.        -2.091E-02   AT NODE        162   DOF  1
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.        8.923E-03   AT NODE         90   DOF  1
   INSTANCE: PART-1                                                                          
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     6


 AVERAGE FORCE                      3.782E-02   TIME AVG. FORCE       2.340E-02
 LARGEST RESIDUAL FORCE             7.010E-02   AT NODE         72   DOF  3
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         1.875E-02   AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.       -1.759E-02   AT NODE         72   DOF  3
   INSTANCE: PART-1                                                                          
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     7


 AVERAGE FORCE                      3.776E-02   TIME AVG. FORCE       2.338E-02
 LARGEST RESIDUAL FORCE             4.557E-03   AT NODE        162   DOF  1
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         1.875E-02   AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.        4.417E-03   AT NODE         72   DOF  1
   INSTANCE: PART-1                                                                          
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     8


 AVERAGE FORCE                      3.773E-02   TIME AVG. FORCE       2.338E-02
 LARGEST RESIDUAL FORCE             4.514E-03   AT NODE         72   DOF  3
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         1.875E-02   AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.        4.658E-03   AT NODE        162   DOF  1
   INSTANCE: PART-1                                                                          
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     9


 AVERAGE FORCE                      3.774E-02   TIME AVG. FORCE       2.338E-02
 LARGEST RESIDUAL FORCE             6.798E-06   AT NODE        162   DOF  3
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         1.875E-02   AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.       -1.749E-04   AT NODE         79   DOF  3
   INSTANCE: PART-1                                                                          
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   9 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  9 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  6.250E-02,  FRACTION OF STEP COMPLETED  0.875    
 STEP TIME COMPLETED       0.875    ,  TOTAL TIME COMPLETED        0.875    


  INCREMENT     6 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  6.250E-02

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.218E-02   TIME AVG. FORCE       2.651E-02
 LARGEST RESIDUAL FORCE            -1.682E-02   AT NODE        638   DOF  3
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         1.875E-02   AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.        4.802E-03   AT NODE        648   DOF  3
   INSTANCE: PART-1                                                                          
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.222E-02   TIME AVG. FORCE       2.652E-02
 LARGEST RESIDUAL FORCE             7.986E-04   AT NODE        584   DOF  3
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         1.875E-02   AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.        2.321E-03   AT NODE        666   DOF  2
   INSTANCE: PART-1                                                                          
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     3


 AVERAGE FORCE                      4.223E-02   TIME AVG. FORCE       2.652E-02
 LARGEST RESIDUAL FORCE             6.692E-06   AT NODE         72   DOF  3
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         1.875E-02   AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.       -1.578E-04   AT NODE        162   DOF  1
   INSTANCE: PART-1                                                                          
          THE FORCE     EQUILIBRIUM EQUATIONS HAVE CONVERGED

 ITERATION SUMMARY FOR THE INCREMENT:   3 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  3 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  6.250E-02,  FRACTION OF STEP COMPLETED  0.938    
 STEP TIME COMPLETED       0.938    ,  TOTAL TIME COMPLETED        0.938    


  INCREMENT     7 STARTS. ATTEMPT NUMBER  1, TIME INCREMENT  6.250E-02

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     1


 AVERAGE FORCE                      4.712E-02   TIME AVG. FORCE       2.946E-02
 LARGEST RESIDUAL FORCE            -1.155E-03   AT NODE        638   DOF  3
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         1.875E-02   AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.        1.603E-03   AT NODE         72   DOF  1
   INSTANCE: PART-1                                                                          
          FORCE     EQUILIBRIUM NOT ACHIEVED WITHIN TOLERANCE.

  NUMBER OF EQUATIONS =      2187     NUMBER OF RHS =    1
  NUMBER OF FLOATING PT. OPERATIONS =  8.65E+07

     CHECK POINT   START OF SOLVER    

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   0.0000    
       ELAPSED SYSTEM TIME (SEC)    =   0.0000    
       ELAPSED TOTAL CPU TIME (SEC) =   0.0000    
       ELAPSED WALLCLOCK TIME (SEC) =          0

               CONVERGENCE CHECKS FOR EQUILIBRIUM ITERATION     2


 AVERAGE FORCE                      4.713E-02   TIME AVG. FORCE       2.947E-02
 LARGEST RESIDUAL FORCE             1.260E-05   AT NODE        720   DOF  3
   INSTANCE: PART-1                                                                          
 LARGEST INCREMENT OF DISP.         1.875E-02   AT NODE          9   DOF  2
   INSTANCE: PART-1                                                                          
 LARGEST CORRECTION TO DISP.        2.271E-04   AT NODE        648   DOF  1
   INSTANCE: PART-1                                                                          
 ESTIMATE OF DISP.    CORRECTION   -2.477E-06
 FORCE     EQUILIB. ACCEPTED BASED ON SMALL RESIDUAL AND ESTIMATED CORRECTION

 ITERATION SUMMARY FOR THE INCREMENT:   2 TOTAL ITERATIONS, OF WHICH
   0 ARE SEVERE DISCONTINUITY ITERATIONS AND  2 ARE EQUILIBRIUM ITERATIONS.

 TIME INCREMENT COMPLETED  6.250E-02,  FRACTION OF STEP COMPLETED   1.00    
 STEP TIME COMPLETED        1.00    ,  TOTAL TIME COMPLETED         1.00    


          THE ANALYSIS HAS BEEN COMPLETED



     ANALYSIS SUMMARY:
     TOTAL OF          7  INCREMENTS
                       2  CUTBACKS IN AUTOMATIC INCREMENTATION
                      40  ITERATIONS INCLUDING CONTACT ITERATIONS IF PRESENT
                      40  PASSES THROUGH THE EQUATION SOLVER OF WHICH 
                      40  INVOLVE MATRIX DECOMPOSITION, INCLUDING
                       0  DECOMPOSITION(S) OF THE MASS MATRIX
                       1  REORDERING OF EQUATIONS TO MINIMIZE WAVEFRONT
                       0  ADDITIONAL RESIDUAL EVALUATIONS FOR LINE SEARCHES
                       0  ADDITIONAL OPERATOR EVALUATIONS FOR LINE SEARCHES
                       1  WARNING MESSAGES DURING USER INPUT PROCESSING
                       9  WARNING MESSAGES DURING ANALYSIS
                       0  ANALYSIS WARNINGS ARE NUMERICAL PROBLEM MESSAGES
                       7  ANALYSIS WARNINGS ARE NEGATIVE EIGENVALUE MESSAGES
                       0  ERROR MESSAGES



     JOB TIME SUMMARY
       USER TIME (SEC)      =   1.5000    
       SYSTEM TIME (SEC)    =  0.30000    
       TOTAL CPU TIME (SEC) =   1.8000    
       WALLCLOCK TIME (SEC) =          3
