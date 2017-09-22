// função PMV low
var PMVlow=function(CLO,MET,WME,TA,TR,VEL,RH,PA=0)
{
  var PA=0
  //---Initial Calculation
  var FNPST = Math.exp(16.6536-4030.183/(TA+235))    //saturated vapor pressure KPa
  if(PA==0) PA=RH*10*FNPST                  //water vapor pressure, Pa
  var ICL=0.155*CLO                              //thermal insulation of the clothing in m2K/W
  var M=MET*58.15                                //metabolic rate in W/m2
  var W=WME*58.15                                //external work in W/m2
  var MW=(M-W)                                   //internal heat production in the human body
  
  var FCL=null
  if(ICL<=0.078){var FCL=1+1.29*ICL} else {var FCL=1.05+0.645*ICL}
    
  //---heat transf. coeff. by forced convection
  var hcf=12.1*Math.sqrt(VEL)
  var taa=TA+273
  var tra=TR+273
  var tcla=taa+(35.5-TA)/(3.5*ICL+0.1)
  var p1=ICL*FCL
  var p2=p1*3.96
  var p3=p1*100
  var p4=p1*taa
  var p5=(308.7-0.028*MW)+(p2*(tra/100)**4)
  var xn=tcla/100
  var xf=tcla/50
  var eps=0.00015
  
  var n=0
  while (Math.abs(xn-xf)>eps)
  {
    var xf=(xf+xn)/2
    var hcn=2.38*(Math.abs(100*xf-taa)**0.25)
    if (hcf>hcn) {var hc=hcf} else {var hc=hcn}
    
    var xn=(p5+p4*hc-p2*(xf**4))/(100+p3*hc)
    
    var n=n+1
    
    if (n>150)
    {
      //print("Maximum number of iterations exceeded in PMV low air speed model n=150...")
      break;
    }
  }
      
  var tcl=100*xn-273
  
  //heat loss diff. through skin 
  var hl1=3.05*0.001*(5733-(6.99*MW)-PA)
  //heat loss by sweating
  if (MW>58.15) {var hl2=0.42*(MW-58.15)} else {var hl2=0}
  //latent respiration heat loss 
  var hl3=1.7*0.00001*M*(5867-PA)
  //dry respiration heat loss
  var hl4=0.0014*M*(34-TA)
  //heat loss by radiation  
  var hl5=3.96*FCL*(xn**4-(tra/100)**4)
  //heat loss by convection
  var hl6=FCL*hc*(tcl-TA)
  
  var ts=0.303*Math.exp(-0.036*M)+0.028
  
  
  //---Predicted Mean Vote - PMV
  PMV=ts*(MW-hl1-hl2-hl3-hl4-hl5-hl6)
  
  //------Predicted Percentage Dissatisfied - PPD
  var PPD=100-95*Math.exp(-0.03353*PMV**4-0.2179*PMV**2)
  var result=[parseFloat(PMV.toFixed(2)),parseFloat(PPD.toFixed(0))]
  
  return(result)
}

// função PMV elevate
var PMVelevated=function(CLO,MET,WME,TA,TR,VEL,RH,PA)
{

  //---Secant method of iteration
  var a1=0
  var b1=100
  var count1=1
  var F1=SETashrae(CLO,MET,WME,TA-a1,TR-a1,0.10,RH)-SETashrae(CLO,MET,WME,TA,TR,VEL,RH)
  if(Math.abs(F1)<0.001) {var ce=a1}
  var F2=SETashrae(CLO,MET,WME,TA-b1,TR-b1,0.10,RH)-SETashrae(CLO,MET,WME,TA,TR,VEL,RH)
  if(Math.abs(F2)<0.001) {var ce=b1}

  while (count1<=101)
  {
    var slope=(F2-F1)/(b1-a1)
    var c1=b1-F2/slope
    
    if(c1>999999 || c1==Infinity || c1==-Infinity)
    {
      var ce=NaN
      break;
    }
    
    var F3=SETashrae(CLO,MET,WME,TA-c1,TR-c1,0.10,RH)-SETashrae(CLO,MET,WME,TA,TR,VEL,RH)
    
    if(Math.abs(F3)<0.001)
    {
      var ce=c1
      break;
    }
    var a1=b1
    var b1=c1
    var F1=F2
    var F2=F3
    
    var count1=count1+1
        
    if(count1==101)
    {
      var ce=NaN
      //print("Maximum number of iteration exceeded in secant method of PMV elevated air speed model...")
      break;
    }
  }
  
  //---Bisection method of iteration
  var a2=-85
  var b2=+94
  var DIFF=0.5
  var count2=1
  if(ce=='NaN')
  {    
    while (Math.abs(DIFF)>0.002)
    { 
      var c2=(a2+b2)/2
      var Eq1=SETashrae(CLO,MET,WME,TA-a2,TR-a2,0.10,RH)-SETashrae(CLO,MET,WME,TA,TR,VEL,RH)
      var Eq2=SETashrae(CLO,MET,WME,TA-b2,TR-b2,0.10,RH)-SETashrae(CLO,MET,WME,TA,TR,VEL,RH)
      var Eq3=SETashrae(CLO,MET,WME,TA-c2,TR-c2,0.10,RH)-SETashrae(CLO,MET,WME,TA,TR,VEL,RH)
      
      if (Eq1*Eq3<0) {var b2=c2}
      if (Eq1*Eq3>0) {var a2=c2}
      if (Math.abs(Eq3)<0.0001) {break}
      
      var DIFF=abs(b2-a2)
      var ce=c2
      var count2=count2+1
            
      if(count2==101)
      {
        //print("Maximum number of iteration exceeded in bisection method of PMV elevated air speed model...")
        break;
      }
      
    }
  } else { 
    var ce=c1
  }

  //---Predicted Mean Vote - PMV
  var PMV=PMVlow(CLO,MET,WME,TA-ce,TR-ce,0.10,RH,PA)[0]
  
  //---Predicted Percentage Dissatisfied - PPD
  var PPD=100-95*Math.exp(-0.03353*PMV**4-0.2179*PMV**2)

  var result=[parseFloat(PMV.toFixed(2)),parseFloat(PPD.toFixed(0))]
  
  return(result)
}

// Função SET
var SETashrae=function (CLO,MET,WME,TA,TR,VEL,RH)
{

  var PA=101.325   //---Default of 1 atm equal to 101.325 kPa
  
  //---vapor pressure function
  var FindSVP=function(Temperature)
  {
  var resultadoSVP=Math.exp(18.6686-4030.183/(Temperature+235))
  return (resultadoSVP)
  }

  //---Constants
  var BodyWeight=69.9
  var BodySurfaceArea=1.8258
  var MetFactor=58.2
  var StefanB=5.6697*10**(-8) //---constante de Stefan Boltzmann [W/m2K4]
  var TTCR=36.49             //---neutral core body temperature
  var TTSK=33.7              //---neutral skin body temperature
  var TTBM=36.49             //---neutral body temperature
  var SKBFN=6.3              //---skin blood flow neutral
  var PressureInAtmospheres=PA*0.0098692 //---presure from kPa to atm

  //---Initial values
  var VaporPressure=RH*FindSVP(TA)/100   
  var AirVelocity=Math.max(VEL,0.1) //---0.1 m/s at least
  var M=MET*MetFactor      
  var RM=MET*MetFactor
  var RCL=0.155*CLO            //---m2K/W
  var FACL=1+0.15*CLO      
  
  var LR=2.2/PressureInAtmospheres
  var alpha=0.1                //---actual skin to total body mass ratio
  var MSHIV=0
  var ESK=0.1*MET

  //---Simulation initial values
  var TSK=TTSK   //---Temperature Skin
  var TCR=TTCR   //---Temperature Core
  var SKBF=SKBFN //---Skin Blood Flow
  
  //---Physiological temperatura regulation controls
  var CDIL=120   //----liters/(m2.h.K)
  var CSTR=0.5   //----dimensionless
  var CSW=170    //----g/(m2.h)
  var SKBFL=90   //----liter/(m2h) max SKBF
  var REGSWL=500
  
  var CHCt=3*(PressureInAtmospheres)**0.53
  var CHCV=8.600001*(AirVelocity*PressureInAtmospheres)**0.53
  var CHC=Math.max(CHCt,CHCV)
  var CHR=4.7              
  var CTC=CHR+CHC          
  var RA=1/(FACL*CTC)   
  var TOP=(CHR*TR+CHC*TA)/CTC   
  var TCL=TOP+(TSK-TOP)/(CTC*(RA+RCL))
  
  if (CLO<=0) {
    var WCRIT=0.38*AirVelocity**(-0.29)     
    var ICL=1
  } else {
    var WCRIT=0.59*AirVelocity**(-0.08)
    var ICL=0.45
  }
  
  ////////////////////////////////////////////////////---Begin 60 minutes simulation---////////////////////////////////////////////////////
  var LTIME=60   //---Timestep equal to 60
  var TIM=1
  var TCLold=TCL
  var flag=true
  while (TIM<=LTIME)
  {
    //---Dry heat balance - Solve TCL and CHR
    var DIFF1=0.5
    if(flag==true)
    {
      var count1=0
      while (DIFF1>0.001)
      {
        var TCLold=TCL
        var CHR=4*0.72*StefanB*((TCL+TR)/2+273.15)**3
        var CTC=CHR+CHC
        var RA=1/(FACL*CTC)
        var TOP=(CHR*TR+CHC*TA)/CTC
        var DIFF1=Math.abs(TCLold-TCL)
        var TCL=(RA*TSK+RCL*TOP)/(RA+RCL)
        var count1=count1+1
        
        //---Iteration from bisection method
        if(DIFF1>99999) //---an alternative control from non-convergente (rarely used)
        {
           var DIFF5=0.5
           var o1=-300
           var o2=+350
           var count2=0
           while(Math.abs(DIFF5)>0.001)
           {
             var o3=(o1+o2)/2
          
             var TCL=o1
             var CHR=4*0.72*StefanB*((TCL+TR)/2+273.15)**3
             var CTC=CHR+CHC
             var RA=1/(FACL*CTC)
             var TOP=(CHR*TR+CHC*TA)/CTC
             var funo1=(RA*TSK+RCL*TOP)/(RA+RCL)-TCL
          
             var TCL=o2
             var CHR=4*0.72*StefanB*((TCL+TR)/2+273.15)**3
             var CTC=CHR+CHC
             var RA=1/(FACL*CTC)
             var TOP=(CHR*TR+CHC*TA)/CTC
             var funo2=(RA*TSK+RCL*TOP)/(RA+RCL)-TCL
          
             var TCL=o3
             var CHR=4*0.72*StefanB*((TCL+TR)/2+273.15)**3
             var CTC=CHR+CHC
             var RA=1/(FACL*CTC)
             var TOP=(CHR*TR+CHC*TA)/CTC
             var funo3=(RA*TSK+RCL*TOP)/(RA+RCL)-TCL
          
             if(funo1*funo3<0) {var o2=o3}
		   if(funo1*funo3>0) {var o1=o3}
             if(Math.abs(funo3)<0.0001) break;
          
             var DIFF5=Math.abs(o2-o1)
             var TCL=o3
             var count2=count2+1
             
             if(count2==101)
             {
              //print("Maximum number of iteration exceeded in TCL inside SET from bisection method...")
              break;
             } 
           }
         }//---End of control of non-cnovergence
      }
    }
    var flag=false
    //---Heat flow from clothing surface to environment (FACL=1 if CLOE used)
    var DRY=(TSK-TOP)/(RA+RCL)                
    var HFCS=(TCR-TSK)*(5.28+1.163*SKBF)      

    //---Dry and latent respiratory heat losses
    var ERES=0.0023*M*(44-VaporPressure)     
    var CRES=0.0014*M*(34-TA)                

    var SCR=M-HFCS-ERES-CRES-WME
    var SSK=HFCS-DRY-ESK
  
    //---Thermal capacities
    var TCSK=0.97*alpha*BodyWeight      
    var TCCR=0.97*(1-alpha)*BodyWeight  
      
    //---Temperature changes in 1 minute
    var DTSK=(SSK*BodySurfaceArea)/(TCSK*60)   
    var DTCR=(SCR*BodySurfaceArea)/(TCCR*60)   
    var TSK=TSK+DTSK                //
    var TCR=TCR+DTCR                
    var TB=alpha*TSK+(1-alpha)*TCR  

    //---Definition of vascular control signals
    if (TSK>TTSK)
    {
      var WARMS=TSK-TTSK
      var COLDS=0
    } else {
      var COLDS=TTSK-TSK
      var WARMS=0
    }
    if (TCR>TTCR)
    {
      var WARMC=TCR-TTCR
      var COLDC=0
    } else {
      var COLDC=TTCR-TCR
      var WARMC=0
    }
    if (TB>TTBM)
    {
      var WARMB=TB-TTBM
      var COLDB=0
    } else {
      var COLDB=TTBM-TB
      var WARMB=0
    }

    //---Control skin blood flow
    var DILAT=CDIL*WARMC              
    var STRIC=CSTR*COLDS             
    var SKBF=(SKBFN+DILAT)/(1+STRIC) 
    
    //---SKBF is never below 0.5 liter/m2h not above SKBFL
    var SKBF=Math.max(0.5,Math.min(SKBFL,SKBF))
    
    //---control of regulatory sweting
    var REGSW=CSW*WARMB*Math.exp(WARMS/10.7)
    var REGSW=Math.min(REGSW,REGSWL)
    var ERSW=0.68*REGSW
      
    //---Mass transfer equation between skin and environment
    var REA=1/(LR*FACL*CHC)
    var RECL=RCL/(LR*ICL)
    var EMAX=(FindSVP(TSK)-VaporPressure)/(REA+RECL)
    var PRSW=ERSW/EMAX

    //---PDIF for nonsweating skin
    var PWET=0.06+0.94*PRSW   
    var EDIF=PWET*EMAX-ERSW   
    var ESK=ERSW+EDIF         

    //---Beginning of dripping (swet not evaporated on skin surface)
    if (PWET>WCRIT)
    {
      var PWET=WCRIT
      var PRSW=WCRIT/0.94          
      var ERSW=PRSW*EMAX
      var EDIF=0.06*(1-PRSW)*EMAX
      var ESK=ERSW+EDIF
    }

    //---When EMAX<0 condensation on skin ocurs
    if (EMAX<0)
    {
      var EDIF=0
      var ERSW=0
      var PWET=WCRIT
      var PRSW=WCRIT
      var ESK=EMAX
    }
  
    //---Adjustment of metabolic heat due to shivering
    var ESK=ERSW+EDIF
    var MSHIV=19.4*COLDS*COLDC
    var M=RM+MSHIV
    
    //---Ratio of skin-core masses change with SKBF
    var alpha=0.0417737+0.7451833/(SKBF+0.585417)
    
    var TIM=TIM+1 //---Next minute
  }
  ////////////////////////////////////////////////////---End 60 minutes simulation---////////////////////////////////////////////////////
  
  var HSK=DRY+ESK
  var RN=M-WME
  var ECOMF=0.42*(RN-(1*MetFactor))
  if (ECOMF<0) {var ECOMF=0}
  var EMAX=EMAX*WCRIT
  var W=PWET
  var PSSK=FindSVP(TSK)
  
  //---Standard environment
  var CHRS=CHR
  
  //---CHCS=standard convective heat transfer coefficient (level walking/still air)
  if (MET<0.85) {
    var CHCS=3
  } else {
    var CHCSt=5.66*(MET-0.85)**0.39
    var CHCS=Math.max(3,CHCSt)
  }

  var KCLO=0.25       
  var CTCS=CHCS+CHRS  
  var RCLOS=1.52/((MET-WME/MetFactor)+0.6944)-0.1835  
  var RCLS=0.155*RCLOS
  var FACLS=1+KCLO*RCLOS
  var FCLS=1/(1+0.155*FACLS*CTCS*RCLOS)
  var IMS=0.45
  var ICLS=IMS*CHCS/CTCS*(1-FCLS)/(CHCS/CTCS-FCLS*IMS)
  var RAS=1/(FACLS*CTCS)
  var REAS=1/(LR*FACLS*CHCS)
  var RECLS=RCLS/(LR*ICLS)
  var HDS=1/(RAS+RCLS)
  var HES=1/(REAS+RECLS)
  
  ////////////////////////////////////////////////////---Standard Effective Temperature SET---////////////////////////////////////////////////////
  
  //---Newton method of iteration
  var a3=TSK-HSK/HDS //---Initial estimation of SET
  var DIFF3=100
  var count3=0
  var delta=0.0001
  var SETOLD=a3
  while (Math.abs(DIFF3)>0.01)
  { 
    var ER1=HSK-HDS*(TSK-SETOLD)-W*HES*(PSSK-0.5*FindSVP(SETOLD))
    var ER2=HSK-HDS*(TSK-SETOLD-delta)-W*HES*(PSSK-0.5*FindSVP(SETOLD+delta))
    var SET=SETOLD-delta*ER1/(ER2-ER1)
    var DIFF3=SET-SETOLD
    var SETOLD=SET
    var count3=count3+1
    
    if(count3==101)
    {
      //print("Maximum number of iteration exceeded in SET iteration for Newton method")
      break;
    }
    
    if(Math.abs(DIFF3)>999999)
    {
      var SET=NULL
      var SETOLD=NULL
      var DIFF2=0.5
      var count4=0
      var a2=-301
      var b2=486
      while (DIFF2>0.001)
      {
        var c2=(b2+a2)/2
        var SETOLD=a2
        var ERR1=HSK-HDS*(TSK-SETOLD)-W*HES*(PSSK-0.5*FindSVP(SETOLD))
        var SETOLD=b2
        var ERR2=HSK-HDS*(TSK-SETOLD)-W*HES*(PSSK-0.5*FindSVP(SETOLD))
        var SETOLD=c2
        var ERR3=HSK-HDS*(TSK-SETOLD)-W*HES*(PSSK-0.5*FindSVP(SETOLD))
        if (ERR1*ERR3<0) {var b2=c2}
	  if (ERR1*ERR3>0) {var a2=c2}
        if (Math.abs(ERR3)<0.0001) break;
        var DIFF2=Math.abs(b2-a2)
        var SETOLD=c2
        var count4=count4+1
        
        if(count4==101)
        {
          //print("Maximum number of iteration exceeded in SET iteration for bisection method")
          break;
        }
        
        var SET=SETOLD //---iterative calculation
      }
    }
  }
    
  var result=parseFloat(SET.toFixed(2))
  return(result)
}

// função PMV both
var PMVboth=function(CLO,MET,WME,TA,TR,VEL,RH,PA)
{

  if (VEL<=0.10)
  {
    var PA=0
    PMV=parseFloat((PMVlow(CLO,MET,WME,TA,TR,VEL,RH,PA)[0]).toFixed(2))
  }
  if (VEL>0.10)
  {
    var PA=101.325
    PMV=parseFloat((PMVelevated(CLO,MET,WME,TA,TR,VEL,RH,PA)[0]).toFixed(2))
  }
    
  var resultado=PMV
  
  //---Predicted Percentage Dissatisfied - PPD
  var PPD=100-95*Math.exp(-0.03353*PMV**4-0.2179*PMV**2)

  var result=[parseFloat(PMV.toFixed(2)),parseFloat(PPD.toFixed(2))]
   
  return(result)
}