import numpy as np
# Initial values related to the state variables
x0 = 0
y0 = np.array([0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]) #Rotor side control simulaton input initialization
y1 = np.array([0.0,0.0,1]) # FEC input initialization ifed,ifeq,Vdc
pastDiff=[0,0,0,0,0,0,0,0];
Pset_i=1;Qset_i=0;

t=1 #simulation end time in seconds
h = 1e-5 #time step

n=(t-x0)/h; #no of iterations

#Empty class of all variables that are to be iterated
class inputs:
    Ls=0;Lr=0;Lm=0;Vrd=0;Vrq=0;Vsd=0;Vsq=0;Vs=0;ird=0;irq=0;isd=0;isq=0;We=0;Wr=0;Rr=0;Rs=0;dFrd=0;
    dFrq=0;dFsd=0;dFsq=0;pole=0;f=0;simTime=0;nt=0;Nr=0;Frd=0;Frq=0;Fsd=0;Fsq=0

inputToODE = inputs #inputToODE will be provided to solver
inputToODE.Vrd=0;inputToODE.Vrq=0;inputToODE.Vsd=0;inputToODE.ird=y0[0];inputToODE.irq=y0[1];
inputToODE.isd=y0[2];inputToODE.isq=y0[3];inputToODE.dFrd=pastDiff[0];inputToODE.dFrq=pastDiff[1];
inputToODE.dFsd=pastDiff[2];inputToODE.dFsq=pastDiff[3];inputToODE.simTime=t;
R_q_PI_I=0;R_d_PI_Io=0;R_d_PI_Ii=0;
#Global variables
Vrd_g=[];Vrq_g=[];ird_g=[];irq_g=[];isd_g=[];isq_g=[];Vsd_g=[];Vsq_g=[];t_g=[];Vfed_g=[0.0001];Vfeq_g=[0.0001];ifed_g=[0.0001];ifeq_g=[0.0001];Vdc_g=[];Vsd_g=[];Vsq_g=[];
dydx1=np.array(np.zeros(8))
dydx=np.array(np.zeros(3))
vsa=[];isa=[];ira=[];ird=[];irq=[];isd=[];isq=[];Wr1=[];Te=[];P=[];Q=[];x1=[];
count_P=0;count_Q=0;count_V=0;
#Rotor side control set points
Pset = [1.0,0.8] #Active Power set point
Qset = [0.0,0.001] #Reactive Power set point

#no of changes in Active Power requested 
changes_no_P=len(Pset)+1
times_P=np.around(np.linspace(0,int(t/h),changes_no_P),decimals=0)

#no of changes in Active Power requested 
changes_no_Q=len(Qset)+1
times_Q=np.around(np.linspace(0,int(t/h),changes_no_Q),decimals=0)

a=[1.0]
# for i in range(1000):
    # a.append(1.0)
# k=0.15
# for i in range(200):
    # a.append(k)
# for i in range (850):
    # a.append(k)
    # k+=0.001
# for i in range(1000):
    # a.append(1.0)

#no of changes in Supply Voltage requested
#Vset=[1.0,1.0,1.0,1.0,1.0,1.0,0.05,0.05,0.05,1.0,1.0,1.0,1.0,0.8,0.9,1.0,0.8,0.85,0.06,0.1,1.0];
Vset = [1,1.001] #a
changes_no_V=len(Vset)+1
times_V=np.around(np.linspace(0,int(t/h),changes_no_V),decimals=0)



#FEC set point
Vdc_set=415 #Volts

Qfed_set=900 # VAr
Ifed_set= Qfed_set/(1*230)/1.732 #pu
F_q_PI_Io=0;F_q_PI_Ii=0;F_d_PI_Ii=0;F_d_PI_Io=0;


#Machine parameters initialization
inputToODE.Rr=0.8;inputToODE.Rs=1.06;inputToODE.Ls=206.5e-3;inputToODE.Lr=81e-3;inputToODE.Lm=66.4e-3;inputToODE.Vs=230;
inputToODE.f=50;inputToODE.Nr=400;inputToODE.pole=6;inputToODE.nt=0;inputToODE.Vrd=0;inputToODE.Vrq=0;
inputToODE.Wr=2*3.1416*(inputToODE.Nr*inputToODE.pole)/120;
inputToODE.We=2*3.1416*inputToODE.f;
inputToODE.Vsq=1.414*inputToODE.Vs;
inputToODE.Lfe=12e-3;
inputToODE.Rfe=0.1;
inputToODE.C=2.4e-3;
inputToODE.Vfeq=0;inputToODE.Vfed=0;
sigma=0;
Prated=7500
Qrated=7500

#initialization for voltage excursion simulations
Vsd_initial = inputToODE.Vsd
Vsq_initial = inputToODE.Vsq
#print(Vsd_initial,Vsq_initial)


inputToDiff = inputToODE
#DFIG model
def dfig(inputFromSolver, store=0):
#     global ird_o,irq_o,isd_o,isq_o;
    #Get all current state of variables during the solving
    Vsd=inputToDiff.Vsd;Vsq=inputToDiff.Vsq;Rr=inputToDiff.Rr;
    Rs=inputToDiff.Rs;We=inputToDiff.We;Wr=inputToDiff.Wr;ird=inputToDiff.ird;irq=inputToDiff.irq;isd=inputToDiff.isd;
    isq=inputToDiff.isq;dFrd=inputFromSolver[0];dFrq=inputFromSolver[1];dFsd=inputFromSolver[2];dFsq=inputFromSolver[3];
    Lm=inputToDiff.Lm;Ls=inputToDiff.Ls;Lr=inputToDiff.Lr;
    [Vrq,Vrd] = rotor_ctrl(inputToDiff,store)   
    
    dydx1[0] = (dFrd-Lm*isd)/Lr;         #ird equation
    dydx1[1] = (dFrq-Lm*isq)/Lr;         #irq equation
    dydx1[2] = (dFsd-Lm*ird)/Ls;         #isd equation
    dydx1[3] = (dFsq-Lm*irq)/Ls;         #isq equation
    
    dydx1[4]= (Vrd-Rr*ird+(We-Wr)*dFrq); #Flux Frd
    dydx1[5]= (Vrq-Rr*irq-(We-Wr)*dFrd); #Flux Frq
    dydx1[6]= (Vsd-Rs*isd+(We)*dFsq);    #Flux Fsd
    dydx1[7]= (Vsq-Rs*isq-(We)*dFsd);    #Flux Fsq 
    
    #store values for analysis
    if (store==1):     
        ird_g.append(dydx1[0]);
        irq_g.append(dydx1[1]);
        isd_g.append(dydx1[2]);
        isq_g.append(dydx1[3]);
        Vsd_g.append(inputToDiff.Vsd);
        Vsq_g.append(inputToDiff.Vsq);
        inputToDiff.ird=dydx1[0];inputToDiff.irq=dydx1[1];inputToDiff.isd=dydx1[2];inputToDiff.isq=dydx1[3];

#     ird_o=dydx1[0];irq_o=dydx1[1];isd_o=dydx1[2];isq_o=dydx1[3];   
    return dydx1[4:8] #Return calculated curents and flux


# In[51]:


#Rotor Control
def rotor_ctrl(inputToDiff,store=0):
    #Vrq equation
    global R_q_PI_I, Vrq_g,count,Pset_i,Qset_i;#Pset_i=1;Qset_i=0
    isq_set = Pset_i/(1.5*inputToDiff.Vsq)
    irq_set = isq_set * inputToDiff.Ls/inputToDiff.Lm
    Kp_rotor_q = 15
    Ki_rotor_q = 0.05
       
    #PI controller
    R_q_PI_P = (irq_set - inputToDiff.irq) * Kp_rotor_q           
    R_q_PI_I +=  (irq_set - inputToDiff.irq) * Ki_rotor_q
    R_q_PI = R_q_PI_P + R_q_PI_I # PI control output variable
    
    #Cross coupling terms
    Vrq_emf = 0
    Vrq_cross =0
    Vrq = R_q_PI - Vrq_emf - Vrq_cross

    
    #Vrd equation

    global R_d_PI_Io,Vrd_g;
    Kp_rotor_do =-0.02 # Outer loop P constant
    Ki_rotor_do =-0.0001 ## Outer loop I constant
    
    #PI controller - outer loop - d axis
    Qact = 1.5 * (inputToDiff.Vsq * inputToDiff.isd - inputToDiff.Vsd * inputToDiff.isq)
    R_d_PI_Po = (Qset_i - Qact) * Kp_rotor_do
    R_d_PI_Io += (Qset_i - Qact) * Ki_rotor_do
    R_d_PI_o = R_d_PI_Po + R_d_PI_Io
    
    ird_set = R_d_PI_o - inputToDiff.isd * inputToDiff.Ls/inputToDiff.Lm
    
    #PI controller - inner loop - d axis
    global R_d_PI_Ii
    Kp_rotor_di = 0.02 # Inner loop P constant
    Ki_rotor_di = 0.0000002 ## Inner loop I constant    
    R_d_PI_Pi = (ird_set - inputToDiff.ird) * Kp_rotor_di
    R_d_PI_Ii += (ird_set - inputToDiff.ird) * Ki_rotor_di
    R_d_PI_i = R_d_PI_Pi + R_d_PI_Ii
    
    Vrd = R_d_PI_i
    
    #store values for analysis    
    if (store==1):    
        Vrq_g.append(Vrq)
        Vrd_g.append(Vrd)
    
    return Vrq,Vrd


# In[52]:


#FEC model
def fec(inputFromSolver,store=0):
    global count,Vfed_g,Vfeq_g,ifed_g,ifeq_g,Vdc_g,Vfed_o,Vfeq_o;
    #Get all current state of variables during the solving
    Vsd=Vsd_g[len(Vsd_g)-1];Vsq=Vsq_g[len(Vsq_g)-1];
    Vrd=Vrd_g[len(Vrd_g)-1];Vrq=Vrq_g[len(Vrq_g)-1]; 
    ird=ird_g[len(ird_g)-1];irq=irq_g[len(irq_g)-1];
#     Vrd=0;Vrq=0;ird=0;irq=0;
    Rfe=inputToDiff.Rfe;Lfe=inputToDiff.Lfe;We=inputToODE.We;C=inputToODE.C;
    Vfed = Vfed_g[len(Vfed_g)-1]; Vfeq = Vfeq_g[len(Vfeq_g)-1];
    ifed = inputFromSolver[0];ifeq = inputFromSolver[1];Vdc = inputFromSolver[2];
    
    [Vfed, Vfeq]=fec_ctrl(Vfed,Vfeq,ifed,ifeq,Vdc,Vsq,Vsd)
    
    dydx[0] = (Rfe/Lfe * (-ifed+(Vsd-Vfed)/Rfe-Lfe/Rfe*We*ifeq)); #ifed equation
    dydx[1] = (Rfe/Lfe * (-ifeq+(Vsq-Vfeq)/Rfe+Lfe/Rfe*We*ifed)); #ifeq equation
    dydx[2] = (1/C * (Vsq/Vdc*ifeq-((Vrq*irq+Vrd*ird)/Vdc)));     #Vdc equation, clipped to max possible values
    
    #store values for analysis
    if (store==1):
        Vfed_g.append(Vfed_o);
        Vfeq_g.append(Vfeq_o);
        ifed_g.append(ifed);
        ifeq_g.append(ifeq);
        Vdc_g.append(Vdc);
    Vfed_o = Vfed
    Vfeq_o = Vfeq
    return np.array([dydx[0], dydx[1],dydx[2]]) #Return calculated FEC currents & DC link Voltage

def fec_ctrl(Vfed,Vfeq,ifed,ifeq,Vdc,Vsq,Vsd):
    
    global F_q_PI_Io,F_q_PI_Ii,F_d_PI_Ii,F_d_PI_Io;
    
    Lfe=inputToDiff.Lfe;We=inputToODE.We;
    #Vfeq calculation    
    #PI controller - outer loop - q axis     
    Kp_fec_qo = 500e-1  # Outer loop P constant
    Ki_fec_qo = 1e-6 # Outer loop I constant       
    F_q_PI_Po = (Vdc_set - Vdc) * Kp_fec_qo
    F_q_PI_Io += (Vdc_set - Vdc) * Ki_fec_qo
    F_q_PI_o = F_q_PI_Po + F_q_PI_Io #Ifeq set point
    
    #PI controller - inner loop - q axis       
    Kp_fec_qi = 500e-1 # inner loop P constant
    Ki_fec_qi = 1e-6    # inner loop I constant    
    F_q_PI_Pi = ((F_q_PI_o - ifeq) * Kp_fec_qi)
    F_q_PI_Ii += (F_q_PI_o - ifeq) * Ki_fec_qi
    F_q_PI = F_q_PI_Pi + (F_q_PI_Ii)        #Vfeq firing point ref
    Vfeq_c = - F_q_PI  + Vsq + We * Lfe * ifed      #Vfeq firing pt minus cross coupling
    
    #Vfed calculation   
    Kp_fec_d = 1 #P constant
    Ki_fec_d = 1e-2 #I constant    
    F_d_PI_Pi = ((Ifed_set - ifed) * Kp_fec_d)
    F_d_PI_Ii += (Ifed_set - ifed) * Ki_fec_d
    F_d_PI = F_d_PI_Pi + F_d_PI_Ii                            #Vfed firing point ref
    Vfed_c = -F_d_PI - Vsd - We * Lfe * ifeq              #Vfed firing pt minus cross coupling
    
    return [Vfed_c, Vfeq_c]


# In[53]:


def controls_grid(i):
    global count_P,count_Q,count_V,Pset_i,Qset_i;
    #calculate the no of actual change requests in Active Power
    if(i*h in times_P):
        Pset_i = Pset[count_P]*Prated;count_P=count_P+1
    #calculate the no of actual change requests in Reactive Power
    if(i*h in times_Q):
        Qset_i = Qset[count_Q]*Qrated;count_Q=count_Q+1
    #calculate the no of actual change requests in supply Voltage
    if(i in times_V):
        inputToDiff.Vsd=Vsd_initial*Vset[count_V];
        inputToDiff.Vsq=Vsq_initial*Vset[count_V];
        Pset_i = Pset[count_P-1]*Prated;
        Pset_i = np.clip(Pset_i,0,Prated*Vset[count_V]); #Limit active power proportional to voltage
        count_V+= 1
#         print(i,count_V-1,count_V)


# In[54]:


#Dormand Prince Runge Kutta solver 
def RKS(dydx,inputVector):
    y = inputVector
    k1 = h * dydx(inputVector) 
    k2 = h * dydx(inputVector + 1/5 * k1) 
    k3 = h * dydx(inputVector + 3/40 * k1 + 9/40 * k2) 
    k4 = h * dydx(inputVector + 44/45 * k1 -56/15 * k2 + 32/9 * k3) 
    k5 = h * dydx(inputVector + 19372/6561 * k1 -25360/2187 * k2 + 64448/6561 * k3 -212/729 * k4)
    k6 = h * dydx(inputVector + 9017/3168 * k1 -355/33 * k2 + 46732/5247 * k3 + 49/176 * k4 - 5103/18656 * k5) 
    y+= 35/384 * k1 + 500/1113 *k3 + 125/192 * k4 -2187/6784 * k5 + 11/84 * k6
    dydx(y,store=1)
    return y