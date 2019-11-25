import numpy as np
from numpy import pi
from pylab import plt
from pylab import imshow
## Benining

## Setting parameters

##%%%% Fundamental constants %%%%
 # reduced Planck's constant
hb = 6.62607015/(2*pi) #10^(-34) J s
# Electron charge
ec = 1.602176634;  # 10^(-19) C
#Flux quantum
QC0 = 2*ec  # 10^(-19) C
QF0 = hb/QC0  # \mu V /GHz
QZ0 = 1e-4*QF0/QC0 # \Om = V/A

#  %%%%% Units %%%%%
# %Capacitance (C):       fF = 10^(-15) s^2 (A^2/J)
# %Inductance (L):        nH = 10^(-9) (J/A^2)  
# %Energy (E):            GHz = 10^(9) Hz
# %Coupling (g0):         (nH)^(-1) = G (A^2/J)
# %Norm coupling (g):     GHz = (10^(3)/2pi) * [g0 / C]
# %Frequency (w):         GHz = (10^(3)/2pi) * [LC]^(-1/2)
# %Flux (\Phi) :          \mu V / GHz = 10^(-2) [2pi \hbar v L]^(1/2)
# %Charge (Q) :           (1e-19) C = 10^(-1) [2pi \hbar v C]^(1/2)
# %Impedance (Z):         V/A = 10^(3) [L/C]^(1/2)


#  %%%% Notation %%%%
# %resonator
# %fluxonium (qubit)

# %%% Full circuit %%%
# %Capacitances
Cr =  20.3 
Cq = 5.3
# %Inductances
Lr0 = 15.6;
Lq0 = 386.0;
Lc0 = 4.5;
# Lr0 = 2
# Lq0 = 500.0
# Lc0 = 1.0

# Josephson energy
E_J = 6.20

# %% %%% Simplified circuit %%%
# %Effective values
Lr_s = ( Lr0*Lc0 +Lc0*Lq0 +Lq0*Lr0 )**2 / \
    ((Lr0*(Lq0+Lc0)**2 +Lc0*Lq0**2) +Lq0*Lc0**2)
Lq_s = ( Lr0*Lc0 +Lc0*Lq0 +Lq0*Lr0 )**2 / \
    (Lq0*(Lr0+Lc0)**2 +Lc0*Lr0**2 +Lr0*Lc0**2)
    
    
# Coupling
g0 = ( Lq0*Lc0*(Lr0 +Lc0) +Lr0*Lc0**2 ) \
    / ( Lr0*Lc0 +Lc0*Lq0 +Lq0*Lr0 )**2
# Normalized coupling (GHz)
g = (10**3/(2*pi))*np.sqrt(g0/(np.sqrt(Cr*Cq)))
#Turns ratio
r_cq = Lc0/Lq0

# %%Energy units coupling (GHz) 
#g_n = g0*((10**4)*(QF0**2)/hb)

#%Ideal values: Lq0 >> Lr0 ~ Lc0
Lr_ideal = (Lr0+Lc0)
g_ideal = r_cq/Lr_ideal

#%Data matrix
#%data_L = [[Lr0,Lr_ideal,Lr,0]';[Lq0,Lq,0]';Lc0]

#%Individial systems frequencies
vr0i = (10**3/(2*pi))/(np.sqrt(Lr0*Cr))
vq0i = (10**3/(2*pi))/(np.sqrt(Lq0*Cq))
## %Uncoupled (joint-system) frecuencies (GHz)
vr0 = (10**3/(2*pi))/(np.sqrt(Lr_s*Cr))
vq0 = (10**3/(2*pi))/(np.sqrt(Lq_s*Cq))

#%% Dressed modes %%
#%Coupled frequencies/(np.sqrt(Lr0*Cr))
vr = np.sqrt( (1/2)*(vr0**2 +vq0**2) +np.sqrt( (1/4)*(vq0**2 - vr0**2)**2 +g**4) )
vq = np.sqrt( (1/2)*(vr0**2 +vq0**2) -np.sqrt( (1/4)*(vq0**2 - vr0**2)**2 +g**4) )

# %Dressed modes inductances
Lr = (10**6/(2*pi)**2)/(Cr*vr**2)
Lq = (10**6/(2*pi)**2)/(Cq*vq**2)

# %Mixing mix = (theta/2)
mix = np.arctan(2*g**2/abs(vr0**2-vq0**2))/2;
# %Unitary transform to dressed modes
Kd = np.array([[np.cos(mix), np.sqrt(Cq/Cr)*np.sin(mix)], [-np.sqrt(Cr/Cq)*np.sin(mix), np.cos(mix)]])

# %Lambda coefficients
lr = Kd[1,0]
lq = Kd[1,1]
1-lq

# %Zero point fluctuations (ZPF)
# %QF0 %Quantum of flux
# %QC0 %Cooper pair charge
# %QZ0 %Quantum of impedance

# %Resonator
# %Flux (\mu V /GHz)
Fr_ZPZ = np.sqrt(hb*2*pi*vr*Lr)/10**2;
# %Charge (10**(-9) C)
Qr_ZPZ = np.sqrt(hb*2*pi*vr*Cr)/10;

# %Normalized flux (QF0)
Fr0 = Fr_ZPZ/QF0
# %Normalized charge (QC0)
Qr0 = Qr_ZPZ/QC0
Fr0Qr0 = Fr0*Qr0

# %Impendance (\Oms)
Zr = (10**3)*np.sqrt(Lr/Cr);
# %Normalized impendance (\Oms)
Zr0 = Zr/QZ0

# %Fluxonium
# %Flux (\mu V /GHz)
Fq_ZPZ = np.sqrt(hb*2*pi*vq*Lq)/10**2;
# %Charge (10**(-9) C)
Qq_ZPZ = np.sqrt(hb*2*pi*vq*Cq)/10;
# %Normalized flux (QF0)
Fq0 = Fq_ZPZ/QF0
# %Normalized charge (QC0)
Qq0 = Qq_ZPZ/QC0
Fq0Qq0 = Fq0*Qq0

# %Impendance (\Oms)
Zq = (10**3)*np.sqrt(Lq/Cq);
# %Normalized impendance (\Oms)
Zq0 = Zq/QZ0

# %% %%% Initial steps for numerical diagonalization %%%
# %z coeficients
z_r = Fr0*abs(lr)
z_q = Fq0*lq

#%%% Hilbert space dimension %%%
#%Antennae
Mr_max = 10
#%Fluxonium
Mq_max = 20

#%%% Normalized energies (vq) %%%
eq = 1
er = vr/vq
eJ = E_J/vq

#%%% Normalized energies (vr) %%%
#% eq = vq/vr
#% er = 1 %vr/vq
#% eJ = E_J/vr

#%%% Energy matrices %%%
#%Resonator
Hr_mat = er*np.diag(range(Mr_max))


#%%
#%Fluxonium
Hq_mat = eq*np.diag(range(Mq_max))

# Dressed hamiltonian %%%
H_dr = np.kron(Hr_mat, np.eye(Mq_max)) + np.kron(np.eye(Mr_max), Hq_mat)

#% figure(3); clf(3)
#% imagesc(H_dr)
#% colorbar

# Coupling matrices: antennae %%%
#clear z
#clear M_max

z = z_r
M_max = Mr_max
#%%
import numpy as np
from math import factorial
#%Interaction matrices
#%Numerical diagonalization of fluxonium and antennae
#%
#%This function calculates the full interaction matrix as well as
#%the symmetric one ( Re[I] = cos() ) and
#%the anti-symmetric one ( Im[I] = sin() ).
#%
#%It depends on the energy offset "m",
#%the relative energy difference "a" between the levels,
#%and the strenth of the Flux ZPF "z_x" 
#
#%Strength of the flux ZPF
#%zr = 0.0540;
#%zq = 2.8903;

#%z = zq;
#
#%Hilbert space dimension
#%M_max = 20;
#



def interaction_function_1(a, z):
    """
    Interaction function 1
    Numerical diagonalization of fluxonium and antennae
    
    This function is independent of the "base" energy quantum
    (energy offset), it only takes into account the relative 
    energy difference between the levels.
    
    Order parameter "a" (difference of energy quanta)
    z coefficient "z" (relative strength of the "zero point fluctiations")
    """
    return ( 1/( np.sqrt( 1.*factorial(a)*(2**a) ) ) )*(z**a)*np.exp(-(1/4)*z**2)

def interaction_function_2(m, a, z):
    """
    Numerical diagonalization of fluxonium and antennae
    This function is depends on both the "base" energy quanta
    (energy offset) and the the relative energy difference 
    between the levels.

    Offset energy quanta "m"
    Order parameter "a" (difference of energy quanta)
    z coefficient "z" (relative strength of the "zero point fluctiations")
    """
    if m==0:
        F_maz = 1
    else:
        Gmma = np.sqrt(1.*factorial(m+a) / ( factorial(a) *factorial(m) ) )
        Pmma = 0
        for vv in range(m+1):
            Pmma = Pmma + ((-1)**vv)* ( factorial(a) *factorial(m) / \
                ( (2**vv) * factorial(vv)*factorial(vv+a)*factorial(m-vv) ) )\
                *( z**(2*vv) )
        F_maz = Gmma*Pmma
    return F_maz

def interaction_matrices(M_max, z):
    Csg_a = 0
    #% Interaction matrix
    M_int = np.zeros((M_max, M_max))
    #%Symmetric interaction matrix "cos()"
    S_int = np.zeros((M_max, M_max))
    #%Anti-symmetric interaction matrix "cos()"
    A_int = np.zeros((M_max, M_max))
    
    #%Interaction function 1
    C1_az = np.zeros(M_max)
    for a in range(M_max):
        F_az = interaction_function_1(a, z)
        if np.mod(a, 2) == 0:
            Csg_a = (-1)**(a/2)
        else:
            Csg_a = (-1)**((a-1)/2)
        C1_az[a] = F_az*Csg_a
    
    for m in range(M_max):
        for a in range(M_max-m):
            print(m, a, z)
            F_maz = interaction_function_2(m, a, z)
            print(F_maz)
            
            #Full interaction matrix
            M_int[m,m+a] = C1_az[a]*F_maz
            M_int[m+a,m] = C1_az[a]*F_maz
            
            if np.mod(a,2) == 0:
                #%Symmetric
                S_int[m,m+a] = C1_az[a]*F_maz
                S_int[m+a,m] = C1_az[a]*F_maz
                #%Anti-symmetric
                A_int[m,m+a] = 0
                A_int[m+a,m] = 0
            else:
                #%Symmetric
                S_int[m,m+a] = 0
                S_int[m+a,m] = 0
                #%Anti-symmetric
                A_int[m,m+a] = C1_az[a]*F_maz
                A_int[m+a,m] = C1_az[a]*F_maz
    return S_int, A_int


Sr_int, Ar_int = interaction_matrices(M_max, z)


#%%% Coupling matrices: fluxonium %%%
z = z_q
M_max = Mq_max

Sq_int, Aq_int = interaction_matrices(M_max, z)

if True:
    plt.figure()
    plt.imshow(Sq_int + Aq_int)

# %%% Interaction matrices: Antennae x Fluxonium %%%
V_SS = np.kron(Sr_int,Sq_int)
V_AA = np.kron(Ar_int,Aq_int)
V_SA = -np.kron(Sr_int,Aq_int)# %Modification * (-1)
V_AS = np.kron(Ar_int,Sq_int)# %Modification * (-1)

#% V_SS = kron(eye(Mr_max),Sq_int);
#% V_AA = kron(zeros(Mr_max),Aq_int);
#% V_SA = -kron(eye(Mr_max),Aq_int); %Modification * (-1)
#% V_AS = kron(zeros(Mr_max),Sq_int); %Modification * (-1)

VS = (V_SS+V_AA)
VA = (V_SA+V_AS)

if True:
    plt.figure()
    imshow(V_SS)
    plt.colorbar()
    
    plt.figure()
    imshow(V_AA)
    plt.colorbar()
    plt.figure()
    imshow(V_SA)
    plt.colorbar()
    plt.figure()
    imshow(V_AS)
    plt.colorbar()
    plt.figure()
    imshow(VS)
    plt.colorbar()
    plt.figure()
    imshow(VA)
    plt.colorbar()


# %% %%% Numerical diagonalization %%%
# %Normalized external flux (QF0)
F_max = 50

NE_max = 10
E_rq = np.zeros((NE_max,F_max))

FF_ext = np.linspace(0, pi, 50)

for jj, F_ext in enumerate(FF_ext):
    # Full interaction matrix %%%
    V_rq = np.cos(F_ext)*VS + np.sin(F_ext)*VA
    V_rq = eJ*V_rq

    # Full Hamiltonian %%%
    H_full = H_dr - V_rq

    if False: pass
        #figure(10); clf(10)
        #%surf(H_full)
        #PlotDensityMatrix(U_ND(:,3)*U_ND(:,3)',1)
        #colorbar
     
    # Eigenvalues and eigenvectors
    E_ND, U_ND = np.linalg.eig(H_full)
    E_ND = sorted(E_ND)
    for ll in range(NE_max):
        E_rq[ll, jj] = E_ND[ll]-E_ND[0]
    
# %%
#%%% Figure 2 for paper // Heat
#fprintf('Axis fontsize \n')
#ffs = 20 %Axis labels
#
#fprintf('Ticks fontsize \n')
#tsx = 15 %Ticks fontsize
#
#fprintf('Markersize \n')
#mks = 4 % Markersize
#
#fprintf('Axis linewidht \n')
#axlwd = 1.2 %Axis linewidht
#
#fprintf('Cartesian-axis width \n')
#crwd = 0.5 %Cartesian-axis width
#
#fprintf('Cartesian-axis line style \n')
#crls = '-' %Cartesian-axis line style
#
#fprintf('Errorbar-cap size \n')
#ebcp = 0.1 %Errorbar-cap size
#%nmx = 14; %Subplot label size
#
#fprintf('Legend fontsize \n')
#lfz=10 %Legend fontsize
#
#fprintf('Results-line width \n')
#%aaf = 12; %Axis labels
#rlwd = 1.1 %Results-line width
#
#fprintf('Errorbars width \n')
#ebwd = 1.4 %Errorbars width
#
#fprintf('Demon readout: \n')
#dmnx = 0% demon ON: 1; OFF:0
#
#
#%Axis
#AxVec = [0.0 0.50 0 20];
#XTh = [0.0:0.05:0.5];
#XThL = {'0',' ','0.1',' ','0.2',' ','0.3',' ','0.4',' ',' '};
#YTh = [0:2:18];
#YThL = {'0',' ','4',' ','8',' ','12',' ','16'};
#
#%xx=2;
#
#figure(12); clf(12);
#AXS0 = gca;
#%surf(E_rq)
plt.figure()
for ll in range(1,10):
    plt.plot(FF_ext/(2*pi),vq*E_rq[ll,:],'-')


#
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%% Fancy axis properties %%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#
#%%%%%%%% Left and down %%%%%%%%%
#%AXS0 = gca;
#%%%%Fig2 position
#set(AXS0, 'Units', 'normalized');
#set(AXS0, 'Position', [0.150, 0.175, 0.825, 0.80]);
#%%%%Range
#axis(AxVec)
#%axis([-0.05 1.25 0.475 0.810])
#%%%%Tick marks and labels
#xticks(XTh)
#xticklabels(XThL)
#%set(AXS0,'XTickLabel',[]);
#yticks(YTh)
#yticklabels(YThL)
#%%%%Tick size
#set(AXS0,'FontSize',tsx);
#%%%%Labels
#xlabel('$\varphi_\textrm{ext}$ ','Interpreter','latex','FontSize',ffs)
#ylabel('$\nu$ (GHz)','Interpreter','latex','FontSize',ffs)
#%%%%Color
#%set(AXS0,'XColor',k)
#%set(AXS0,'YColor',k)
#%%%%Axis width
#set(AXS0,'linewidth',axlwd)
#
#%%%%%%%% Up and right %%%%%%%%
#ax1_pos = get(AXS0,'Position'); % position of first axes
#ax2 = axes('Position',ax1_pos,...
#    'XAxisLocation','top',...
#    'YAxisLocation','right',...
#    'Color','none');
#%%%%Range
#axis(AxVec)
#%axis([-0.05 2.00 0.475 0.810])
#%%%%Tick marks and labels
#xticks(XTh)
#%xticklabels({'0',' ','0.1',' ','0.2',' ','0.3',' ','0.4',' ','0.5'})
#set(ax2,'XTickLabel',[]);
#yticks(YTh)
#%yticklabels({'0.0',' ','5.0',' ','10.0',' ','15.0',' ','20.0'}
#set(ax2,'YTickLabel',[]);
#%%%%Tick size
#set(ax2,'FontSize',tsx);
#%%%%Labels
#%xlabel('$t_{\textrm{inj}}$ (ms)','Interpreter','latex','FontSize',ffs)
#%ylabel('QND transfer','Interpreter','latex','FontSize',ffs)
#%%%%Color
#%set(ax2,'XColor',k)
#%set(ax2,'YColor',k)
#%%%%Axis width
#set(ax2,'linewidth',axlwd)
#
#%title('Heat vs inverse temperature','Interpreter','latex','FontSize',14)
#
#% Saving Figure 2
#% if dmnx ==1
#%     verx = 'ON'; %Version
#% else
#%     verx = 'OFF'; %Version
#% end
#% verx = 'v3';
#Figx = 'Fluxonium coupled to antennae'; %Figure
#set(gcf, 'color','white');
#%print(gcf, '-r600', '-djpeg', [PicDir Figx filesep Figx '_' verx '_HD']);
#print(gcf, '-r600', '-dpng', ['.\' Figx '_HD']);
#savefig(['.\' Figx ])