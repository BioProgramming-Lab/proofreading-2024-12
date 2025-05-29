from scipy.integrate import odeint
import numpy as np

def dc_dt(c, t, derivs_0, derivs_L, diff_coeff_fun, diff_coeff_params, rxn_fun, rxn_params, n_species, h):
    """
    Time derivative of concentrations in an R-D system for constant flux BCs.
    
    Parameters
    c: ndarray, shape(n_species * n_gridpoints)
    The concentration of the chemical species interleaved in a Numpy array. The interleaving allows us to take
    advantage of the banded structure of the Jacobian when using the Hindmarsh algorithm for integrating in time.
    
    t: Float, time.
    
    derivs_0: ndarray, shape(n_species), derivs_0[i] is dc_i/dt at x = 0.
    
    derivs_L: ndarray, shape(n_species), derivs_L[i] is dc_i/dt at x = L.
    
    diff_coeff_fun: funtion of the form diff_coeff_fun(c_tuple, t, *diff_coeff_params), computing the diffusion term
    with a tuple, where entry i is a numpy array containing the ith species at the grid points.
    
    diff_coeff_params: arbitrary.
    rxn_fun: function of the form rxn_fun(c_tuple, t, *rxn_params). Return a tuple where entry i is a numpy array
    contatining the ith species at the grid points.
    
    rxn_params: arbitrary.
    
    n_species: int, number of chemical species (Bmp4, Chordin).
    
    h: float, Grid spacing (constant)
    
    Returns
    dc_dt: ndarray, shape(n_species * n_gridpoints) The time derivatives of the concentrations of the chemical species
    at grid points interleaved in a NumPy array.
    """
    # Tuple of concentrations
    c_tuple = tuple([c[i::n_species] for i in range(n_species)])
    
    # Compute diffusion terms
    D_tuple = diff_coeff_fun(c_tuple, t, *diff_coeff_params)
    
    # Compute reaction terms
    rxn_tuple = rxn_fun(c_tuple, t, *rxn_params)
    
    # Initiate return array
    conc_deriv = np.empty_like(c)
    
    # Empty array for storing concentrations
    da_dt = np.empty(len(c_tuple[0]))
    
    # Useful term
    h2 = h**2
    
    # Compute diffusion terms (central differencing w/ reflective BCs)
    for i in range(n_species):
        # View of concentrations of convenience
        a = np.copy(c_tuple[i])
        
        # Time derivative at left boundary
        # da_dt[0] = D_tuple[i][0] / h2 * 2 * (a[1] - a[0] - h * derivs_0[i])
        da_dt[0] = derivs_0[i] + D_tuple[i][0] / h2 * (a[1] - a[0])
        
        # Time derivative for middle grid points
        da_dt[1:-1] = D_tuple[i][1:-1] / h2 * np.diff(a, 2)
        
        # Time derivative for right boundary
        # da_dt[-1] = D_tuple[i][-1] / h2 * 2 * (a[-2] - a[-1] - h * derivs_L[i])
        da_dt[-1] = derivs_L[i] + D_tuple[i][-1] / h2 * (a[-2] - a[-1])
        
        # Store in output array with reaction terms
        conc_deriv[i::n_species] = da_dt + rxn_tuple[i]
    
    return conc_deriv

def RD_solve(c_0_tuple, t, L=1, derivs_0=0, derivs_L=0, diff_coeff_fun=None, diff_coeff_params=(), rxn_fun=None,
             rxn_params=(), rtol=1.49012e-6, atol=1.49012e-6):
    """
    Parameters
    ----------
    c_0_tuple : tuple
        c_0_tuple[i] is a NumPy array of length n_gridpoints with the 
        initial concentrations of chemical species i at the grid points.
    t : ndarray
        An array of time points for which the solution is desired.
    L : float
        Total length of the x-domain.
    derivs_0 : ndarray, shape (n_species)
        derivs_0[i] is the value of dc_i/dx at x = 0.
    derivs_L : ndarray, shape (n_species)
        derivs_L[i] is the value of dc_i/dx at x = L, the rightmost
        boundary of the domain of x.
    diff_coeff_fun : function
        Function of the form diff_coeff_fun(c_tuple, t, *diff_coeff_params).
        Returns an tuple where entry i is a NumPy array containing
        the diffusion coefficient of species i at the grid points.
        c_tuple[i] is a NumPy array containing the concentrations of
        species i at the grid poitns.
    diff_coeff_params : arbitrary
        Tuple of parameters to be passed into diff_coeff_fun.
    rxn_fun : function
        Function of the form rxn_fun(c_tuple, t, *rxn_params).
        Returns an tuple where entry i is a NumPy array containing
        the net rate of production of species i by chemical reaction
        at the grid points.  c_tuple[i] is a NumPy array containing 
        the concentrations of species i at the grid poitns.
    rxn_params : arbitrary
        Tuple of parameters to be passed into rxn_fun.
    rtol : float
        Relative tolerance for solver.  Default as odeint's default.
    atol : float
        Absolute tolerance for solver.  Default as odeint's default.
        
    Returns
    -------
    c_tuple : tuple
        c_tuple[i] is a NumPy array of shape (len(t), n_gridpoints)
        with the initial concentrations of chemical species i at 
        the grid points over time.
        
    Notes
    -----
    .. When intergrating for long times near a steady state, you
       may need to lower the absolute tolerance (atol) because the
       solution does not change much over time and it may be difficult
       for the solver to maintain tight tolerances.
    """
    # Number of grid points
    n_gridpoints = len(c_0_tuple[0])
    
    # Number of chemical species
    n_species = len(c_0_tuple)

    # Grid spacing
    h = (L / (n_gridpoints - 1.0))
    
    # Set up boundary conditions
    if np.isscalar(derivs_0):
        derivs_0 = np.array(n_species * [derivs_0])
    if np.isscalar(derivs_L):
        derivs_L = np.array(n_species * [derivs_L])

    # Set up parameters to be passed in to dc_dt
    params = (derivs_0, derivs_L, diff_coeff_fun, diff_coeff_params, 
              rxn_fun, rxn_params, n_species, h)

    # Set up initial condition
    c0 = np.empty(n_species * n_gridpoints)
    for i in range(n_species):
        c0[i::n_species] = c_0_tuple[i]

    # Solve using odeint, taking advantage of banded structure
    sol = odeint(
               dc_dt, c0, t, args=params, ml=n_species, mu=n_species,
               rtol=rtol, atol=atol)
    
    '''# Test whether the model can reach a steady state.
    precision = 1e-11
    sol1=sol[-1, :]
    sol2=sol[-2, :]
    time_0 = 60 #h
    j=0
    while np.mean(np.square(sol1/sol2-1))>precision : 
        j+=1
        t = np.linspace((j-1)*time_0,j*time_0,60)
        sol = odeint(
                    dc_dt, sol1, t, args=params, ml=n_species, mu=n_species,
                    rtol=rtol, atol=atol)
        sol1=sol[-1, :]
        sol2=sol[-2, :]
        if j>100 : 
            print("no steady state")'''    
    return tuple([sol[:,i::n_species] for i in range(n_species)])

class DIFFUSION(object):
    # Diffusion coefficient
    def __init__(self, **kwargs):
        """
        D_A: diffusion rate of A
        D_B: diffusion rate of B
        D_R: diffusion rate of R
        D_C: diffusion rate of A'
        D_complex: diffusion rate of complex
        """
        self.D_A = 10
        self.D_B = 10  
        self.D_C = 10 
        self.D_complex = 10

        
        # Put in params that were specified in input
        for entry in kwargs:
            setattr(self, entry, kwargs[entry])

def Diff_fun(c_tuple, t, diff_coeffs):
    # Compute D
    '''
    [A][B][C][R][AC][BC][AR][BR]
    '''
    D = np.zeros((len(c_tuple), len(c_tuple[0])))
    D[0] = diff_coeffs.D_A * np.ones(len(c_tuple[0])) # free A
    D[1] = diff_coeffs.D_B * np.ones(len(c_tuple[0])) # free B
    D[2] = diff_coeffs.D_C * np.ones(len(c_tuple[0])) # free A'
    D[3] = np.zeros(len(c_tuple[0])) # unbinded receptor
    D[4] = diff_coeffs.D_complex * np.ones(len(c_tuple[0])) # AA'
    D[5] = diff_coeffs.D_complex * np.ones(len(c_tuple[0])) # BA'
    D[6] = np.zeros(len(c_tuple[0])) # A_receptor complex
    D[7] = np.zeros(len(c_tuple[0])) # B_receptor complex

    return tuple(D)

def RD_rxn(c_tuple, t, rxn_coeffs, production_rate):
        """
        c_A: Free A
        c_B: Free B
        c_C: Free A'
        c_R: Receptor 
        c_AC: AA'
        c_BC: BA'
        c_AD: AD
        c_BD: BD
        c_AR: A-Receptor complex
        c_BR: B-Receptor complex
        """

        def hill(x):
            return (x/rxn_coeffs.k)**rxn_coeffs.n / (1 + (x/rxn_coeffs.k)**rxn_coeffs.n)
        
        def hill_ar(a, r):
            return (a/rxn_coeffs.k_ac)**rxn_coeffs.n_ac / \
                (
                    1 + 
                    (a/rxn_coeffs.k_ac)**rxn_coeffs.n_ac + 
                    (r/rxn_coeffs.k_rp)**rxn_coeffs.n_rp
                )

        # Unpack
        c_A, c_B, c_C, c_R, c_AC, c_BC, c_AR, c_BR= c_tuple
        
        #j_a, j_b: when the feedback respond to [AR] [BR]；j_ac,j_bc: when the feedback respond to [AC] [BC]
        j_A, j_B, j_C, j_a, j_b,j_c, j_ac, j_bc, j_R, j_ac_rp= production_rate
        
        # A reaction rate
        A_rate = j_A - rxn_coeffs.k_AC * c_A * c_C + rxn_coeffs.r_AC * c_AC \
                     - rxn_coeffs.k_AR * c_A * c_R + rxn_coeffs.r_AR * c_AR \
                     - rxn_coeffs.deg * c_A\
                     + j_a * hill(c_AR) + j_ac * hill(c_AC) + j_ac_rp * hill_ar(c_BC, c_AC)
                    
                    
        # B reaction rate
        B_rate = j_B - rxn_coeffs.k_BC * c_B * c_C + rxn_coeffs.r_BC * c_BC \
                     - rxn_coeffs.k_BR * c_B * c_R + rxn_coeffs.r_BR * c_BR \
                     - rxn_coeffs.deg * c_B \
                     + j_b * hill(c_BR) + j_bc * hill(c_BC)  + j_ac_rp * hill_ar(c_AC, c_BC)

        # A' reaction rate
        C_rate = j_C - rxn_coeffs.k_AC * c_A * c_C + rxn_coeffs.r_AC * c_AC \
                 - rxn_coeffs.k_BC * c_B * c_C + rxn_coeffs.r_BC * c_BC \
                 - rxn_coeffs.deg * c_C\
                 + j_c * hill(c_AR) + j_ac * hill(c_AC) + j_c * hill(c_BR) + j_bc * hill(c_BC)  + j_ac_rp * (hill_ar(c_BC, c_AC)+hill_ar(c_AC, c_BC))
         


        # Receptor reaction rate (assume constant turnover)
        '''
        R_rate = j_R - rxn_coeffs.k_AR * c_A * c_R + rxn_coeffs.r_AR * c_AR + rxn_coeffs.gamma1 * c_AR \
                 - rxn_coeffs.k_BR * c_B * c_R + rxn_coeffs.r_BR * c_BR + rxn_coeffs.gamma1 * c_BR
        '''   
        R_rate = j_R- rxn_coeffs.k_AR * c_A * c_R + rxn_coeffs.r_AR * c_AR + rxn_coeffs.gamma * c_AR \
                 - rxn_coeffs.k_BR * c_B * c_R + rxn_coeffs.r_BR * c_BR + rxn_coeffs.gamma * c_BR - rxn_coeffs.deg * c_R

        # AA' reaction rate
        ac_rate = rxn_coeffs.k_AC * c_A * c_C - rxn_coeffs.r_AC * c_AC - rxn_coeffs.deg * c_AC 
        
        # BA' reaction rate
        bc_rate = rxn_coeffs.k_BC * c_B * c_C - rxn_coeffs.r_BC * c_BC - rxn_coeffs.deg * c_BC 

        # A-Receptor reaction rate (assume constant turnover)
        ar_rate = rxn_coeffs.k_AR * c_A * c_R - (rxn_coeffs.r_AR + rxn_coeffs.gamma) * c_AR - rxn_coeffs.deg * c_AR            

        # B-Receptor reaction rate (assume constant turnover)
        br_rate = rxn_coeffs.k_BR * c_B * c_R - (rxn_coeffs.r_BR + rxn_coeffs.gamma) * c_BR - rxn_coeffs.deg * c_BR

        
        return (A_rate, B_rate, C_rate, R_rate, ac_rate, bc_rate, ar_rate, br_rate)

class RXN_params(object):
    """
    Container for reaction parameters
    """
    def __init__(self, **kwargs):
        """
        k_AC: association rate for A and Receptor
        r_AC: dissociation rate for AA'
        k_AR: association rate for A and Receptor
        r_AR: dissociation rate for A-Receptor complex
        k_BC: association rate for B and Receptor
        r_BC: dissociation rate for BA'
        k_BR: association rate for B and Receptor
        r_BR: dissociation rate for B-Receptor complex
        k_AD: association rate for B and D
        r_AD: dissociation rate for BD
        k_BD: association rate for B and D
        r_BD: dissociation rate for BD
        RTotal: Total receptor number
        gamma: receptor turnover rate ([X-R] --> R)
        deg: universal degradation rate for free proteins
        """
        # RTotal: Total receptor number (nM) from fitting
        #self.RTotal = 2.7 # nM
        # k (nM-1*s-1) and r (s-1) of A+A' <-> AA' 
        self.k_AC = 1e-4  # /nM/min
        self.r_AC = 1e-4  # /min
        # k (nM-1*s-1) and r (s-1) of B+A' <-> BA'
        self.k_BC = 1e-4  # /nM/min
        self.r_BC = 1e-3  # /min
        # k (nM-1*s-1) and r (s-1) of A+Receptor <-> A-Receptor complex
        self.k_AR = 4.5e-4  # /nM/min
        self.r_AR = 1e-3  # /min
        # k (nM-1*s-1) and r (s-1) of B+Receptor <-> B-Receptor complex
        self.k_BR = 4.5e-4  # /nM/min
        self.r_BR = 1e-3  # /min
        # Recycling rate (s-1) of A-Receptor and B-Receptor
        self.gamma = 4e-4  #/nM/min
        #/nM/min
        # Not considering basal degradation of species
        self.deg = 2e-5

        # Hill function
        self.k = 0
        self.n = 0
        
        self.ratio = self.r_BC/self.r_AC
        
        # Put in params that were specified in input
        for entry in kwargs:
            setattr(self, entry, kwargs[entry])