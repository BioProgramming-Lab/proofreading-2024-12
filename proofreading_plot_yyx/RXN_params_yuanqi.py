class RXN_params_yuanqi(object):
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
        RTotal: Total receptor number
        gamma: receptor turnover rate ([X-R] --> R)
        deg: universal degradation rate for free proteins
        """
        # RTotal: Total receptor number (nM)
        # self.RTotal = 2.7 # nM  0.0375 – 2.7 nM (corresponding to 18 –1300 molecules/µm2 and a 800 µm Matrigel layer on the cell surface)
        # k (nM-1*s-1) and r (s-1) of A+A' <-> AA'
        self.k_AC = 1e-4   # /nM/min
        self.r_AC = 1e-4  # /min
        # k (nM-1*s-1) and r (s-1) of B+A' <-> BA'
        self.k_BC = 1e-4  # nM/smin
        self.r_BC = 1e-3  # /min
        # k (nM-1*s-1) and r (s-1) of A+Receptor <-> A-Receptor complex
        self.k_AR = 4.5e-4  # nM/min
        self.r_AR = 1e-3  # /s
        # k (nM-1*s-1) and r (s-1) of B+Receptor <-> B-Receptor complex
        self.k_BR = 4.5e-4  # nM/min
        self.r_BR = 1e-3  # /min
        # Recycling rate (s-1) of A-Receptor and B-Receptor
        self.gamma = 4e-4

        self.ratio = self.r_BC / self.r_AC

        # Not considering basal degradation of species
        self.deg = 2e-5

        # Hill function
        self.k = 1
        self.n = 3

        # Put in params that were specified in input
        for entry in kwargs:
            setattr(self, entry, kwargs[entry])
