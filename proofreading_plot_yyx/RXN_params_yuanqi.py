class RXN_params_yuanqi(object):
    """
    Container for reaction parameters
    """

    def __init__(self, **kwargs):
        """
        k_kin: kinase rate on membrane boundary(x=0)
        k_p: phosphotase rate in cytosol
        k_CAp: association rate for Ap and A'
        r_CAp: dissociation rate for A'-Ap complex
        k_CBp: association rate for Ap and A'
        r_CBp: dissociation rate for A'-Bp complex
        """
        self.k_kin = 0.2  # nM/s
        self.k_p = 5  # /s
        self.k_CAp = 0.1  # /nM/s = 10^5 /M/s
        self.r_CAp = 0.1  # /s
        self.k_CBp = 0.1  # /nM/s = 10^5 /M/s
        self.r_CBp = 1  # /s

        self.ratio = self.r_CBp / self.r_CAp

        # Put in params that were specified in input
        for entry in kwargs:
            setattr(self, entry, kwargs[entry])
