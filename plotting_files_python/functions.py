


# Importing
import numpy as np
import bisect
import random
import time
import matplotlib.pyplot as plt
from scipy.optimize import fsolve


### THIS FILE DEFINES FUNCTIONS

start = time.time()


# Function to compute final speed
def simulate_speed(
        N: int,
        lam: float,
        omega: float,
        T: float
) -> float:
    """
        Simulates the speed for the finite N miner Blockchain CTMC model.

        Parameters
        ----------
        N : int
            Number of miners
        lam : float
            Mining rate
        omega : float
            Communication rate
        T : float
            Final time of simulation

        Returns
        -------
        the simulated estimated speed defined by v_N = (X_T - X_T0)/(T - T0),
        where:
         - T0 = 0.5T (to measure average speed between T0 and T, removes bias for early time)
         - X_T = position of front at time T
        """

    levels=[0]
    counts=[N]

    t=0
    mean=0

    T0 = round(0.5 * T)
    mean_T0 = 0
    T0_check = -1

    while t<T:

        n=np.array(counts)


        Rb=lam*N
        Ri=omega/(2*N)*(N**2-n.dot(n))
        R=Rb+Ri

        t+=np.random.exponential(1/R)

        if random.random()<Rb/R:

            # birth
            i=random.choices(range(len(levels)),weights=counts)[0]
            k=levels[i]

            counts[i]-=1

            if counts[i]==0:
                del levels[i]
                del counts[i]

            pos=bisect.bisect_left(levels,k+1)

            if pos<len(levels) and levels[pos]==k+1:
                counts[pos]+=1
            else:
                levels.insert(pos,k+1)
                counts.insert(pos,1)

            mean+=1/N

        else:

            n=np.array(counts)
            C=np.cumsum(n)

            weights=n*np.concatenate(([0],C[:-1]))

            j=random.choices(range(len(levels)),weights=weights)[0]
            ell=levels[j]

            i=random.choices(range(j),weights=n[:j])[0]
            k=levels[i]

            counts[i]-=1
            counts[j]+=1

            if counts[i]==0:
                del levels[i]
                del counts[i]

            mean+=(ell-k)/N

        if t>T0 and T0_check == -1:

            mean_T0 = mean
            T0_check = 1

    return (mean-mean_T0)/(t-T0)



# Function to compute the Derrida speed prediction v0
def predict_speed(lam, omega):
    def f(g0, lam, omega):
        return lam * (np.exp(g0) - 1) + omega - lam * g0 * np.exp(g0)

    g0_initial = 1.0
    g0_solution, infodict, ier, mesg = fsolve(
        lambda x: f(x, lam, omega), g0_initial, full_output=True
    )
    g0 = g0_solution[0]
    v0 = lam * np.exp(g0)
    return g0, v0
