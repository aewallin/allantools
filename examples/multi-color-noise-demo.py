#!/usr/bin/env python


import matplotlib.pyplot as plt
import allantools as at
# note that latex requires e.g. a texlive installation and dvipng
# plt.rc('text', usetex=True)  # for latex
# plt.rc('font', family='serif')

"""

multi-colored noise example startinng from the ADEV


"""


def main():

    # Example Usage (Orolia RAFS)
    taus =[10, 30, 120, 480, 1.92e3, 1.54e4, 6.14e4, 1.23e5, 2.46e5, 4.92e5, 9.83e5 ]
    adevs = [1e-12, 3.99e-13, 1.66e-13, 8.85e-14, 4.15e-14, 1.62e-14, 8.99e-15, 9.73e-15, 1.17e-14, 1.26e-14, 1.36e-14]

    # Example Usage (Accubeat USO)
    #taus = [4, 8, 16, 32, 64, 128, 256, 512, 1.02e3, 2.05e3, 4.1e3, 8.19e3, 1.64e4, 3.28e4, 6.55e4]
    #adevs = [1.07e-13, 1.04e-13, 9.84e-14, 1.04e-13, 1.11e-13, 1.21e-13, 1.33e-13, 1.49e-13, 1.6e-13, 2.17e-13, 3.4e-13, 6.6e-13, 1.12e-12, 1.97e-12, 3.27e-12]
    
    # 1) ADEV -> PSD parameters
    f_nodes, Sy_nodes, h, alpha = at.adev2psd_piecewise_approx(adevs, taus, vartype="adev")
    print("adev2psd_piecewise_approx() done")
    
    # 2) PSD -> ADEV
    adev_from_psd = at.psd_piecewise_to_adev(h, alpha, f_nodes, taus)
    print("psd_piecewise_to_adev() done")
    
    # 3) PSD -> noise 
    duration = taus[-1] * 100
    timestep = taus[0] / 2
    noise = at.noise.timmer_koenig_from_psd(f_nodes, h, alpha, duration, timestep, output='phase', seed=1)
    print("timmer_koening noise generation done")
    
    # 4) noise -> ADEV
    taus_noise, adev_noise, _, _ = at.adev(noise, rate=1/timestep, data_type='phase', taus=taus)
    
    # Plot
    plt.figure(figsize=(10, 6))
    plt.loglog(taus, adevs, label="Original ADEV")
    plt.loglog(taus, adev_from_psd, label="ADEV -> PSD -> ADEV")
    plt.loglog(taus_noise, adev_noise, label="ADEV -> PSD -> noise -> ADEV")
    plt.grid(which="both")
    plt.ylabel("ADEV [-]")
    plt.xlabel(r"$\tau$ [s]") 
    plt.legend()
    plt.show()



if __name__ == "__main__":
    main()
