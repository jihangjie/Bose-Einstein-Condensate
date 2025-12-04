***************************************************************
************ Explanation of the supplemental files ************
***************************************************************


The supplemental folder contains numerical data and source code that we use to generate the numerical data. In particular, the whole supplemental package contains:


Folder "helper"
    * Contains helper functions for all c code to simulate the PDE
    * All c code in Folders "figure1", "figure2", "figure3", "figure4ab", and "figure4cd" needs the files saved in "helper"


Folder "figure1" (The solution dynamics starting from a Gaussian initial profile that evolves into a finite time singularity near x = 0.)
    * "BE_fig1.c": source code to simulate data used in Figure 1
    * "out_pde.dat": numerical data produced by BE_fig1.c and used to plot Figure 1
    * "out1_pde.dat": additional numerical data produced by BE_fig1.c 
    * "readme.txt": doc file to explain the format of data files


Folder "figure2" (The solution dynamics starting from a smoothed piecewise-constant initial profile with a dip in the middle that evolves into a finite time singularity near x = 0.)
    * "BE_fig2.c": source code to simulate data used in Figure 2
    * "out_pde.dat": numerical data produced by BE_fig2.c and used to plot Figure 2
    * "out1_pde.dat": additional numerical data produced by BE_fig2.c 
    * "readme.txt": doc file to explain the format of data files

Folder "figure3" (The solution dynamics starting from a mixture of two Gaussian-shaped humps that evolves into ridge moving toward x = 0.)
    * "BE_fig3.c": source code to simulate data used in Figure 3
    * "out_pde.dat": numerical data produced by BE_fig3.c and used to plot Figure 3
    * "out1_pde.dat": additional numerical data produced by BE_fig3.c 
    * "readme.txt": doc file to explain the format of data files

Folder "figure4ab" (The solution dynamics starting from a mixture of two Gaussian-shaped humps that evolves into ridge moving toward x = 0. The parameters n = 2, alpha = 8, and beta = 0.5)
    * "BE_fig4ab.c": source code to simulate data used in Figure 4(a-b)
    * "out_pde.dat": numerical data produced by BE_fig4ab.c and used to plot Figure 4(a-b)
    * "out1_pde.dat": additional numerical data produced by BE_fig4ab.c 
    * "readme.txt": doc file to explain the format of data files

Folder "figure4cd" (The solution dynamics starting from a mixture of two Gaussian-shaped humps that evolves into ridge moving toward x = 0. The parameters n = 4, alpha = 8, and beta = 0.5)
    * "BE_fig4cd.c": source code to simulate data used in Figure 4(c-d)
    * "out_pde.dat": numerical data produced by BE_fig4cd.c and used to plot Figure 4(c-d)
    * "out1_pde.dat": additional numerical data produced by BE_fig4cd.c 
    * "readme.txt": doc file to explain the format of data files



