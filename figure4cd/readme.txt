BE_fig4cd.c: source code to generate data for figure 4(c-d)

* This simulates the PDE (2.2) - (2.4) starting from a mixture of two Gaussian-shaped humps that evolves into ridge moving toward x = 0.
* The parameters n = 4, alpha = 8, and beta = 0.5
* This file depends on all functions in the folder helper/


out_pde.dat: numerical data produced by BE_fig4cd.c and used to plot figure 4(c-d) structured in data blocks separated by two empty lines, where the i-th data block contains the following information

* The 1st-5th columns contain information at time t_i in the order of <x>, <u>, <t>, <umax>, <xmax>
	1st column of data: values of <x> (grid point positions) at time t0 
	2nd column of data: values of <u> (u(x)) at time t0
	3rd column of data: values of <t> (time) at time t0
	4th column of data: values of <umax> (max_x(u(x))) at time t0
	5th column of data: values of <xmax> (x coordinate where umax is attained) at time t0

out1_pde.dat: numerical data produced by BE_fig4cd.c and used to plot figure 4(c-d) containing the following information

* The 1st-9th columns contain information in the order of <t>, <umax>, <xmax>, <uxxmax>, <xmin>, <umin>, <mass>, <contactLineLocation_right>, <contactLineLocation_left>

	1st column of data: values of <t> (time)
	2nd column of data: values of <umax> (max_x(u(x)) at time t) 
	3rd column of data: values of <xmax> (x coordinate where the maximum of u at time t is attained)
	4th column of data: values of <uxxmax> (the second-order derivative u_xx at xmax)
	5th column of data: values of <xmin> (x coordinate where the minimum of u at time t is attained)
	6th column of data: values of <umin> (min_x(u(x)) at time t) 
	7th column of data: values of <mass> (mass of solution at time t)
	8th column of data: values of <contactLineLocation_right> (x coordinate of the right contact line of the u profile at time t)
	9th column of data: values of <contactLineLocation_left> (x coordinate of the left contact line of the u profile at time t)

	
	
	
