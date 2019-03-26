=====================================================
MemCurv
=====================================================

:Authors:       Ramachandra M. Bhaskara
:Year:          2019
:Licence:       GPLv3
:Copyright:     2019 Ramachandra M. Bhaskara, Gerhard Hummer
:References:    [Bhaskara et al, 2019] Bhaskara RM, Grumati P, Garcia-Padro J, Kalayil S, Covurrbais-Pinto A, Chen W, Kudrashyev M, Dikic I, Hummer G. Curvature induction and sensing by FAM134B Reticulon homology domain assist selective-ER phagy.

Description
=============

The MemCurv software package contains three python programs used to describe membrane shapes and quantify their curvature properties from MD smulations using the coarse-grained MARTINI model. The individula python programs

1. Compute Mean (Constant) curvature of the discontinuous membrane system by fitting membrane shape to spherical surfaces.
2. a. Fit shape profiles of curved or buckled membranes under periodic boundary conditions to obtain optimized Fourier representation of the curved surface.
   b. Compute local curvature properties of protein inclusions in curved or buckled membranes from optimized height profiles.

For a detailed description of the procedures and the algorithm, we refer to [Bhaskara et al, 2019].


Overview of the software
------------------------

This software contains two Python packages:

* "Fit_bicelle_sph_cap.py" provides algorithm for description of shape changes of a discontinuous membrane pathch. For.e.g. bicelle-to-vesicle conversion by patch-closure. The program uses least square fits for curved membrane-patches to spherical cup/cap shaped surfaces, to compute the constant mean curvature. 

.. image:: fam_bicelle.png

* "Get_mem_profile_hcf_abcd.py" and "Cal_local_curv_props_abcd.py" are part of a combined package which provide algorithms for description of curved membrane surfaces under periodic boundary conditions. They can be used independently with a trajectory of curved membrane to quantify the membrane shapes and to compute local curvature properties of protein inclusion. 

.. image:: fam_buckle.png
Help
====

Please, if you have an issue with the software, open an issue here on the github repository https://github.com/bio-phys/MemCurv/issues .

Dependencies and Software Requirements
=========================================

* Python 2.7
* Python packages: numpy, scipy, MDAnalysis, pandas

Installation
============
No instalation required. Directly run the Python programs if the dependecies are installed before.

Usage
=====

1.      Example input files from MD simulations of a discontinuous membrane system using MARTINI force-field, for e.g. bicelle/bilayer patches (GROMACS formatted .gro and .xtc).

* fam_bicelle.gro, a configuration file contining FAM134B embedded in bicelle (DMPC+DHPC)
* fam_bicelle.xtc, a processed GROMACS trajectory of the bicelle.

.. code-block:: bash 

        python Fit_bicelle_sph_cap.py

The output file contains

2.      Example input files from MD simulations of a curved membrane under periodic boundary conditions using MARTINI force-field (GROMACS formatted .gro and .xtc).

* fam_buckle.gro, a configuration file of the curved bilayer with protein (FAM134B) inclusion.
* fam_buckle_100ns.xtc, a processed GROMACS trajectory of the curved membrane simulation.

.. code-block:: bash

        python Get_mem_profile_hcf_abcd.py
        python Cal_local_curv_props_abcd.py

Other options and settings
--------------------------
Depending on membrane composition, the selection of atom type constituting the mid-plane of the bilayer needs to be selected. 

In curvature computations, If protein coordinates is provided, the local curvature is computed at the centre of mass of the protein atoms. 

Default settings
----------------

Output
------
Three output files are generated with default names.

Spherical surfaces of discontinuous membrane systems. 
-----------------------------------------
"Fit_bicelle_sph_cap.py" generates an output file "bicelle_curv_ts.dat" which contains the following tab delimited data.

======  ======
Column  Description
======  ======
1       Trajectory frame / Time [ns]
2       Radius of curvature, R_c [Å]
3       Mean curvature, H [1/Å]
4-6     Sphere center coordinates [x_c, y_c, z_c]
7       Residual (least square fit)
8       Iteration
======  ======

Shape profile of curved bilayer under PBC. 
-------------------------------------------
"Get_mem_profile_hcf_abcd.py" generates an output file "popf_1ns_k3_abcd_fam_buckled.dat" which contains optimized height coefficients describing the membrane shape profile. 

======  ======
Column  Description
======  ======
1       Trajectory frame / Time [ns]
2-last  Height coefficients. They can be written as four real valued (k x k) matrices.
======  ======

Local curvarure of protein inclusion in curved bilayer under PBC.
---------------------------------------------------------------------------
"Cal_local_curv_props_abcd.py" generates an output file "fam_1ns_abcd_k3_curv_props.dat" which contains the local curvature properties of the sampled protein curvatures (tab delimited file).

======  ======
Column  Description
======  ======
1       Trajectory frame / Time [ns]
2-7     Local curvature properties at protein center-of-mass position (COM)
2       Local Gaussian curvature, K_G(x,y)[1/Å^2]
3       Local Mean curvature, H(x,y) [1/Å]
4       Local Principal directional curvature, k_1(x,y) [1/Å]
5       Local Principal directional curvature, k_2(x,y) [1/Å]
6       Angle (deg) between local k_1(x,y) and e_x
7       Angle (deg) between local k_2(x,y) and e_x
8-13    Local curvature properties at amphipathic helix-1 COM (AH_1)
14-19   Local curvature properties at amphipathic helix 2 COM (AH_2)
20      Angle (deg) between AH_1 and e_x
21      Angle (deg) between AH_2 and e_x
22      Angle (deg) between AH_1 and AH_2
23      Angle (deg) between AH_1 and k_1(x,y)
24      Angle (deg) between AH_2 and k_1(x,y)
25-26   Protein COM position along the membrane [Px, Py]
======  ======

FAQs
====
Q: How is the sign of the bicelle curvature decided?

A: The sign of bicelle curvature depends on the identity of leaflets. Lipid flip-flop at the open edge merges the two leaflet, therefore, we use the identity as defined at the start frame. If more than 50 % of upper leaflet lipids remain above the fitted mid-plane spherical surface, then the curvature is positive, else the sign of the curvature is negative.

Q: My curvature values are negative instead of positive. What is going on?

A: The sign of the curvature is dictated by directionality-convention. In the Monge representation, using h(x,y) representation of the membrane profile gives the positive curvature for the lower leaflet and negative curvature for the upper leaflet. To compute curvatues along upper leaflet, as in the manuscript, we compute the negative of the Shape operator matrix (S=-S).
