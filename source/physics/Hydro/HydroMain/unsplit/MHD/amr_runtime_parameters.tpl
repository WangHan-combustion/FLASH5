## Lines starting with ## are comments inside template file
## All other lines including empty lines are non-comments
## 
## This file is a template for generating the amr_runtime_parameters
## file, which is used by Paramesh4.0 at runtime if compiled in 
## "ParameshLibraryMode", i.e., with LIBRARY defined, and by Paramesh4dev.
## For syntax of this file see "Readme.template"
##
## You probably do not want to change any lines in this files that
## already express values by using %(VARNAME)s syntax.
## Other entries may need manual adjustment for some simulations,
## in accordance with changes that would be made to paramesh_preprocessor.fh
## if not using ParameshLibraryMode.
##
##
## VALID VARIABLE NAMES FOR THIS TEMPLATE
##
##
## nguard            -> # of guard cells
## fixedblocksize    -> 0/1 (are we using fixed block size?) -- probably not useful here
## npg               -> F/T (are we using NO_PERMANENT_GUARDCELLS mode?)
## curvilinear       -> F/T (Enable Paramesh4 curvilinear code?)
## cartesian         -> F/T (Was the setup command asked to configure for "cartesian" geometry?)
## noncartesian      -> F/T (Was the setup command asked to configure for a non-Cartesian geometry?)
## cartesian_eff     -> F/T (Will run be effectively Cartesian, based on setup? - The negation of noncartesian.)
## ndim              -> # of dimensions (user specified)
## nxb,nyb,nzb       -> zones per block in each dimension (user specified)
## maxblocks         -> User specified value for MAXBLOCKS
## nflux             -> # of fluxes excluding those for named species and mass scalars -- probably not useful here
## nfacevars         -> # of face variables
## nvars             -> # of declared named variables -- probably not useful here
## nspecies          -> # of species -- probably not useful here
## nmass_scalars     -> # of mass scalars -- probably not useful here
##
## nunk_vars         -> # of vars including species and mass scalars
## nfluxes           -> # of fluxes including species and mass scalars
## nboundaries       -> # of boundary regions where boundary conditions are to be applied, normally 2*ndim
##
## NO other variables may be used with %(VARNAME)s syntax
##
## For the meaning of these AMR runtime parameters, cf. Paramesh4's amr_set_runtime_parameters routine.
##
## NO empty lines after theses initial comments!
##
%(maxblocks|-1)s   , maxblocks
%(ndim|2)s   , ndim
           1 , l2p5d  --- adjusted manually for StaggeredMesh MHD
%(nxb|-1)s   , nxb
%(nyb|-1)s   , nyb
%(nzb|-1)s   , nzb
%(nunk_vars)s   , nvar (including any species and mass scalars)
%(nfacevars)s   , nfacevar
           0 , nvaredge
           0 , nvarcorn
           1 , nvar_work
%(nguard)s   , nguard
%(nguard)s   , nguard_work
%(nfluxes)s   , nfluxvar (including any species and mass scalars)
           1 , nedgevar1  --- adjusted manually for StaggeredMesh MHD
           0 , iface_off
           1 , mflags
           1 , nfield_divf  --- adjusted manually for StaggeredMesh MHD
%(nboundaries)s   , nboundaries
 T , diagonals
 F , amr_error_checking
%(npg)s , no_permanent_guardcells
 F , advance_all_levels
 F , force_consistency
%(noncartesian)s , consv_fluxes
%(cartesian)s , consv_flux_densities
%(cartesian)s , edge_value
%(noncartesian)s , edge_value_integ
 F , var_dt
 F , pred_corr
 F , empty_cells
 F , conserve
%(cartesian_eff)s , divergence_free  --- adjusted manually for StaggeredMesh MHD
%(curvilinear)s , curvilinear
%(curvilinear)s , curvilinear_conserve
%(cartesian)s , cartesian_pm
%(noncartesian)s , cylindrical_pm
 F , spherical_pm
 F , polar_pm
 F , lsingular_line
 F , timing_mpi
 F , timing_mpix
'./' , output_dir (for log files etc. generated by Paramesh)
 ===== END OF RUNTIME PARAMETERS FROM TEMPLATE =====
