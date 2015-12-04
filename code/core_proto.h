#include "core_allvars.h"

size_t myfread(void  *ptr,  size_t  size,  size_t  nmemb,  FILE *stream);
size_t myfwrite(void  *ptr,  size_t  size,  size_t  nmemb,  FILE *stream);
int myfseek(FILE *stream, long offset, int whence);

void construct_galaxies(int halonr, int tree);
void evolve_galaxies(int halonr, int ngal, int tree);
int  join_galaxies_of_progenitors(int halonr, int nstart);
void init(void);
void set_units(void);

void load_tree_table(int filenr);
void load_tree(int filenr, int nr);
void save_galaxies(int filenr, int tree);

void prepare_galaxy_for_output(int filenr, int tree, struct GALAXY *g, struct GALAXY_OUTPUT *o);

void free_galaxies_and_tree(void);
void free_tree_table(void);
void print_allocated(void);

void read_parameter_file(char *fname);
void *mymalloc(size_t n);
void myfree(void *p);
void myexit(int signum);

void finalize_galaxy_file(int filenr);

void starformation_and_feedback(int p, int centralgal, double time, double dt, int halonr, int step);
void add_galaxies_together(int t, int p, double mass_ratio, double *disc_mass_ratio, int centralgal, double dt, double *PostRetroGas);
void init_galaxy(int p, int halonr);
double infall_recipe(int centralgal, int ngal, double Zcurr);
void add_infall_to_hot(int centralgal, double infallingGas);
double cooling_recipe(int centralgal, double dt);
void cool_gas_onto_galaxy(int p, int centralgal, double coolingGas, double dt, int step);
void reincorporate_gas(int centralgal, double dt);
double estimate_merging_time(int prog, int mother_halo, int ngal);
void deal_with_galaxy_merger(int p, int merger_centralgal, int centralgal, double time, double dt, int halonr, int step);
double dmax(double x, double y);
double do_reionization(int centralgal, double Zcurr);
double do_AGN_heating(double coolingGas, int centralgal, double dt, double x, double rcool);
void collisional_starburst_recipe(double disc_mass_ratio[N_BINS], int merger_centralgal, int centralgal, double time, double dt, int halonr, int mode, int step, double mass_ratio);
void update_from_star_formation(int p, double stars, double metallicity, int i);
void update_from_feedback(int p, int centralgal, double reheated_mass, double metallicity, int i);
void update_from_ejection(int centralgal, double ejected_mass);
void stars_to_bulge(int t, int p);
double grow_black_hole(int merger_centralgal, double* disc_mass_ratio);
void check_disk_instability(int p, int centralgal, double time, double dt, int step);

// Adam's functions (most of the others above have been heavily modified too)
void retrograde_gas_collision(int p, double RetroGas[N_BINS], double cos_angle_halo_new, double cos_angle_disc_new, double J_disc, double J_cool);
void combine_stellar_discs(int p, double NewStars[N_BINS], double NewStarsMetals[N_BINS]);
void project_disc(double DiscMass[N_BINS], double cos_angle, int p, double *NewDisc);
double deal_with_unstable_gas(double unstable_gas, int p, int i, double V_rot, double metallicity, int centralgal, int direct_to_BH, double r_inner, double r_outer);
double get_disc_gas(int halonr);
double get_disc_stars(int halonr);
double get_disc_ang_mom(int halonr, int type);
void precess_gas(int p, double dt, int halonr);
void update_HI_H2(int p);
void update_disc_radii(int p);
//

void strip_from_satellite(int halonr, int centralgal, int gal);
void disrupt_satellite_to_ICS(int centralgal, int gal);
void quasar_mode_wind(int gal, float BHaccrete);

double get_metallicity(double gas, double metals);
double get_virial_velocity(int halonr);
double get_virial_radius(int halonr);
double get_virial_mass(int halonr);
double get_disk_radius(int halonr, int p);


void read_output_snaps(void);
void read_snap_list(void);
void read_cooling_functions(void);
double get_metaldependent_cooling_rate(double logTemp, double logZ);
double get_rate(int tab, double logTemp);

double time_to_present(double z);
double integrand_time_to_present(double a, void *param);

double metallicity_dependent_star_formation(int p);
double Z_dependent_SF(float lower_limit, float upper_limit, float Sigma_c0, float Xi, float gamma);
double integrand_Z_dependent_SF(double q, void *p);




