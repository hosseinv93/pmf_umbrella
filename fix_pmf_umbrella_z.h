/* -*- c++ -*- ----------------------------------------------------------
   Custom LAMMPS fix: pmf/umbrella/z

   Harmonic umbrella windows along the z-coordinate of a group COM.
   The fix applies a mass-weighted harmonic restraint to the group COM,
   automatically advances through a sequence of z0 windows, and writes
   per-window histograms and summary statistics for WHAM/MBAR/umbrella
   integration analysis.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(pmf/umbrella/z,FixPMFUmbrellaZ);
// clang-format on
#else

#ifndef LMP_FIX_PMF_UMBRELLA_Z_H
#define LMP_FIX_PMF_UMBRELLA_Z_H

#include "fix.h"

namespace LAMMPS_NS {

class FixPMFUmbrellaZ : public Fix {
 public:
  FixPMFUmbrellaZ(class LAMMPS *, int, char **);
  ~FixPMFUmbrellaZ() override;

  int setmask() override;
  void init() override;
  void setup(int) override;
  void post_force(int) override;
  double compute_scalar() override;
  double compute_vector(int) override;

 private:
  // input parameters
  double k_spring;
  double z_start, z_stop, dz_window;
  double dz_signed;
  bigint nequil, nprod, nevery;
  double bin_width;
  char *prefix;

  // derived quantities
  int nwindow;
  int nbins;
  double hist_zmin, hist_zmax;
  bigint steps_per_window;
  bigint startstep;
  double masstotal;

  // current state
  int current_window;
  int active_window;
  int phase;                 // 0 equilibration, 1 production, 2 finished
  double z0_current;
  double zcom_current;
  double f_spring_z;
  double f_system_z;
  double bias_energy;

  // accumulated production statistics for one window
  bigint nsamples;
  bigint out_of_range_low;
  bigint out_of_range_high;
  double sum_z;
  double sum_z2;
  double sum_fspring;
  double sum_fsystem;
  double sum_bias;
  double *hist;

  // output
  FILE *fsum;
  FILE *fhist;

  void reset_accumulators();
  void write_window(int);
  double target_z(int) const;
  int bin_index(double) const;
  void sample_current_state();
};

}    // namespace LAMMPS_NS

#endif
#endif
