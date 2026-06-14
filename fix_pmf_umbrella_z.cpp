// clang-format off
/* ----------------------------------------------------------------------
   Custom LAMMPS fix: pmf/umbrella/z

   This fix applies harmonic umbrella windows to the center of mass of
   the fix group along z:

       U_bias(z) = 0.5*K*(z_COM - z0)^2
       F_spring,z = -K*(z_COM - z0)

   It advances automatically through a list of window centers from zstart
   to zstop with spacing dz, and writes histograms for WHAM/MBAR plus
   summary mean-force quantities.
------------------------------------------------------------------------- */

#include "fix_pmf_umbrella_z.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "memory.h"
#include "respa.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixPMFUmbrellaZ::FixPMFUmbrellaZ(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg), prefix(nullptr), hist(nullptr), fsum(nullptr), fhist(nullptr)
{
  // Syntax:
  // fix ID group-ID pmf/umbrella/z K zstart zstop dz nequil nprod nevery binwidth file prefix
  if (narg != 13) error->all(FLERR,"Illegal fix pmf/umbrella/z command");

  k_spring = utils::numeric(FLERR,arg[3],false,lmp);
  z_start = utils::numeric(FLERR,arg[4],false,lmp);
  z_stop = utils::numeric(FLERR,arg[5],false,lmp);
  dz_window = utils::numeric(FLERR,arg[6],false,lmp);
  nequil = utils::bnumeric(FLERR,arg[7],false,lmp);
  nprod = utils::bnumeric(FLERR,arg[8],false,lmp);
  nevery = utils::bnumeric(FLERR,arg[9],false,lmp);
  bin_width = utils::numeric(FLERR,arg[10],false,lmp);

  if (strcmp(arg[11],"file") != 0)
    error->all(FLERR,"Illegal fix pmf/umbrella/z command: expected keyword 'file'");

  prefix = utils::strdup(arg[12]);

  if (k_spring <= 0.0) error->all(FLERR,"K must be > 0 for fix pmf/umbrella/z");
  if (dz_window <= 0.0) error->all(FLERR,"dz must be > 0 for fix pmf/umbrella/z");
  if (nequil < 0 || nprod <= 0 || nevery <= 0)
    error->all(FLERR,"Bad nequil/nprod/nevery for fix pmf/umbrella/z");
  if (bin_width <= 0.0) error->all(FLERR,"binwidth must be > 0 for fix pmf/umbrella/z");

  dz_signed = (z_stop >= z_start) ? dz_window : -dz_window;
  nwindow = static_cast<int>(std::floor(std::fabs(z_stop-z_start)/dz_window + 0.5)) + 1;
  if (nwindow < 1) nwindow = 1;
  steps_per_window = nequil + nprod;

  // histogram range: one additional window spacing outside the scan interval
  hist_zmin = MIN(z_start,z_stop) - dz_window;
  hist_zmax = MAX(z_start,z_stop) + dz_window;
  nbins = static_cast<int>(std::ceil((hist_zmax-hist_zmin)/bin_width));
  if (nbins < 1) error->all(FLERR,"Bad histogram range for fix pmf/umbrella/z");
  hist_zmax = hist_zmin + nbins*bin_width;

  memory->create(hist,nbins,"pmf/umbrella/z:hist");
  reset_accumulators();

  // This fix contributes potential energy via compute_scalar() if enabled by fix_modify energy yes.
  scalar_flag = 1;
  global_freq = 1;
  extscalar = 0;
  energy_global_flag = 1;

  // Vector output:
  // 1 z0, 2 z_COM, 3 F_spring_z, 4 F_system_z, 5 U_bias, 6 window, 7 phase, 8 nsamples
  vector_flag = 1;
  size_vector = 8;
  extvector = 0;

  restart_global = 0;
  virial_global_flag = virial_peratom_flag = 0;

  current_window = -1;
  active_window = -1;
  phase = 0;
  z0_current = z_start;
  zcom_current = 0.0;
  f_spring_z = 0.0;
  f_system_z = 0.0;
  bias_energy = 0.0;

  if (comm->me == 0) {
    char fname1[1024], fname2[1024];
    snprintf(fname1,1024,"%s.summary.dat",prefix);
    snprintf(fname2,1024,"%s.hist.dat",prefix);
    fsum = fopen(fname1,"w");
    fhist = fopen(fname2,"w");
    if (fsum == nullptr || fhist == nullptr)
      error->one(FLERR,"Could not open output files for fix pmf/umbrella/z");

    fprintf(fsum,"# fix pmf/umbrella/z summary\n");
    fprintf(fsum,"# K %.16g zstart %.16g zstop %.16g dz %.16g nequil %lld nprod %lld nevery %lld binwidth %.16g\n",
            k_spring,z_start,z_stop,dz_window,(long long)nequil,(long long)nprod,(long long)nevery,bin_width);
    fprintf(fsum,"# columns: window z0 nsamples <z_COM> std_z <F_spring_z> <F_system_z> <U_bias> out_low out_high\n");
    fflush(fsum);

    fprintf(fhist,"# fix pmf/umbrella/z histograms for WHAM/MBAR\n");
    fprintf(fhist,"# K %.16g zmin %.16g zmax %.16g binwidth %.16g nbins %d\n",
            k_spring,hist_zmin,hist_zmax,bin_width,nbins);
    fprintf(fhist,"# columns: window z0 K bin_center count\n");
    fflush(fhist);
  }
}

/* ---------------------------------------------------------------------- */

FixPMFUmbrellaZ::~FixPMFUmbrellaZ()
{
  // write last active window if the fix is destroyed before switching windows
  if (active_window >= 0 && active_window < nwindow && nsamples > 0) write_window(active_window);

  if (comm->me == 0) {
    if (fsum) fclose(fsum);
    if (fhist) fclose(fhist);
  }
  memory->destroy(hist);
  delete [] prefix;
}

/* ---------------------------------------------------------------------- */

int FixPMFUmbrellaZ::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixPMFUmbrellaZ::init()
{
  masstotal = group->mass(igroup);
  if (masstotal <= 0.0) error->all(FLERR,"Group mass is zero in fix pmf/umbrella/z");
}

/* ---------------------------------------------------------------------- */

void FixPMFUmbrellaZ::setup(int vflag)
{
  startstep = update->ntimestep;
  current_window = -1;
  active_window = -1;
  phase = 0;
  reset_accumulators();
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixPMFUmbrellaZ::post_force(int /*vflag*/)
{
  bigint elapsed = update->ntimestep - startstep;
  if (elapsed < 0) elapsed = 0;

  int win = static_cast<int>(elapsed / steps_per_window);
  bigint local_step = elapsed - static_cast<bigint>(win)*steps_per_window;

  if (win >= nwindow) {
    if (active_window >= 0 && active_window < nwindow && nsamples > 0) {
      write_window(active_window);
      reset_accumulators();
    }
    current_window = nwindow;
    active_window = -1;
    phase = 2;
    f_spring_z = f_system_z = bias_energy = 0.0;
    return;
  }

  if (win != active_window) {
    if (active_window >= 0 && active_window < nwindow && nsamples > 0) write_window(active_window);
    reset_accumulators();
    active_window = win;
    current_window = win;
  }

  phase = (local_step < nequil) ? 0 : 1;
  z0_current = target_z(win);

  // group COM
  double xcm[3];
  group->xcm(igroup,masstotal,xcm);
  zcom_current = xcm[2];

  // harmonic force on the group COM
  double dz = zcom_current - z0_current;
  f_spring_z = -k_spring * dz;
  f_system_z = -f_spring_z;
  bias_energy = 0.5 * k_spring * dz * dz;

  // distribute force to atoms by mass fraction
  double **f = atom->f;
  int *mask = atom->mask;
  int *type = atom->type;
  double *mass = atom->mass;
  double *rmass = atom->rmass;
  int nlocal = atom->nlocal;
  double massfrac;

  if (rmass) {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        massfrac = rmass[i] / masstotal;
        f[i][2] += f_spring_z * massfrac;
      }
    }
  } else {
    for (int i = 0; i < nlocal; i++) {
      if (mask[i] & groupbit) {
        massfrac = mass[type[i]] / masstotal;
        f[i][2] += f_spring_z * massfrac;
      }
    }
  }

  if (phase == 1) {
    bigint prod_step = local_step - nequil;
    if ((prod_step % nevery) == 0) sample_current_state();
  }
}

/* ---------------------------------------------------------------------- */

double FixPMFUmbrellaZ::compute_scalar()
{
  return bias_energy;
}

/* ---------------------------------------------------------------------- */

double FixPMFUmbrellaZ::compute_vector(int n)
{
  if (n == 0) return z0_current;
  if (n == 1) return zcom_current;
  if (n == 2) return f_spring_z;
  if (n == 3) return f_system_z;
  if (n == 4) return bias_energy;
  if (n == 5) return static_cast<double>(current_window+1);
  if (n == 6) return static_cast<double>(phase);
  if (n == 7) return static_cast<double>(nsamples);
  return 0.0;
}

/* ---------------------------------------------------------------------- */

void FixPMFUmbrellaZ::reset_accumulators()
{
  nsamples = 0;
  out_of_range_low = 0;
  out_of_range_high = 0;
  sum_z = sum_z2 = 0.0;
  sum_fspring = sum_fsystem = sum_bias = 0.0;
  for (int i = 0; i < nbins; i++) hist[i] = 0.0;
}

/* ---------------------------------------------------------------------- */

void FixPMFUmbrellaZ::write_window(int win)
{
  if (comm->me != 0 || nsamples <= 0) return;

  double z0 = target_z(win);
  double mean_z = sum_z / static_cast<double>(nsamples);
  double mean_z2 = sum_z2 / static_cast<double>(nsamples);
  double var_z = mean_z2 - mean_z*mean_z;
  if (var_z < 0.0) var_z = 0.0;
  double std_z = sqrt(var_z);
  double mean_fspring = sum_fspring / static_cast<double>(nsamples);
  double mean_fsystem = sum_fsystem / static_cast<double>(nsamples);
  double mean_bias = sum_bias / static_cast<double>(nsamples);

  fprintf(fsum,"%d %.16g %lld %.16g %.16g %.16g %.16g %.16g %lld %lld\n",
          win+1,z0,(long long)nsamples,mean_z,std_z,mean_fspring,mean_fsystem,mean_bias,
          (long long)out_of_range_low,(long long)out_of_range_high);
  fflush(fsum);

  for (int i = 0; i < nbins; i++) {
    double zbin = hist_zmin + (i+0.5)*bin_width;
    fprintf(fhist,"%d %.16g %.16g %.16g %.16g\n",win+1,z0,k_spring,zbin,hist[i]);
  }
  fflush(fhist);
}

/* ---------------------------------------------------------------------- */

double FixPMFUmbrellaZ::target_z(int win) const
{
  if (win <= 0) return z_start;
  if (win >= nwindow-1) return z_stop;
  return z_start + static_cast<double>(win)*dz_signed;
}

/* ---------------------------------------------------------------------- */

int FixPMFUmbrellaZ::bin_index(double z) const
{
  int ibin = static_cast<int>(std::floor((z - hist_zmin)/bin_width));
  return ibin;
}

/* ---------------------------------------------------------------------- */

void FixPMFUmbrellaZ::sample_current_state()
{
  nsamples++;
  sum_z += zcom_current;
  sum_z2 += zcom_current*zcom_current;
  sum_fspring += f_spring_z;
  sum_fsystem += f_system_z;
  sum_bias += bias_energy;

  int ibin = bin_index(zcom_current);
  if (ibin < 0) out_of_range_low++;
  else if (ibin >= nbins) out_of_range_high++;
  else hist[ibin] += 1.0;
}

/* ---------------------------------------------------------------------- */
