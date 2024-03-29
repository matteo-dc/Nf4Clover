#include "read_input.hpp"
#include "aliases.hpp"
#include "global.hpp"
#include "operations.hpp"
#include "read.hpp"

#define DEFAULT_STR_VAL "null"
#define DEFAULT_INT_VAL -1
#define DEFAULT_DOUBLE_VAL 1.2345

// define global variables
int nconfs, njacks, nr, ntypes, nhits, Nf, Nc, UseSigma1, UseEffMass, nbeta, only_basic,
    compute_mpcac, load_ave, load_chir, load, QCD_on_the_right, sub_boosted, sub_ptilde,
    subtraction;
int clust_size, nbil, combo;
vector<double> beta;
vector<int> nm_Sea;
int nm_Sea_max;
vector<vector<int>> SeaMasses_label; // SeaMasses_label[Nbeta][NSeaMass]
vector<vector<double>> SeaMasses;    // SeaMasses[Nbeta][NSeaMass]
int L, T;
vector<double> ainv;
vector<double> csw;
int conf_init, conf_step, nm, neq, neq2, nmr, delta_tmin, delta_tmax;
double kappa, mu_sea, plaquette, LambdaQCD, p2min, p2max, thresh, p2ref, stepfunc_min, stepfunc_max;
vector<double> mass_val;
vector<double> p2max_M3_M4, p2min_M3_M4;
string mom_path, action, path_folder, scheme, BC, out_hadr, analysis, clover, path_ensemble,
    an_suffix, chir_ansatz_val, chir_ansatz_sea;
string load_label;
vector<string> path_analysis;
vector<string> beta_label;   // beta_label[Nbeta]
vector<string> volume_label; // volume_label[Nbeta]
bool inte_analysis;
bool clover_analysis;

coords_t sizeV;

char tok[128];

TK_glb_t get_TK_glb(FILE *fin)
{
  // read a token
  int rc = fscanf(fin, "%s", tok);
  if (rc != 1)
  {
    if (feof(fin))
      return FEOF_GLB_TK;
    else
    {
      fprintf(stderr, "Getting %d while reading token\n", rc);
      exit(FAILED_READ);
    }
  }

  // parse the token
  if (strcasecmp(tok, nconfs_tag) == 0)
    return NCONFS_TK;
  if (strcasecmp(tok, njacks_tag) == 0)
    return NJACKS_TK;
  if (strcasecmp(tok, act_tag) == 0)
    return ACT_TK;
  if (strcasecmp(tok, path_folder_tag) == 0)
    return PATH_TK;
  if (strcasecmp(tok, nbeta_tag) == 0)
    return NBETA_TK;
  if (strcasecmp(tok, beta_tag) == 0)
    return BETA_TK;
  if (strcasecmp(tok, beta_label_tag) == 0)
    return BETA_LAB_TK;
  if (strcasecmp(tok, volume_label_tag) == 0)
    return VOLUME_LAB_TK;
  if (strcasecmp(tok, nm_Sea_tag) == 0)
    return NM_SEA_TK;
  if (strcasecmp(tok, SeaMasses_label_tag) == 0)
    return SEAMASSES_LAB_TK;
  if (strcasecmp(tok, SeaMasses_tag) == 0)
    return SEAMASSES_TK;
  if (strcasecmp(tok, nr_tag) == 0)
    return NR_TK;
  if (strcasecmp(tok, ntypes_tag) == 0)
    return NTYPES_TK;
  if (strcasecmp(tok, nhits_tag) == 0)
    return NHITS_TK;
  if (strcasecmp(tok, UseSigma1_tag) == 0)
    return USE_SIGMA_TK;
  if (strcasecmp(tok, UseEffMass_tag) == 0)
    return USE_EFF_MASS_TK;
  if (strcasecmp(tok, scheme_tag) == 0)
    return SCHEME_TK;
  if (strcasecmp(tok, Nc_tag) == 0)
    return NC_TK;
  if (strcasecmp(tok, Nf_tag) == 0)
    return NF_TK;
  if (strcasecmp(tok, ainv_tag) == 0)
    return AINV_TK;
  if (strcasecmp(tok, csw_tag) == 0)
    return CSW_TK;
  if (strcasecmp(tok, LambdaQCD_tag) == 0)
    return LAMBDAQCD_TK;
  if (strcasecmp(tok, BC_tag) == 0)
    return BC_TK;
  if (strcasecmp(tok, p2min_tag) == 0)
    return P2MIN_TK;
  if (strcasecmp(tok, p2max_tag) == 0)
    return P2MAX_TK;
  if (strcasecmp(tok, thresh_tag) == 0)
    return THRESH_TK;
  if (strcasecmp(tok, out_hadr_tag) == 0)
    return OUT_HADR_TK;
  if (strcasecmp(tok, only_basic_tag) == 0)
    return ONLY_BASIC_TK;
  if (strcasecmp(tok, compute_mpcac_tag) == 0)
    return COMPUTE_MPCAC_TK;
  if (strcasecmp(tok, analysis_tag) == 0)
    return ANALYSIS_TK;
  if (strcasecmp(tok, clover_tag) == 0)
    return CLOVER_TK;
  if (strcasecmp(tok, p2ref_tag) == 0)
    return P2REF_TK;
  if (strcasecmp(tok, load_ave_tag) == 0)
    return LOAD_AVE_TK;
  if (strcasecmp(tok, load_chir_tag) == 0)
    return LOAD_CHIR_TK;
  if (strcasecmp(tok, an_suffix_tag) == 0)
    return SUFFIX_TK;
  if (strcasecmp(tok, QCD_on_the_right_tag) == 0)
    return QCDONTHERIGHT_TK;
  if (strcasecmp(tok, sub_boosted_tag) == 0)
    return SUB_BOOSTED_TK;
  if (strcasecmp(tok, sub_ptilde_tag) == 0)
    return SUB_PTILDE_TK;
  if (strcasecmp(tok, subtraction_tag) == 0)
    return SUBTRACTION_TK;
  if (strcasecmp(tok, chir_ansatz_val_tag) == 0)
    return CHIR_ANSATZ_VAL_TK;
  if (strcasecmp(tok, chir_ansatz_sea_tag) == 0)
    return CHIR_ANSATZ_SEA_TK;
  if (strcasecmp(tok, stepfunc_min_tag) == 0)
    return STEPFUNC_MIN_TK;
  if (strcasecmp(tok, stepfunc_max_tag) == 0)
    return STEPFUNC_MAX_TK;
  if (strcasecmp(tok, p2max_M3_M4_tag) == 0)
    return P2MAX_M3_M4_TK;
  if (strcasecmp(tok, p2min_M3_M4_tag) == 0)
    return P2MIN_M3_M4_TK;
  if (strcasecmp(tok, load_label_tag) == 0)
    return LOAD_LABEL_TK;

  return VALUE_GLB_TK;
}

TK_t get_TK(FILE *fin)
{
  // read a token
  int rc = fscanf(fin, "%s", tok);
  if (rc != 1)
  {
    if (feof(fin))
      return FEOF_TK;
    else
    {
      fprintf(stderr, "Getting %d while reading token\n", rc);
      exit(FAILED_READ);
    }
  }

  // parse the token
  if (strcasecmp(tok, mom_list_tag) == 0)
    return MOM_LIST_TK;
  if (strcasecmp(tok, L_tag) == 0)
    return L_TK;
  if (strcasecmp(tok, T_tag) == 0)
    return T_TK;
  if (strcasecmp(tok, conf_init_tag) == 0)
    return CONF_INIT_TK;
  if (strcasecmp(tok, conf_step_tag) == 0)
    return CONF_STEP_TK;
  if (strcasecmp(tok, mu_sea_tag) == 0)
    return MU_SEA_TK;
  if (strcasecmp(tok, kappa_tag) == 0)
    return KAPPA_TK;
  if (strcasecmp(tok, nm_tag) == 0)
    return NM_VAL_TK;
  if (strcasecmp(tok, mass_val_tag) == 0)
    return VALMASSES_TK;
  if (strcasecmp(tok, plaquette_tag) == 0)
    return PLAQ_TK;
  if (strcasecmp(tok, delta_tmin_tag) == 0)
    return DELTA_TMIN_TK;
  if (strcasecmp(tok, delta_tmax_tag) == 0)
    return DELTA_TMAX_TK;

  return VALUE_TK;
}

// parse the value string (glb input)
template <class T>
void _get_value_glb(FILE *fin, T &ret, const char *t)
{
  TK_glb_t tk = get_TK_glb(fin);
  if (tk != VALUE_GLB_TK)
  {
    fprintf(stderr, "Getting token %s in the wrong place\n", tok);
    exit(MISPLACED_TK);
  }

  int rc = sscanf(tok, t, &ret);

  if (rc != 1)
  {
    fprintf(stderr, "Converting %s to %s failed\n", tok, t);
    exit(FAILED_CONVERSION);
  }
}

void get_value_glb(FILE *fin, double &out) { return _get_value_glb(fin, out, "%lg"); }

void get_value_glb(FILE *fin, int &out) { return _get_value_glb(fin, out, "%d"); }

void get_value_glb(FILE *fin, string &out)
{
  char temp[1024];
  _get_value_glb(fin, temp, "%s");
  out = string(temp);
}

// parse the value string
template <class T>
void _get_value(FILE *fin, T &ret, const char *t)
{
  TK_t tk = get_TK(fin);
  if (tk != VALUE_TK)
  {
    fprintf(stderr, "Getting token %s in the wrong place\n", tok);
    exit(MISPLACED_TK);
  }

  int rc = sscanf(tok, t, &ret);

  if (rc != 1)
  {
    fprintf(stderr, "Converting %s to %s failed\n", tok, t);
    exit(FAILED_CONVERSION);
  }
}

void get_value(FILE *fin, double &out) { return _get_value(fin, out, "%lg"); }

void get_value(FILE *fin, int &out) { return _get_value(fin, out, "%d"); }

void get_value(FILE *fin, string &out)
{
  char temp[1024];
  _get_value(fin, temp, "%s");
  out = string(temp);
}

// check
void check_str_par(const string str, const char *name)
{
  if (str.compare(DEFAULT_STR_VAL) == 0)
  {
    fprintf(stderr, "%s not initialized\n", name);
    exit(UNINITIALIZED_PAR);
  }
}

void check_int_par(const int val, const char *name)
{
  if (val == DEFAULT_INT_VAL)
  {
    fprintf(stderr, "%s not initialized\n", name);
    exit(UNINITIALIZED_PAR);
  }
}

void check_double_par(const double val, const char *name)
{
  if (val == DEFAULT_DOUBLE_VAL)
  {
    fprintf(stderr, "%s not initialized\n", name);
    exit(UNINITIALIZED_PAR);
  }
}

// factorial
int fact(int n)
{
  if (n > 1)
    return n * fact(n - 1);
  else
    return 1;
}

// print line of text
void printLine(char chLeft, char chMiddle, char chRight, int length)
{
  printf("%c", chLeft);
  for (int i = 0; i < length - 2; i++)
    printf("%c", chMiddle);
  printf("%c", chRight);
  printf("\n");
}

// print the box
void printBox(const char *content)
{
  int contentLength = strlen(content);
  int lineLength = 2 * contentLength + 2;

  printLine('*', '-', '*', lineLength);
  string padding(contentLength / 2, ' ');
  cout << "|" << padding << content << padding << "|" << endl;
  printLine('*', '-', '*', lineLength);
  cout << endl;
}

// reads the input file
void read_input_glb(const char path[])
{
  FILE *fin = fopen(path, "r");
  if (not fin)
  {
    fprintf(stderr, "Failed to open \"%s\"\n", path);
    exit(FAILED_OPEN);
  }

  action = DEFAULT_STR_VAL;
  scheme = DEFAULT_STR_VAL;
  path_folder = DEFAULT_STR_VAL;
  BC = DEFAULT_STR_VAL;
  nconfs = DEFAULT_INT_VAL;
  njacks = DEFAULT_INT_VAL;
  nr = DEFAULT_INT_VAL;
  ntypes = DEFAULT_INT_VAL;
  nhits = DEFAULT_INT_VAL;
  Nc = DEFAULT_INT_VAL;
  Nf = DEFAULT_INT_VAL;
  UseSigma1 = DEFAULT_INT_VAL;
  UseEffMass = DEFAULT_INT_VAL;
  LambdaQCD = DEFAULT_DOUBLE_VAL;
  p2min = DEFAULT_DOUBLE_VAL;
  p2max = DEFAULT_DOUBLE_VAL;
  thresh = DEFAULT_DOUBLE_VAL;
  out_hadr = DEFAULT_STR_VAL;
  only_basic = DEFAULT_INT_VAL;
  compute_mpcac = DEFAULT_INT_VAL;
  analysis = DEFAULT_STR_VAL;
  clover = DEFAULT_STR_VAL;
  p2ref = DEFAULT_DOUBLE_VAL;
  load_ave = DEFAULT_INT_VAL;
  load_chir = DEFAULT_INT_VAL;
  an_suffix = DEFAULT_STR_VAL;
  QCD_on_the_right = DEFAULT_INT_VAL;
  sub_boosted = DEFAULT_INT_VAL;
  sub_ptilde = DEFAULT_INT_VAL;
  subtraction = DEFAULT_INT_VAL;
  chir_ansatz_val = DEFAULT_STR_VAL;
  chir_ansatz_sea = DEFAULT_STR_VAL;
  stepfunc_min = DEFAULT_DOUBLE_VAL;
  stepfunc_max = DEFAULT_DOUBLE_VAL;
  load_label = DEFAULT_STR_VAL;

  while (not feof(fin))
  {
    TK_glb_t tk = get_TK_glb(fin);
    switch (tk)
    {
    case VALUE_GLB_TK:
      fprintf(stderr, "Invalid token %s found\n", tok);
      exit(1);
      break;
    case NCONFS_TK:
      get_value_glb(fin, nconfs);
      break;
    case NJACKS_TK:
      get_value_glb(fin, njacks);
      break;
    case BC_TK:
      get_value_glb(fin, BC);
      break;
    case NBETA_TK:
      get_value_glb(fin, nbeta);
      break;
    case BETA_TK:
      beta.resize(nbeta);
      for (int b = 0; b < nbeta; b++)
        get_value_glb(fin, beta[b]);
      break;
    case BETA_LAB_TK:
      beta_label.resize(nbeta);
      for (int b = 0; b < nbeta; b++)
        get_value_glb(fin, beta_label[b]);
      break;
    case VOLUME_LAB_TK:
      volume_label.resize(nbeta);
      for (int b = 0; b < nbeta; b++)
        get_value_glb(fin, volume_label[b]);
      break;
    case NM_SEA_TK:
      nm_Sea.resize(nbeta);
      for (int b = 0; b < nbeta; b++)
        get_value_glb(fin, nm_Sea[b]);
      break;
    case SEAMASSES_LAB_TK:
      SeaMasses_label.resize(nbeta);
      for (int b = 0; b < nbeta; b++)
      {
        SeaMasses_label[b].resize(nm_Sea[b]);
        for (int m = 0; m < nm_Sea[b]; m++)
          get_value_glb(fin, SeaMasses_label[b][m]);
      }
      break;
    case SEAMASSES_TK:
      SeaMasses.resize(nbeta);
      for (int b = 0; b < nbeta; b++)
      {
        SeaMasses[b].resize(nm_Sea[b]);
        for (int m = 0; m < nm_Sea[b]; m++)
          get_value_glb(fin, SeaMasses[b][m]);
      }
      break;
    case ACT_TK:
      get_value_glb(fin, action);
      break;
    case PATH_TK:
      get_value_glb(fin, path_folder);
      break;
    case NR_TK:
      get_value_glb(fin, nr);
      break;
    case NTYPES_TK:
      get_value_glb(fin, ntypes);
      break;
    case NHITS_TK:
      get_value_glb(fin, nhits);
      break;
    case USE_SIGMA_TK:
      get_value_glb(fin, UseSigma1);
      break;
    case USE_EFF_MASS_TK:
      get_value_glb(fin, UseEffMass);
      break;
    case SCHEME_TK:
      get_value_glb(fin, scheme);
      break;
    case NC_TK:
      get_value_glb(fin, Nc);
      break;
    case NF_TK:
      get_value_glb(fin, Nf);
      break;
    case AINV_TK:
      ainv.resize(nbeta);
      for (int b = 0; b < nbeta; b++)
        get_value_glb(fin, ainv[b]);
      break;
    case CSW_TK:
      csw.resize(nbeta);
      for (int b = 0; b < nbeta; b++)
        get_value_glb(fin, csw[b]);
      break;
    case LAMBDAQCD_TK:
      get_value_glb(fin, LambdaQCD);
      break;
    case P2MIN_TK:
      get_value_glb(fin, p2min);
      break;
    case P2MAX_TK:
      get_value_glb(fin, p2max);
      break;
    case THRESH_TK:
      get_value_glb(fin, thresh);
      break;
    case OUT_HADR_TK:
      get_value_glb(fin, out_hadr);
      break;
    case ONLY_BASIC_TK:
      get_value_glb(fin, only_basic);
      break;
    case COMPUTE_MPCAC_TK:
      get_value_glb(fin, compute_mpcac);
      break;
    case ANALYSIS_TK:
      get_value_glb(fin, analysis);
      break;
    case CLOVER_TK:
      get_value_glb(fin, clover);
      break;
    case P2REF_TK:
      get_value_glb(fin, p2ref);
      break;
    case LOAD_AVE_TK:
      get_value_glb(fin, load_ave);
      break;
    case LOAD_CHIR_TK:
      get_value_glb(fin, load_chir);
      break;
    case SUFFIX_TK:
      get_value_glb(fin, an_suffix);
      break;
    case QCDONTHERIGHT_TK:
      get_value_glb(fin, QCD_on_the_right);
      break;
    case SUB_BOOSTED_TK:
      get_value_glb(fin, sub_boosted);
      break;
    case SUB_PTILDE_TK:
      get_value_glb(fin, sub_ptilde);
      break;
    case SUBTRACTION_TK:
      get_value_glb(fin, subtraction);
      break;
    case CHIR_ANSATZ_VAL_TK:
      get_value_glb(fin, chir_ansatz_val);
      break;
    case CHIR_ANSATZ_SEA_TK:
      get_value_glb(fin, chir_ansatz_sea);
      break;
    case STEPFUNC_MIN_TK:
      get_value_glb(fin, stepfunc_min);
      break;
    case STEPFUNC_MAX_TK:
      get_value_glb(fin, stepfunc_max);
      break;
    case P2MAX_M3_M4_TK:
      p2max_M3_M4.resize(nbeta);
      for (int b = 0; b < nbeta; b++)
        get_value_glb(fin, p2max_M3_M4[b]);
      break;
    case P2MIN_M3_M4_TK:
      p2min_M3_M4.resize(nbeta);
      for (int b = 0; b < nbeta; b++)
        get_value_glb(fin, p2min_M3_M4[b]);
      break;
    case LOAD_LABEL_TK:
      get_value_glb(fin, load_label);
      break;

    case FEOF_GLB_TK:
      break;
    }
  }

  // check initialization
  check_int_par(nconfs, nconfs_tag);
  check_int_par(njacks, njacks_tag);
  check_str_par(path_folder, path_folder_tag);
  check_str_par(scheme, scheme_tag);
  check_str_par(BC, BC_tag);
  check_str_par(action, act_tag);
  check_int_par(Nc, Nc_tag);
  check_int_par(Nf, Nf_tag);
  check_int_par(nr, nr_tag);
  check_int_par(ntypes, ntypes_tag);
  check_int_par(nhits, nhits_tag);
  check_int_par(UseSigma1, UseSigma1_tag);
  check_int_par(UseEffMass, UseEffMass_tag);
  for (auto &b : beta)
    check_double_par(b, beta_tag);
  for (auto &bl : beta_label)
    check_str_par(bl.c_str(), beta_label_tag);
  for (auto &vl : volume_label)
    check_str_par(vl.c_str(), volume_label_tag);
  for (auto &a : ainv)
    check_double_par(a, ainv_tag);
  for (auto &c : csw)
    check_double_par(c, csw_tag);
  for (auto &ms : nm_Sea)
    check_int_par(ms, nm_Sea_tag);
  for (auto &im : SeaMasses_label)
    for (auto &jm : im)
      check_int_par(jm, SeaMasses_label_tag);
  for (auto &im : SeaMasses)
    for (auto &jm : im)
      check_double_par(jm, SeaMasses_tag);
  check_double_par(LambdaQCD, LambdaQCD_tag);
  check_double_par(p2min, p2min_tag);
  check_double_par(p2max, p2max_tag);
  check_double_par(thresh, thresh_tag);
  check_str_par(out_hadr, out_hadr_tag);
  check_int_par(only_basic, only_basic_tag);
  check_int_par(compute_mpcac, compute_mpcac_tag);
  check_str_par(analysis, analysis_tag);
  check_str_par(clover, clover_tag);
  check_double_par(p2ref, p2ref_tag);
  check_int_par(load_ave, load_ave_tag);
  check_int_par(load_chir, load_chir_tag);
  check_str_par(an_suffix, an_suffix_tag);
  check_int_par(QCD_on_the_right, QCD_on_the_right_tag);
  check_int_par(sub_boosted, sub_boosted_tag);
  check_int_par(sub_ptilde, sub_ptilde_tag);
  check_int_par(subtraction, subtraction_tag);
  check_str_par(chir_ansatz_val, chir_ansatz_val_tag);
  check_str_par(chir_ansatz_sea, chir_ansatz_sea_tag);
  check_double_par(stepfunc_min, stepfunc_min_tag);
  check_double_par(stepfunc_max, stepfunc_max_tag);
  for (auto &p : p2max_M3_M4)
    check_double_par(p, p2max_M3_M4_tag);
  for (auto &p : p2min_M3_M4)
    check_double_par(p, p2min_M3_M4_tag);
  check_str_par(load_label, load_label_tag);

  fclose(fin);

  clover_analysis = false;
  inte_analysis = false;

  if (strcmp(clover.c_str(), "no") == 0)
  {
    if (strcmp(analysis.c_str(), "inte") == 0)
    {
      path_analysis = {"Nf4" + an_suffix};

      inte_analysis = true;
    }
    else
    {
      cout << "Choose the analysis: 'inte'." << endl;
      exit(0);
    }
  }
  else if (strcmp(clover.c_str(), "yes") == 0)
  {

    clover_analysis = true;

    if (strcmp(analysis.c_str(), "inte") == 0)
    {
      path_analysis = {"Nf4" + an_suffix};

      inte_analysis = true;
    }
    else
    {
      cout << "Only 'inte' analysis implemented." << endl;
      exit(0);
    }
  }
  else
  {
    cout << "Specify if the analysis is Clover: 'yes' or 'no'." << endl;
    exit(0);
  }

  load = (load_ave or load_chir);

  if (load and only_basic)
  {
    cout << "Cannot load saved quantities in the only_basic mode." << endl;
    exit(0);
  }

  if (ntypes != 1)
  {
    cout << "Only ntypes=1 implemented." << endl;
    exit(0);
  }

  if (strcmp(chir_ansatz_val.c_str(), "linear") != 0 and
      strcmp(chir_ansatz_val.c_str(), "constant") != 0 and
      strcmp(chir_ansatz_val.c_str(), "quadratic") != 0)
  {
    cout << "Choose the valence chiral fit Ansatz among: \"linear/constant/quadratic\"." << endl;
    exit(0);
  }
  if (strcmp(chir_ansatz_sea.c_str(), "linear") != 0 and
      strcmp(chir_ansatz_sea.c_str(), "constant") != 0 and
      strcmp(chir_ansatz_sea.c_str(), "quadratic") != 0)
  {
    cout << "Choose the sea chiral fit Ansatz among: \"linear/constant/quadratic\"." << endl;
    exit(0);
  }

  // this is the path to the directory which contains 'print', 'plots', ecc.
  string full_path = path_folder + path_analysis[0] + "/";

  // evaluate max number of sea masses for the ensembles
  nm_Sea_max = *max_element(nm_Sea.begin(), nm_Sea.end());

  // print input parameters
  const char *globalConf = "Global configuration";
  printBox(globalConf);

  printf(" %s = %s", analysis_tag, analysis.c_str()); // free, inte, ratio
  if (clover_analysis)
  {
    printf(" with Clover\n\n");
  }
  else
  {
    printf("\n\n");
  }

  printf(" %s = %s\n", scheme_tag, scheme.c_str());
  printf("    with BC: %s \n\n", BC.c_str());

  printf(" %s = %.2lf\n", thresh_tag, thresh);
  printf(" Continuum limit range: p2 = [%.1lf,%.1lf]\n\n", p2min, p2max);

  printf(" Evolution at the scale: p2ref = %.1lf\n", p2ref);
  printf(" Step scaling test: (%.2lf,%.2lf) GeV^2 (Z[p2_max]/Z[p2_min])\n", stepfunc_min,
         stepfunc_max);
  printf(" p2-range for M3 and M4 methods: ");
  for (int b = 0; b < nbeta; b++)
    printf("[%.1lf-%.1lf] GeV2  ", p2min_M3_M4[b], p2max_M3_M4[b]);
  printf("\n\n");

  printf(" %s = %s  --  %s = %d  -- %s = %d -- %s = %.3lf \n", act_tag, action.c_str(), Nf_tag, Nf,
         Nc_tag, Nc, LambdaQCD_tag, LambdaQCD);
  printf(" %s = %d  (%d njacks) \n", nconfs_tag, nconfs, njacks);
  printf(" %s = %s \n\n", path_folder_tag, full_path.c_str());

  printf(" %s = %d\n", nr_tag, nr);
  printf(" %s = %d\n", ntypes_tag, ntypes);
  printf(" %s = %d\n\n", nhits_tag, nhits);

  printf(" Working with %d beta: \n", nbeta);
  for (int b = 0; b < nbeta; b++)
  {
    printf("    beta = %.3lf : ainv = %.2lf", beta[b], ainv[b]);
    if (clover_analysis)
    {
      printf("   -  csw = %.4lf\n", csw[b]);
    }
    else
    {
      printf("\n");
    }
    printf("                  Ensembles: ");
    if (!clover_analysis)
      for (int m = 0; m < nm_Sea[b]; m++)
        printf("%s%d ", beta_label[b].c_str(), SeaMasses_label[b][m]);
    else
      for (int m = 0; m < nm_Sea[b]; m++)
        printf("%s.d.%d.%s ", beta_label[b].c_str(), SeaMasses_label[b][m],
               volume_label[b].c_str());
    printf("\n");
  }

  printf(" Using Zq according to ");
  if (UseSigma1)
    printf("RI'-MOM variant (eq.33 of 1004.1115). \n");
  else
    printf("RI'-MOM prescription. \n");

  printf(" Using Z^{QCD} factorized on the ");
  if (QCD_on_the_right)
    printf("RIGHT. \n");
  else
    printf("LEFT. \n");

  if (subtraction == 0)
    printf(" No perturbative subtraction.\n");
  else if (subtraction > 0 and subtraction < 3)
  {
    printf(" Perturbative subtraction done at order ");
    if (subtraction == 1)
      printf("O(g2a2)");
    else if (subtraction == 2)
      printf("O(g2ainf)");
    printf(" with ");
    if (sub_boosted)
      printf("g^2[boosted] = g0^2/<Plaq>  (g0^2 = 6/beta)");
    else
      printf("g0^2 = 6/beta");
    printf(" and using ");
    if (sub_ptilde)
      printf("p2tilde.\n");
    else
      printf("p2.\n");
  }
  else
  {
    printf(" Input parameter not valid. Choose between:\n");
    printf("   0: no subtraction\n");
    printf("   1: O(g2a2)\n");
    printf("   2: O(g2ainf)\n\n");
    exit(0);
  }
  printf(" Using [%s] Ansatz for valence chiral extrapolation.\n", chir_ansatz_val.c_str());
  printf(" Using [%s] Ansatz for sea chiral extrapolation.\n", chir_ansatz_sea.c_str());

  printf("\n");

  if (only_basic)
    printf(" Computing only basic quantities. \n");

  if (load_ave)
    printf(" Loading averaged quantities: '%s' \n", load_label.c_str());
  else
    printf(" Saving averaged quantities as '%s' \n", load_label.c_str());

  // define global variables from input
  clust_size = nconfs / njacks;

  nbil = 5;

  // slightly increment thresh to include border
  thresh *= 1 + 1e-10;

  printf("\n\n");
}

void read_input(const string &path_to_ens, const string &name)
{
  string path_to_input = path_to_ens + "input.txt";

  FILE *fin = fopen(path_to_input.c_str(), "r");
  if (not fin)
  {
    fprintf(stderr, "Failed to open \"%s\"\n", path_to_input.c_str());
    exit(FAILED_OPEN);
  }

  mom_path = DEFAULT_STR_VAL;

  int L = DEFAULT_INT_VAL;
  int T = DEFAULT_INT_VAL;

  conf_init = DEFAULT_INT_VAL;
  conf_step = DEFAULT_INT_VAL;
  nm = DEFAULT_INT_VAL;
  delta_tmin = DEFAULT_INT_VAL;
  delta_tmax = DEFAULT_INT_VAL;
  kappa = DEFAULT_DOUBLE_VAL;
  mu_sea = DEFAULT_DOUBLE_VAL;
  plaquette = DEFAULT_DOUBLE_VAL;
  //    for(auto &m : mass_val) m=DEFAULT_DOUBLE_VAL;

  while (not feof(fin))
  {
    TK_t tk = get_TK(fin);
    switch (tk)
    {
    case VALUE_TK:
      fprintf(stderr, "Invalid token %s found\n", tok);
      exit(1);
      break;
    case MOM_LIST_TK:
      get_value(fin, mom_path);
      break;
    case L_TK:
      get_value(fin, L);
      break;
    case T_TK:
      get_value(fin, T);
      break;
    case CONF_INIT_TK:
      get_value(fin, conf_init);
      break;
    case CONF_STEP_TK:
      get_value(fin, conf_step);
      break;
    case PLAQ_TK:
      get_value(fin, plaquette);
      break;
    case KAPPA_TK:
      get_value(fin, kappa);
      break;
    case MU_SEA_TK:
      get_value(fin, mu_sea);
      break;
    case NM_VAL_TK:
      get_value(fin, nm);
      break;
    case VALMASSES_TK:
      mass_val.resize(nm);
      for (int i = 0; i < nm; i++)
        get_value(fin, mass_val[i]);
      break;
    case DELTA_TMIN_TK:
      get_value(fin, delta_tmin);
      break;
    case DELTA_TMAX_TK:
      get_value(fin, delta_tmax);
      break;

    case FEOF_TK:
      break;
    }
  }

  check_str_par(mom_path.c_str(), mom_list_tag);
  check_int_par(L, L_tag);
  check_int_par(T, T_tag);
  check_double_par(plaquette, plaquette_tag);
  check_double_par(kappa, kappa_tag);
  check_int_par(conf_init, conf_init_tag);
  check_int_par(conf_step, conf_step_tag);
  check_int_par(nm, nm_tag);
  for (int i = 0; i < nm; i++)
    check_double_par(mass_val[i], mass_val_tag);
  check_int_par(delta_tmin, delta_tmin_tag);
  check_int_par(delta_tmax, delta_tmax_tag);

  fclose(fin);

  if (plaquette == 0.0)
    plaquette = read_plaquette(path_to_ens);

  sizeV = {T, L, L, L};

  // print ensemble details
  string ensemble = "Ensemble " + name;
  const char *ensembleString = ensemble.c_str();
  printBox(ensembleString);

  printf(" %s = \"%s\"\n", mom_list_tag, mom_path.c_str());
  printf(" Dimensions = %dc%d\n", L, T);
  printf(" %s = %lf -- %s = %.6lf -- %s = %.4lf \n", plaquette_tag, plaquette, kappa_tag, kappa,
         mu_sea_tag, mu_sea);

  printf(" %s = %d  [from %d to %d]\n", nconfs_tag, nconfs, conf_init,
         conf_init + (nconfs - 1) * conf_step);

  printf(" Fit range for deltam_cr: [%d,%d]\n", delta_tmin, delta_tmax);

  printf(" %s = ", mass_val_tag);
  for (int i = 0; i < nm; i++)
    printf("%.4lf  ", mass_val[i]);
  printf("\n\n");

  nmr = nm * nr;
  combo = nm * nr * ntypes * nhits * nconfs;
  neq = fact(nm + nr - 1) / fact(nr) / fact(nm - 1);
  neq2 = nm;
}
