// Generated by dust (version 0.12.9) - do not edit
#include <cpp11.hpp>

[[cpp11::register]]
cpp11::sexp dust_InfectionSplitModelOD_capabilities();

[[cpp11::register]]
cpp11::sexp dust_InfectionSplitModelOD_gpu_info();
[[cpp11::register]]
SEXP dust_cpu_InfectionSplitModelOD_alloc(cpp11::list r_pars, bool pars_multi, cpp11::sexp r_time,
                         cpp11::sexp r_n_particles, int n_threads,
                         cpp11::sexp r_seed, bool deterministic,
                         cpp11::sexp gpu_config, cpp11::sexp ode_control);

[[cpp11::register]]
SEXP dust_cpu_InfectionSplitModelOD_run(SEXP ptr, cpp11::sexp r_time_end);

[[cpp11::register]]
SEXP dust_cpu_InfectionSplitModelOD_simulate(SEXP ptr, cpp11::sexp time_end);

[[cpp11::register]]
SEXP dust_cpu_InfectionSplitModelOD_set_index(SEXP ptr, cpp11::sexp r_index);

[[cpp11::register]]
SEXP dust_cpu_InfectionSplitModelOD_update_state(SEXP ptr, SEXP r_pars, SEXP r_state,
                                           SEXP r_time, SEXP r_set_initial_state,
                                           SEXP index, SEXP reset_step_size);

[[cpp11::register]]
SEXP dust_cpu_InfectionSplitModelOD_state(SEXP ptr, SEXP r_index);

[[cpp11::register]]
SEXP dust_cpu_InfectionSplitModelOD_time(SEXP ptr);

[[cpp11::register]]
void dust_cpu_InfectionSplitModelOD_reorder(SEXP ptr, cpp11::sexp r_index);

[[cpp11::register]]
SEXP dust_cpu_InfectionSplitModelOD_resample(SEXP ptr, cpp11::doubles r_weights);

[[cpp11::register]]
SEXP dust_cpu_InfectionSplitModelOD_rng_state(SEXP ptr, bool first_only, bool last_only);

[[cpp11::register]]
SEXP dust_cpu_InfectionSplitModelOD_set_rng_state(SEXP ptr, cpp11::raws rng_state);

[[cpp11::register]]
SEXP dust_cpu_InfectionSplitModelOD_set_data(SEXP ptr, cpp11::list data, bool shared);

[[cpp11::register]]
SEXP dust_cpu_InfectionSplitModelOD_compare_data(SEXP ptr);

[[cpp11::register]]
SEXP dust_cpu_InfectionSplitModelOD_filter(SEXP ptr, SEXP time_end,
                                     bool save_trajectories,
                                     cpp11::sexp time_snapshot,
                                     cpp11::sexp min_log_likelihood);

[[cpp11::register]]
void dust_cpu_InfectionSplitModelOD_set_n_threads(SEXP ptr, int n_threads);

[[cpp11::register]]
int dust_cpu_InfectionSplitModelOD_n_state(SEXP ptr);

[[cpp11::register]]
void dust_cpu_InfectionSplitModelOD_set_stochastic_schedule(SEXP ptr, SEXP time);

[[cpp11::register]]
void dust_cpu_InfectionSplitModelOD_ode_statistics(SEXP ptr);
#include <dust/r/dust.hpp>

// Generated by odin.dust (version 0.2.29) - do not edit
template <typename real_type, typename container>
__host__ __device__ real_type odin_sum1(const container x, size_t from, size_t to);
template <typename real_type, typename container>
__host__ __device__ real_type odin_sum2(const container x, int from_i, int to_i, int from_j, int to_j, int dim_x_1);
template <typename real_type, typename T, typename U>
__host__ __device__ real_type fmodr(T x, U y) {
  real_type tmp = std::fmod(static_cast<real_type>(x),
                            static_cast<real_type>(y));
  if (tmp * y < 0) {
    tmp += y;
  }
  return tmp;
}

// These exist to support the model on the gpu, as in C++14 std::min
// and std::max are constexpr and error without --expt-relaxed-constexpr
template <typename T>
__host__ __device__ T odin_min(T x, T y) {
  return x < y ? x : y;
}

template <typename T>
__host__ __device__ T odin_max(T x, T y) {
  return x > y ? x : y;
}
// [[dust::class(InfectionSplitModelOD)]]
// [[dust::param(Cas0, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(dP1_all, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(dP2_all, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(dt, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(Exp0, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(FOI_spillover, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(Inf0, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(N_age, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(n_years, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(R0, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(Rec0, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(Sus0, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(Vac0, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(vacc_rate_annual, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(vaccine_efficacy, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(year0, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
class InfectionSplitModelOD {
public:
  using real_type = double;
  using rng_state_type = dust::random::generator<real_type>;
  using data_type = dust::no_data;
  struct shared_type {
    real_type beta;
    std::vector<real_type> Cas0;
    int dim_C_sylvatic;
    int dim_C_urban;
    int dim_Cas0;
    int dim_dP1;
    int dim_dP1_all;
    int dim_dP1_all_1;
    int dim_dP1_all_2;
    int dim_dP2;
    int dim_dP2_all;
    int dim_dP2_all_1;
    int dim_dP2_all_2;
    int dim_E_new_sylvatic;
    int dim_E_new_urban;
    int dim_E_sylvatic;
    int dim_E_urban;
    int dim_Exp0;
    int dim_I_new_sylvatic;
    int dim_I_new_urban;
    int dim_I_sylvatic;
    int dim_I_urban;
    int dim_Inf0;
    int dim_inv_P;
    int dim_inv_P_nV;
    int dim_P;
    int dim_P_nV;
    int dim_R;
    int dim_R_new_sylvatic;
    int dim_R_new_urban;
    int dim_Rec0;
    int dim_S;
    int dim_Sus0;
    int dim_V;
    int dim_Vac0;
    int dim_vacc_rate;
    int dim_vacc_rate_annual;
    int dim_vacc_rate_annual_1;
    int dim_vacc_rate_annual_2;
    std::vector<real_type> dP1_all;
    std::vector<real_type> dP2_all;
    real_type dt;
    std::vector<real_type> Exp0;
    real_type FOI_max;
    real_type FOI_spillover;
    real_type FOI_sylvatic;
    std::vector<real_type> Inf0;
    std::vector<real_type> initial_C_sylvatic;
    std::vector<real_type> initial_C_urban;
    real_type initial_day;
    std::vector<real_type> initial_E_sylvatic;
    std::vector<real_type> initial_E_urban;
    real_type initial_FOI_sylv;
    real_type initial_FOI_urb;
    std::vector<real_type> initial_I_sylvatic;
    std::vector<real_type> initial_I_urban;
    std::vector<real_type> initial_R;
    std::vector<real_type> initial_S;
    real_type initial_time;
    std::vector<real_type> initial_V;
    real_type initial_year;
    int N_age;
    int n_years;
    int offset_variable_C_sylvatic;
    int offset_variable_C_urban;
    int offset_variable_E_sylvatic;
    int offset_variable_E_urban;
    int offset_variable_I_sylvatic;
    int offset_variable_I_urban;
    int offset_variable_R;
    int offset_variable_V;
    real_type Pmin;
    real_type R0;
    real_type rate1;
    real_type rate2;
    std::vector<real_type> Rec0;
    std::vector<real_type> Sus0;
    real_type t_incubation;
    real_type t_infectious;
    real_type t_latent;
    std::vector<real_type> Vac0;
    std::vector<real_type> vacc_rate_annual;
    real_type vaccine_efficacy;
    real_type year0;
  };
  struct internal_type {
    std::vector<real_type> dP1;
    std::vector<real_type> dP2;
    std::vector<real_type> E_new_sylvatic;
    std::vector<real_type> E_new_urban;
    std::vector<real_type> I_new_sylvatic;
    std::vector<real_type> I_new_urban;
    std::vector<real_type> inv_P;
    std::vector<real_type> inv_P_nV;
    std::vector<real_type> P;
    std::vector<real_type> P_nV;
    std::vector<real_type> R_new_sylvatic;
    std::vector<real_type> R_new_urban;
    std::vector<real_type> vacc_rate;
  };
  InfectionSplitModelOD(const dust::pars_type<InfectionSplitModelOD>& pars) :
    shared(pars.shared), internal(pars.internal) {
  }
  size_t size() {
    return shared->dim_C_sylvatic + shared->dim_C_urban + shared->dim_E_sylvatic + shared->dim_E_urban + shared->dim_I_sylvatic + shared->dim_I_urban + shared->dim_R + shared->dim_S + shared->dim_V + 5;
  }
  std::vector<real_type> initial(size_t step) {
    std::vector<real_type> state(shared->dim_C_sylvatic + shared->dim_C_urban + shared->dim_E_sylvatic + shared->dim_E_urban + shared->dim_I_sylvatic + shared->dim_I_urban + shared->dim_R + shared->dim_S + shared->dim_V + 5);
    state[0] = shared->initial_time;
    state[1] = shared->initial_day;
    state[2] = shared->initial_year;
    state[3] = shared->initial_FOI_sylv;
    state[4] = shared->initial_FOI_urb;
    std::copy(shared->initial_S.begin(), shared->initial_S.end(), state.begin() + 5);
    std::copy(shared->initial_E_sylvatic.begin(), shared->initial_E_sylvatic.end(), state.begin() + shared->offset_variable_E_sylvatic);
    std::copy(shared->initial_E_urban.begin(), shared->initial_E_urban.end(), state.begin() + shared->offset_variable_E_urban);
    std::copy(shared->initial_I_sylvatic.begin(), shared->initial_I_sylvatic.end(), state.begin() + shared->offset_variable_I_sylvatic);
    std::copy(shared->initial_I_urban.begin(), shared->initial_I_urban.end(), state.begin() + shared->offset_variable_I_urban);
    std::copy(shared->initial_R.begin(), shared->initial_R.end(), state.begin() + shared->offset_variable_R);
    std::copy(shared->initial_V.begin(), shared->initial_V.end(), state.begin() + shared->offset_variable_V);
    std::copy(shared->initial_C_sylvatic.begin(), shared->initial_C_sylvatic.end(), state.begin() + shared->offset_variable_C_sylvatic);
    std::copy(shared->initial_C_urban.begin(), shared->initial_C_urban.end(), state.begin() + shared->offset_variable_C_urban);
    return state;
  }
  void update(size_t step, const real_type * state, rng_state_type& rng_state, real_type * state_next) {
    const real_type day = state[1];
    const real_type * S = state + 5;
    const real_type * E_sylvatic = state + shared->offset_variable_E_sylvatic;
    const real_type * E_urban = state + shared->offset_variable_E_urban;
    const real_type * I_sylvatic = state + shared->offset_variable_I_sylvatic;
    const real_type * I_urban = state + shared->offset_variable_I_urban;
    const real_type * R = state + shared->offset_variable_R;
    const real_type * V = state + shared->offset_variable_V;
    state_next[1] = day + shared->dt;
    state_next[0] = (step + 1) * shared->dt;
    real_type year_i = dust::math::floor((step * shared->dt) / (real_type) 365) + 1;
    state_next[3] = shared->FOI_sylvatic;
    state_next[2] = year_i + shared->year0 - 1;
    for (int i = 1; i <= shared->N_age; ++i) {
      internal.E_new_sylvatic[i - 1] = dust::random::binomial<real_type>(rng_state, static_cast<int>(S[i - 1]), shared->FOI_sylvatic);
    }
    for (int i = 1; i <= shared->N_age; ++i) {
      internal.I_new_sylvatic[i - 1] = E_sylvatic[i - 1] * shared->rate1;
    }
    for (int i = 1; i <= shared->N_age; ++i) {
      internal.I_new_urban[i - 1] = E_urban[i - 1] * shared->rate1;
    }
    for (int i = 1; i <= shared->N_age; ++i) {
      internal.P[i - 1] = S[i - 1] + E_sylvatic[i - 1] + E_urban[i - 1] + I_sylvatic[i - 1] + I_urban[i - 1] + R[i - 1] + V[i - 1];
    }
    for (int i = 1; i <= shared->N_age; ++i) {
      internal.P_nV[i - 1] = S[i - 1] + R[i - 1];
    }
    for (int i = 1; i <= shared->N_age; ++i) {
      internal.R_new_sylvatic[i - 1] = I_sylvatic[i - 1] * shared->rate2;
    }
    for (int i = 1; i <= shared->N_age; ++i) {
      internal.R_new_urban[i - 1] = I_urban[i - 1] * shared->rate2;
    }
    for (int i = 1; i <= shared->N_age; ++i) {
      internal.dP1[i - 1] = shared->dP1_all[shared->dim_dP1_all_1 * (static_cast<int>(year_i) - 1) + i - 1] * shared->dt;
    }
    for (int i = 1; i <= shared->N_age; ++i) {
      internal.dP2[i - 1] = shared->dP2_all[shared->dim_dP2_all_1 * (static_cast<int>(year_i) - 1) + i - 1] * shared->dt;
    }
    for (int i = 1; i <= shared->N_age; ++i) {
      internal.inv_P[i - 1] = 1 / (real_type) internal.P[i - 1];
    }
    for (int i = 1; i <= shared->N_age; ++i) {
      internal.inv_P_nV[i - 1] = 1 / (real_type) internal.P_nV[i - 1];
    }
    real_type P_tot = odin_sum1<real_type>(internal.P.data(), 0, shared->dim_P);
    for (int i = 1; i <= shared->N_age; ++i) {
      state_next[shared->offset_variable_C_sylvatic + i - 1] = internal.I_new_sylvatic[i - 1];
    }
    for (int i = 1; i <= shared->N_age; ++i) {
      state_next[shared->offset_variable_C_urban + i - 1] = internal.I_new_urban[i - 1];
    }
    for (int i = 1; i <= shared->N_age; ++i) {
      state_next[shared->offset_variable_E_sylvatic + i - 1] = dust::math::max(shared->Pmin, E_sylvatic[i - 1] + internal.E_new_sylvatic[i - 1] - internal.I_new_sylvatic[i - 1]);
    }
    for (int i = 1; i <= shared->N_age; ++i) {
      state_next[shared->offset_variable_I_sylvatic + i - 1] = dust::math::max(shared->Pmin, I_sylvatic[i - 1] + internal.I_new_sylvatic[i - 1] - internal.R_new_sylvatic[i - 1]);
    }
    for (int i = 1; i <= shared->N_age; ++i) {
      state_next[shared->offset_variable_I_urban + i - 1] = dust::math::max(shared->Pmin, I_urban[i - 1] + internal.I_new_urban[i - 1] - internal.R_new_urban[i - 1]);
    }
    for (int i = 1; i <= shared->N_age; ++i) {
      internal.vacc_rate[i - 1] = shared->vacc_rate_annual[shared->dim_vacc_rate_annual_1 * (static_cast<int>(year_i) - 1) + i - 1] * shared->vaccine_efficacy * shared->dt * internal.P[i - 1];
    }
    real_type FOI_urban = dust::math::min(shared->FOI_max, shared->beta * ((odin_sum1<real_type>(I_sylvatic, 0, shared->dim_I_sylvatic) + odin_sum1<real_type>(I_urban, 0, shared->dim_I_urban)) / (real_type) P_tot));
    {
       int i = 1;
       state_next[shared->offset_variable_R + i - 1] = dust::math::max(shared->Pmin, R[0] + internal.R_new_sylvatic[0] + internal.R_new_urban[0] - internal.vacc_rate[0] * R[0] * internal.inv_P_nV[0] - (internal.dP2[0] * R[0] * internal.inv_P[0]));
    }
    for (int i = 2; i <= shared->N_age; ++i) {
      state_next[shared->offset_variable_R + i - 1] = dust::math::max(shared->Pmin, R[i - 1] + internal.R_new_sylvatic[i - 1] + internal.R_new_urban[i - 1] - internal.vacc_rate[i - 1] * R[i - 1] * internal.inv_P_nV[i - 1] + (internal.dP1[i - 1] * R[i - 1 - 1] * internal.inv_P[i - 1 - 1]) - (internal.dP2[i - 1] * R[i - 1] * internal.inv_P[i - 1]));
    }
    {
       int i = 1;
       state_next[shared->offset_variable_V + i - 1] = dust::math::max(shared->Pmin, V[0] + internal.vacc_rate[0] - (internal.dP2[0] * V[0] * internal.inv_P[0]));
    }
    for (int i = 2; i <= shared->N_age; ++i) {
      state_next[shared->offset_variable_V + i - 1] = dust::math::max(shared->Pmin, V[i - 1] + internal.vacc_rate[i - 1] + (internal.dP1[i - 1] * V[i - 1 - 1] * internal.inv_P[i - 1 - 1]) - (internal.dP2[i - 1] * V[i - 1] * internal.inv_P[i - 1]));
    }
    for (int i = 1; i <= shared->N_age; ++i) {
      internal.E_new_urban[i - 1] = dust::random::binomial<real_type>(rng_state, static_cast<int>(S[i - 1]), FOI_urban);
    }
    state_next[4] = FOI_urban;
    for (int i = 1; i <= shared->N_age; ++i) {
      state_next[shared->offset_variable_E_urban + i - 1] = dust::math::max(shared->Pmin, E_urban[i - 1] + internal.E_new_urban[i - 1] - internal.I_new_urban[i - 1]);
    }
    {
       int i = 1;
       state_next[5 + i - 1] = dust::math::max(shared->Pmin, S[0] - internal.E_new_sylvatic[0] - internal.E_new_urban[0] - internal.vacc_rate[0] * S[0] * internal.inv_P_nV[0] + internal.dP1[0] - (internal.dP2[0] * S[0] * internal.inv_P[0]));
    }
    for (int i = 2; i <= shared->N_age; ++i) {
      state_next[5 + i - 1] = dust::math::max(shared->Pmin, S[i - 1] - internal.E_new_sylvatic[i - 1] - internal.E_new_urban[i - 1] - internal.vacc_rate[i - 1] * S[i - 1] * internal.inv_P_nV[i - 1] + (internal.dP1[i - 1] * S[i - 1 - 1] * internal.inv_P[i - 1 - 1]) - (internal.dP2[i - 1] * S[i - 1] * internal.inv_P[i - 1]));
    }
  }
private:
  std::shared_ptr<const shared_type> shared;
  internal_type internal;
};
template <typename real_type, typename container>
__host__ __device__ real_type odin_sum2(const container x, int from_i, int to_i, int from_j, int to_j, int dim_x_1) {
  real_type tot = 0.0;
  for (int j = from_j; j < to_j; ++j) {
    int jj = j * dim_x_1;
    for (int i = from_i; i < to_i; ++i) {
      tot += x[i + jj];
    }
  }
  return tot;
}
#include <array>
#include <cpp11/R.hpp>
#include <cpp11/sexp.hpp>
#include <cpp11/doubles.hpp>
#include <cpp11/integers.hpp>
#include <cpp11/list.hpp>
#include <cpp11/strings.hpp>
#include <memory>
#include <vector>

template <typename T>
inline bool is_na(T x);

template <>
inline bool is_na(int x) {
  return x == NA_INTEGER;
}

template <>
inline bool is_na(double x) {
  return ISNA(x);
}

inline size_t object_length(cpp11::sexp x) {
  return ::Rf_xlength(x);
}

template <typename T>
void user_check_value(T value, const char *name, T min, T max) {
  if (is_na(value)) {
    cpp11::stop("'%s' must not be NA", name);
  }
  if (!is_na(min) && value < min) {
    cpp11::stop("Expected '%s' to be at least %g", name, (double) min);
  }
  if (!is_na(max) && value > max) {
    cpp11::stop("Expected '%s' to be at most %g", name, (double) max);
  }
}

template <typename T>
void user_check_array_value(const std::vector<T>& value, const char *name,
                            T min, T max) {
  for (auto& x : value) {
    user_check_value(x, name, min, max);
  }
}

inline size_t user_get_array_rank(cpp11::sexp x) {
  if (!::Rf_isArray(x)) {
    return 1;
  } else {
    cpp11::integers dim = cpp11::as_cpp<cpp11::integers>(x.attr("dim"));
    return dim.size();
  }
}

template <size_t N>
void user_check_array_rank(cpp11::sexp x, const char *name) {
  size_t rank = user_get_array_rank(x);
  if (rank != N) {
    if (N == 1) {
      cpp11::stop("Expected a vector for '%s'", name);
    } else if (N == 2) {
      cpp11::stop("Expected a matrix for '%s'", name);
    } else {
      cpp11::stop("Expected an array of rank %d for '%s'", N, name);
    }
  }
}

template <size_t N>
void user_check_array_dim(cpp11::sexp x, const char *name,
                          const std::array<int, N>& dim_expected) {
  cpp11::integers dim = cpp11::as_cpp<cpp11::integers>(x.attr("dim"));
  for (size_t i = 0; i < N; ++i) {
    if (dim[(int)i] != dim_expected[i]) {
      Rf_error("Incorrect size of dimension %d of '%s' (expected %d)",
               i + 1, name, dim_expected[i]);
    }
  }
}

template <>
inline void user_check_array_dim<1>(cpp11::sexp x, const char *name,
                                    const std::array<int, 1>& dim_expected) {
  if ((int)object_length(x) != dim_expected[0]) {
    cpp11::stop("Expected length %d value for '%s'", dim_expected[0], name);
  }
}

template <size_t N>
void user_set_array_dim(cpp11::sexp x, const char *name,
                        std::array<int, N>& dim) {
  cpp11::integers dim_given = cpp11::as_cpp<cpp11::integers>(x.attr("dim"));
  std::copy(dim_given.begin(), dim_given.end(), dim.begin());
}

template <>
inline void user_set_array_dim<1>(cpp11::sexp x, const char *name,
                                  std::array<int, 1>& dim) {
  dim[0] = object_length(x);
}

template <typename T>
T user_get_scalar(cpp11::list user, const char *name,
                  const T previous, T min, T max) {
  T ret = previous;
  cpp11::sexp x = user[name];
  if (x != R_NilValue) {
    if (object_length(x) != 1) {
      cpp11::stop("Expected a scalar numeric for '%s'", name);
    }
    // TODO: when we're getting out an integer this is a bit too relaxed
    if (TYPEOF(x) == REALSXP) {
      ret = cpp11::as_cpp<T>(x);
    } else if (TYPEOF(x) == INTSXP) {
      ret = cpp11::as_cpp<T>(x);
    } else {
      cpp11::stop("Expected a numeric value for %s", name);
    }
  }

  if (is_na(ret)) {
    cpp11::stop("Expected a value for '%s'", name);
  }
  user_check_value<T>(ret, name, min, max);
  return ret;
}

template <>
inline float user_get_scalar<float>(cpp11::list user, const char *name,
                                    const float previous, float min, float max) {
  double value = user_get_scalar<double>(user, name, previous, min, max);
  return static_cast<float>(value);
}

template <typename T>
std::vector<T> user_get_array_value(cpp11::sexp x, const char * name,
                                    T min, T max) {
  std::vector<T> ret = cpp11::as_cpp<std::vector<T>>(x);
  user_check_array_value<T>(ret, name, min, max);
  return ret;
}

template <typename T, size_t N>
std::vector<T> user_get_array_fixed(cpp11::list user, const char *name,
                                    const std::vector<T> previous,
                                    const std::array<int, N>& dim,
                                    T min, T max) {
  cpp11::sexp x = user[name];
  if (x == R_NilValue) {
    if (previous.size() == 0) {
      cpp11::stop("Expected a value for '%s'", name);
    }
    return previous;
  }

  user_check_array_rank<N>(x, name);
  user_check_array_dim<N>(x, name, dim);

  return user_get_array_value<T>(x, name, min, max);
}

template <typename T, size_t N>
std::vector<T> user_get_array_variable(cpp11::list user, const char *name,
                                       std::vector<T> previous,
                                       std::array<int, N>& dim,
                                       T min, T max) {
  cpp11::sexp x = user[name];
  if (x == R_NilValue) {
    if (previous.size() == 0) {
      cpp11::stop("Expected a value for '%s'", name);
    }
    return previous;
  }

  user_check_array_rank<N>(x, name);
  user_set_array_dim<N>(x, name, dim);

  return user_get_array_value<T>(x, name, min, max);
}

template <>
inline std::vector<float> user_get_array_value(cpp11::sexp x, const char * name,
                                               float min, float max) {
  // NOTE: possible under/overflow here for min/max because we've
  // downcast this.
  std::vector<double> value = user_get_array_value<double>(x, name, min, max);
  std::vector<float> ret(value.size());
  std::copy(value.begin(), value.end(), ret.begin());
  return ret;
}

// This is sum with inclusive "from", exclusive "to", following the
// same function in odin
template <typename real_type, typename container>
__host__ __device__
real_type odin_sum1(const container x, size_t from, size_t to) {
  real_type tot = 0.0;
  for (size_t i = from; i < to; ++i) {
    tot += x[i];
  }
  return tot;
}

inline cpp11::writable::integers integer_sequence(size_t from, size_t len) {
  cpp11::writable::integers ret(len);
  int* data = INTEGER(ret);
  for (size_t i = 0, j = from; i < len; ++i, ++j) {
    data[i] = j;
  }
  return ret;
}
namespace dust {
template<>
dust::pars_type<InfectionSplitModelOD> dust_pars<InfectionSplitModelOD>(cpp11::list user) {
  using real_type = typename InfectionSplitModelOD::real_type;
  auto shared = std::make_shared<InfectionSplitModelOD::shared_type>();
  InfectionSplitModelOD::internal_type internal;
  shared->FOI_max = static_cast<real_type>(0.5);
  shared->initial_day = 0;
  shared->initial_FOI_urb = 0;
  shared->initial_time = 0;
  shared->Pmin = 0;
  shared->t_incubation = 5;
  shared->t_infectious = 5;
  shared->t_latent = 5;
  shared->dt = NA_REAL;
  shared->FOI_spillover = NA_REAL;
  shared->N_age = NA_INTEGER;
  shared->n_years = NA_INTEGER;
  shared->R0 = NA_REAL;
  shared->vaccine_efficacy = NA_REAL;
  shared->year0 = NA_REAL;
  shared->dt = user_get_scalar<real_type>(user, "dt", shared->dt, NA_REAL, NA_REAL);
  shared->FOI_spillover = user_get_scalar<real_type>(user, "FOI_spillover", shared->FOI_spillover, NA_REAL, NA_REAL);
  shared->N_age = user_get_scalar<int>(user, "N_age", shared->N_age, NA_INTEGER, NA_INTEGER);
  shared->n_years = user_get_scalar<int>(user, "n_years", shared->n_years, NA_INTEGER, NA_INTEGER);
  shared->R0 = user_get_scalar<real_type>(user, "R0", shared->R0, NA_REAL, NA_REAL);
  shared->vaccine_efficacy = user_get_scalar<real_type>(user, "vaccine_efficacy", shared->vaccine_efficacy, NA_REAL, NA_REAL);
  shared->year0 = user_get_scalar<real_type>(user, "year0", shared->year0, NA_REAL, NA_REAL);
  shared->beta = (shared->R0 * shared->dt) / (real_type) shared->t_infectious;
  shared->dim_C_sylvatic = shared->N_age;
  shared->dim_C_urban = shared->N_age;
  shared->dim_Cas0 = shared->N_age;
  shared->dim_dP1 = shared->N_age;
  shared->dim_dP1_all_1 = shared->N_age;
  shared->dim_dP1_all_2 = shared->n_years;
  shared->dim_dP2 = shared->N_age;
  shared->dim_dP2_all_1 = shared->N_age;
  shared->dim_dP2_all_2 = shared->n_years;
  shared->dim_E_new_sylvatic = shared->N_age;
  shared->dim_E_new_urban = shared->N_age;
  shared->dim_E_sylvatic = shared->N_age;
  shared->dim_E_urban = shared->N_age;
  shared->dim_Exp0 = shared->N_age;
  shared->dim_I_new_sylvatic = shared->N_age;
  shared->dim_I_new_urban = shared->N_age;
  shared->dim_I_sylvatic = shared->N_age;
  shared->dim_I_urban = shared->N_age;
  shared->dim_Inf0 = shared->N_age;
  shared->dim_inv_P = shared->N_age;
  shared->dim_inv_P_nV = shared->N_age;
  shared->dim_P = shared->N_age;
  shared->dim_P_nV = shared->N_age;
  shared->dim_R = shared->N_age;
  shared->dim_R_new_sylvatic = shared->N_age;
  shared->dim_R_new_urban = shared->N_age;
  shared->dim_Rec0 = shared->N_age;
  shared->dim_S = shared->N_age;
  shared->dim_Sus0 = shared->N_age;
  shared->dim_V = shared->N_age;
  shared->dim_Vac0 = shared->N_age;
  shared->dim_vacc_rate = shared->N_age;
  shared->dim_vacc_rate_annual_1 = shared->N_age;
  shared->dim_vacc_rate_annual_2 = shared->n_years;
  shared->FOI_sylvatic = dust::math::min(shared->FOI_max, shared->FOI_spillover * shared->dt);
  shared->initial_FOI_sylv = shared->FOI_spillover;
  shared->initial_year = shared->year0 - 1;
  shared->rate1 = shared->dt / (real_type) (shared->t_incubation + shared->t_latent);
  shared->rate2 = shared->dt / (real_type) shared->t_infectious;
  internal.dP1 = std::vector<real_type>(shared->dim_dP1);
  internal.dP2 = std::vector<real_type>(shared->dim_dP2);
  internal.E_new_sylvatic = std::vector<real_type>(shared->dim_E_new_sylvatic);
  internal.E_new_urban = std::vector<real_type>(shared->dim_E_new_urban);
  internal.I_new_sylvatic = std::vector<real_type>(shared->dim_I_new_sylvatic);
  internal.I_new_urban = std::vector<real_type>(shared->dim_I_new_urban);
  shared->initial_C_sylvatic = std::vector<real_type>(shared->dim_C_sylvatic);
  shared->initial_C_urban = std::vector<real_type>(shared->dim_C_urban);
  shared->initial_E_sylvatic = std::vector<real_type>(shared->dim_E_sylvatic);
  shared->initial_E_urban = std::vector<real_type>(shared->dim_E_urban);
  shared->initial_I_sylvatic = std::vector<real_type>(shared->dim_I_sylvatic);
  shared->initial_I_urban = std::vector<real_type>(shared->dim_I_urban);
  shared->initial_R = std::vector<real_type>(shared->dim_R);
  shared->initial_S = std::vector<real_type>(shared->dim_S);
  shared->initial_V = std::vector<real_type>(shared->dim_V);
  internal.inv_P = std::vector<real_type>(shared->dim_inv_P);
  internal.inv_P_nV = std::vector<real_type>(shared->dim_inv_P_nV);
  internal.P = std::vector<real_type>(shared->dim_P);
  internal.P_nV = std::vector<real_type>(shared->dim_P_nV);
  internal.R_new_sylvatic = std::vector<real_type>(shared->dim_R_new_sylvatic);
  internal.R_new_urban = std::vector<real_type>(shared->dim_R_new_urban);
  internal.vacc_rate = std::vector<real_type>(shared->dim_vacc_rate);
  shared->Cas0 = user_get_array_fixed<real_type, 1>(user, "Cas0", shared->Cas0, {shared->dim_Cas0}, NA_REAL, NA_REAL);
  shared->dim_dP1_all = shared->dim_dP1_all_1 * shared->dim_dP1_all_2;
  shared->dim_dP2_all = shared->dim_dP2_all_1 * shared->dim_dP2_all_2;
  shared->dim_vacc_rate_annual = shared->dim_vacc_rate_annual_1 * shared->dim_vacc_rate_annual_2;
  shared->Exp0 = user_get_array_fixed<real_type, 1>(user, "Exp0", shared->Exp0, {shared->dim_Exp0}, NA_REAL, NA_REAL);
  shared->Inf0 = user_get_array_fixed<real_type, 1>(user, "Inf0", shared->Inf0, {shared->dim_Inf0}, NA_REAL, NA_REAL);
  for (int i = 1; i <= shared->N_age; ++i) {
    shared->initial_C_urban[i - 1] = 0;
  }
  for (int i = 1; i <= shared->N_age; ++i) {
    shared->initial_E_urban[i - 1] = 0;
  }
  for (int i = 1; i <= shared->N_age; ++i) {
    shared->initial_I_urban[i - 1] = 0;
  }
  shared->offset_variable_C_sylvatic = shared->dim_E_sylvatic + shared->dim_E_urban + shared->dim_I_sylvatic + shared->dim_I_urban + shared->dim_R + shared->dim_S + shared->dim_V + 5;
  shared->offset_variable_C_urban = shared->dim_C_sylvatic + shared->dim_E_sylvatic + shared->dim_E_urban + shared->dim_I_sylvatic + shared->dim_I_urban + shared->dim_R + shared->dim_S + shared->dim_V + 5;
  shared->offset_variable_E_sylvatic = shared->dim_S + 5;
  shared->offset_variable_E_urban = shared->dim_E_sylvatic + shared->dim_S + 5;
  shared->offset_variable_I_sylvatic = shared->dim_E_sylvatic + shared->dim_E_urban + shared->dim_S + 5;
  shared->offset_variable_I_urban = shared->dim_E_sylvatic + shared->dim_E_urban + shared->dim_I_sylvatic + shared->dim_S + 5;
  shared->offset_variable_R = shared->dim_E_sylvatic + shared->dim_E_urban + shared->dim_I_sylvatic + shared->dim_I_urban + shared->dim_S + 5;
  shared->offset_variable_V = shared->dim_E_sylvatic + shared->dim_E_urban + shared->dim_I_sylvatic + shared->dim_I_urban + shared->dim_R + shared->dim_S + 5;
  shared->Rec0 = user_get_array_fixed<real_type, 1>(user, "Rec0", shared->Rec0, {shared->dim_Rec0}, NA_REAL, NA_REAL);
  shared->Sus0 = user_get_array_fixed<real_type, 1>(user, "Sus0", shared->Sus0, {shared->dim_Sus0}, NA_REAL, NA_REAL);
  shared->Vac0 = user_get_array_fixed<real_type, 1>(user, "Vac0", shared->Vac0, {shared->dim_Vac0}, NA_REAL, NA_REAL);
  shared->dP1_all = user_get_array_fixed<real_type, 2>(user, "dP1_all", shared->dP1_all, {shared->dim_dP1_all_1, shared->dim_dP1_all_2}, NA_REAL, NA_REAL);
  shared->dP2_all = user_get_array_fixed<real_type, 2>(user, "dP2_all", shared->dP2_all, {shared->dim_dP2_all_1, shared->dim_dP2_all_2}, NA_REAL, NA_REAL);
  for (int i = 1; i <= shared->N_age; ++i) {
    shared->initial_C_sylvatic[i - 1] = shared->Cas0[i - 1];
  }
  for (int i = 1; i <= shared->N_age; ++i) {
    shared->initial_E_sylvatic[i - 1] = shared->Exp0[i - 1];
  }
  for (int i = 1; i <= shared->N_age; ++i) {
    shared->initial_I_sylvatic[i - 1] = shared->Inf0[i - 1];
  }
  for (int i = 1; i <= shared->N_age; ++i) {
    shared->initial_R[i - 1] = shared->Rec0[i - 1];
  }
  for (int i = 1; i <= shared->N_age; ++i) {
    shared->initial_S[i - 1] = shared->Sus0[i - 1];
  }
  for (int i = 1; i <= shared->N_age; ++i) {
    shared->initial_V[i - 1] = shared->Vac0[i - 1];
  }
  shared->vacc_rate_annual = user_get_array_fixed<real_type, 2>(user, "vacc_rate_annual", shared->vacc_rate_annual, {shared->dim_vacc_rate_annual_1, shared->dim_vacc_rate_annual_2}, NA_REAL, NA_REAL);
  return dust::pars_type<InfectionSplitModelOD>(shared, internal);
}
template <>
cpp11::sexp dust_info<InfectionSplitModelOD>(const dust::pars_type<InfectionSplitModelOD>& pars) {
  const std::shared_ptr<const InfectionSplitModelOD::shared_type> shared = pars.shared;
  cpp11::writable::strings nms({"time", "day", "year", "FOI_sylv", "FOI_urb", "S", "E_sylvatic", "E_urban", "I_sylvatic", "I_urban", "R", "V", "C_sylvatic", "C_urban"});
  cpp11::writable::list dim(14);
  dim[0] = cpp11::writable::integers({1});
  dim[1] = cpp11::writable::integers({1});
  dim[2] = cpp11::writable::integers({1});
  dim[3] = cpp11::writable::integers({1});
  dim[4] = cpp11::writable::integers({1});
  dim[5] = cpp11::writable::integers({shared->dim_S});
  dim[6] = cpp11::writable::integers({shared->dim_E_sylvatic});
  dim[7] = cpp11::writable::integers({shared->dim_E_urban});
  dim[8] = cpp11::writable::integers({shared->dim_I_sylvatic});
  dim[9] = cpp11::writable::integers({shared->dim_I_urban});
  dim[10] = cpp11::writable::integers({shared->dim_R});
  dim[11] = cpp11::writable::integers({shared->dim_V});
  dim[12] = cpp11::writable::integers({shared->dim_C_sylvatic});
  dim[13] = cpp11::writable::integers({shared->dim_C_urban});
  dim.names() = nms;
  cpp11::writable::list index(14);
  index[0] = cpp11::writable::integers({1});
  index[1] = cpp11::writable::integers({2});
  index[2] = cpp11::writable::integers({3});
  index[3] = cpp11::writable::integers({4});
  index[4] = cpp11::writable::integers({5});
  index[5] = integer_sequence(6, shared->dim_S);
  index[6] = integer_sequence(shared->offset_variable_E_sylvatic + 1, shared->dim_E_sylvatic);
  index[7] = integer_sequence(shared->offset_variable_E_urban + 1, shared->dim_E_urban);
  index[8] = integer_sequence(shared->offset_variable_I_sylvatic + 1, shared->dim_I_sylvatic);
  index[9] = integer_sequence(shared->offset_variable_I_urban + 1, shared->dim_I_urban);
  index[10] = integer_sequence(shared->offset_variable_R + 1, shared->dim_R);
  index[11] = integer_sequence(shared->offset_variable_V + 1, shared->dim_V);
  index[12] = integer_sequence(shared->offset_variable_C_sylvatic + 1, shared->dim_C_sylvatic);
  index[13] = integer_sequence(shared->offset_variable_C_urban + 1, shared->dim_C_urban);
  index.names() = nms;
  size_t len = shared->offset_variable_C_urban + shared->dim_C_urban;
  using namespace cpp11::literals;
  return cpp11::writable::list({
           "dim"_nm = dim,
           "len"_nm = len,
           "index"_nm = index});
}
}

cpp11::sexp dust_InfectionSplitModelOD_capabilities() {
  return dust::r::dust_capabilities<InfectionSplitModelOD>();
}

cpp11::sexp dust_InfectionSplitModelOD_gpu_info() {
  return dust::gpu::r::gpu_info();
}
using model_cpu = dust::dust_cpu<InfectionSplitModelOD>;

SEXP dust_cpu_InfectionSplitModelOD_alloc(cpp11::list r_pars, bool pars_multi, cpp11::sexp r_time,
                             cpp11::sexp r_n_particles, int n_threads,
                             cpp11::sexp r_seed, bool deterministic,
                             cpp11::sexp gpu_config, cpp11::sexp ode_control) {
  return dust::r::dust_cpu_alloc<InfectionSplitModelOD>(r_pars, pars_multi, r_time, r_n_particles,
                                        n_threads, r_seed, deterministic,
                                        gpu_config, ode_control);
}

SEXP dust_cpu_InfectionSplitModelOD_run(SEXP ptr, cpp11::sexp r_time_end) {
  return dust::r::dust_run<model_cpu>(ptr, r_time_end);
}

SEXP dust_cpu_InfectionSplitModelOD_simulate(SEXP ptr, cpp11::sexp r_time_end) {
  return dust::r::dust_simulate<model_cpu>(ptr, r_time_end);
}

SEXP dust_cpu_InfectionSplitModelOD_set_index(SEXP ptr, cpp11::sexp r_index) {
  dust::r::dust_set_index<model_cpu>(ptr, r_index);
  return R_NilValue;
}

SEXP dust_cpu_InfectionSplitModelOD_update_state(SEXP ptr, SEXP r_pars, SEXP r_state,
                                           SEXP r_time, SEXP r_set_initial_state, SEXP index, SEXP reset_step_size) {
  return dust::r::dust_update_state<model_cpu>(ptr, r_pars, r_state, r_time,
                                                      r_set_initial_state, index, reset_step_size);
}

SEXP dust_cpu_InfectionSplitModelOD_state(SEXP ptr, SEXP r_index) {
  return dust::r::dust_state<model_cpu>(ptr, r_index);
}

SEXP dust_cpu_InfectionSplitModelOD_time(SEXP ptr) {
  return dust::r::dust_time<model_cpu>(ptr);
}

void dust_cpu_InfectionSplitModelOD_reorder(SEXP ptr, cpp11::sexp r_index) {
  return dust::r::dust_reorder<model_cpu>(ptr, r_index);
}

SEXP dust_cpu_InfectionSplitModelOD_resample(SEXP ptr, cpp11::doubles r_weights) {
  return dust::r::dust_resample<model_cpu>(ptr, r_weights);
}

SEXP dust_cpu_InfectionSplitModelOD_rng_state(SEXP ptr, bool first_only, bool last_only) {
  return dust::r::dust_rng_state<model_cpu>(ptr, first_only, last_only);
}

SEXP dust_cpu_InfectionSplitModelOD_set_rng_state(SEXP ptr, cpp11::raws rng_state) {
  dust::r::dust_set_rng_state<model_cpu>(ptr, rng_state);
  return R_NilValue;
}

SEXP dust_cpu_InfectionSplitModelOD_set_data(SEXP ptr, cpp11::list data,
                                       bool shared) {
  dust::r::dust_set_data<model_cpu>(ptr, data, shared);
  return R_NilValue;
}

SEXP dust_cpu_InfectionSplitModelOD_compare_data(SEXP ptr) {
  return dust::r::dust_compare_data<model_cpu>(ptr);
}

SEXP dust_cpu_InfectionSplitModelOD_filter(SEXP ptr, SEXP time_end,
                                     bool save_trajectories,
                                     cpp11::sexp time_snapshot,
                                     cpp11::sexp min_log_likelihood) {
  return dust::r::dust_filter<model_cpu>(ptr, time_end,
                                                save_trajectories,
                                                time_snapshot,
                                                min_log_likelihood);
}

void dust_cpu_InfectionSplitModelOD_set_n_threads(SEXP ptr, int n_threads) {
  return dust::r::dust_set_n_threads<model_cpu>(ptr, n_threads);
}

int dust_cpu_InfectionSplitModelOD_n_state(SEXP ptr) {
  return dust::r::dust_n_state<model_cpu>(ptr);
}

void dust_cpu_InfectionSplitModelOD_set_stochastic_schedule(SEXP ptr, SEXP time) {
  dust::r::dust_set_stochastic_schedule<model_cpu>(ptr, time);
}

void dust_cpu_InfectionSplitModelOD_ode_statistics(SEXP ptr) {
  dust::r::dust_ode_statistics<model_cpu>(ptr);
}
