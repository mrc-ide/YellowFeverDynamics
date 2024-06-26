// Generated by dust (version 0.15.2) - do not edit
#include <cpp11.hpp>

[[cpp11::register]]
cpp11::sexp dust_SEIRVModelVarFR_gpu_info();
[[cpp11::register]]
SEXP dust_cpu_SEIRVModelVarFR_alloc(cpp11::list r_pars, bool pars_multi, cpp11::sexp r_time,
                         cpp11::sexp r_n_particles, int n_threads,
                         cpp11::sexp r_seed, bool deterministic,
                         cpp11::sexp gpu_config, cpp11::sexp ode_control);

[[cpp11::register]]
cpp11::sexp dust_cpu_SEIRVModelVarFR_capabilities();

[[cpp11::register]]
SEXP dust_cpu_SEIRVModelVarFR_run(SEXP ptr, cpp11::sexp r_time_end);

[[cpp11::register]]
SEXP dust_cpu_SEIRVModelVarFR_simulate(SEXP ptr, cpp11::sexp time_end);

[[cpp11::register]]
SEXP dust_cpu_SEIRVModelVarFR_run_adjoint(SEXP ptr);

[[cpp11::register]]
SEXP dust_cpu_SEIRVModelVarFR_set_index(SEXP ptr, cpp11::sexp r_index);

[[cpp11::register]]
SEXP dust_cpu_SEIRVModelVarFR_update_state(SEXP ptr, SEXP r_pars, SEXP r_state,
                                           SEXP r_time, SEXP r_set_initial_state,
                                           SEXP index, SEXP reset_step_size);

[[cpp11::register]]
SEXP dust_cpu_SEIRVModelVarFR_state(SEXP ptr, SEXP r_index);

[[cpp11::register]]
SEXP dust_cpu_SEIRVModelVarFR_time(SEXP ptr);

[[cpp11::register]]
void dust_cpu_SEIRVModelVarFR_reorder(SEXP ptr, cpp11::sexp r_index);

[[cpp11::register]]
SEXP dust_cpu_SEIRVModelVarFR_resample(SEXP ptr, cpp11::doubles r_weights);

[[cpp11::register]]
SEXP dust_cpu_SEIRVModelVarFR_rng_state(SEXP ptr, bool first_only, bool last_only);

[[cpp11::register]]
SEXP dust_cpu_SEIRVModelVarFR_set_rng_state(SEXP ptr, cpp11::raws rng_state);

[[cpp11::register]]
SEXP dust_cpu_SEIRVModelVarFR_set_data(SEXP ptr, cpp11::list data, bool shared);

[[cpp11::register]]
SEXP dust_cpu_SEIRVModelVarFR_compare_data(SEXP ptr);

[[cpp11::register]]
SEXP dust_cpu_SEIRVModelVarFR_filter(SEXP ptr, SEXP time_end,
                                     bool save_trajectories,
                                     cpp11::sexp time_snapshot,
                                     cpp11::sexp min_log_likelihood);

[[cpp11::register]]
void dust_cpu_SEIRVModelVarFR_set_n_threads(SEXP ptr, int n_threads);

[[cpp11::register]]
int dust_cpu_SEIRVModelVarFR_n_state(SEXP ptr);

[[cpp11::register]]
void dust_cpu_SEIRVModelVarFR_set_stochastic_schedule(SEXP ptr, SEXP time);

[[cpp11::register]]
SEXP dust_cpu_SEIRVModelVarFR_ode_statistics(SEXP ptr);
#include <dust/r/dust.hpp>

// Generated by odin.dust (version 0.3.11) - do not edit
template <typename real_type, typename container>
__host__ __device__ real_type odin_sum1(const container x, size_t from, size_t to);
template <typename real_type, typename container>
__host__ __device__ real_type odin_sum2(const container x, int from_i, int to_i, int from_j, int to_j, int dim_x_1);
template <typename real_type, typename T, typename U>
__host__ __device__ real_type fintdiv(T x, U y) {
  return std::floor(x / y);
}

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

template <typename T>
__host__ __device__ T odin_sign(T x) {
  return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}
// [[dust::class(SEIRVModelVarFR)]]
// [[dust::time_type(discrete)]]
// [[dust::param(dP1_all, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(dP2_all, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(dt, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(E_0, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(FOI_spillover, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(I_0, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(N_age, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(n_t_pts, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(n_years, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(R_0, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(R0, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(S_0, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(t_incubation, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(t_infectious, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(t_latent, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(V_0, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(vacc_rate_daily, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(vaccine_efficacy, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(year0, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
class SEIRVModelVarFR {
public:
  using real_type = double;
  using rng_state_type = dust::random::generator<real_type>;
  using data_type = dust::no_data;
  struct shared_type {
    int dim_C;
    int dim_dP1;
    int dim_dP1_all;
    int dim_dP1_all_1;
    int dim_dP1_all_2;
    int dim_dP2;
    int dim_dP2_all;
    int dim_dP2_all_1;
    int dim_dP2_all_2;
    int dim_E;
    int dim_E_0;
    int dim_E_new;
    int dim_FOI_spillover;
    int dim_I;
    int dim_I_0;
    int dim_I_new;
    int dim_inv_P;
    int dim_inv_P_nV;
    int dim_P;
    int dim_P_nV;
    int dim_R;
    int dim_R_0;
    int dim_R_new;
    int dim_R0;
    int dim_S;
    int dim_S_0;
    int dim_V;
    int dim_V_0;
    int dim_vacc_rate;
    int dim_vacc_rate_daily;
    int dim_vacc_rate_daily_1;
    int dim_vacc_rate_daily_2;
    std::vector<real_type> dP1_all;
    std::vector<real_type> dP2_all;
    real_type dt;
    std::vector<real_type> E_0;
    real_type FOI_max;
    std::vector<real_type> FOI_spillover;
    std::vector<real_type> I_0;
    std::vector<real_type> initial_C;
    std::vector<real_type> initial_E;
    real_type initial_FOI_total;
    std::vector<real_type> initial_I;
    std::vector<real_type> initial_R;
    std::vector<real_type> initial_S;
    real_type initial_time;
    std::vector<real_type> initial_V;
    real_type initial_year;
    int N_age;
    int n_t_pts;
    int n_years;
    int offset_variable_C;
    int offset_variable_E;
    int offset_variable_I;
    int offset_variable_R;
    int offset_variable_V;
    real_type Pmin;
    std::vector<real_type> R_0;
    std::vector<real_type> R0;
    real_type rate1;
    real_type rate2;
    std::vector<real_type> S_0;
    real_type t_incubation;
    real_type t_infectious;
    real_type t_latent;
    std::vector<real_type> V_0;
    std::vector<real_type> vacc_rate_daily;
    real_type vaccine_efficacy;
    real_type year0;
  };
  struct internal_type {
    std::vector<real_type> dP1;
    std::vector<real_type> dP2;
    std::vector<real_type> E_new;
    std::vector<real_type> I_new;
    std::vector<real_type> inv_P;
    std::vector<real_type> inv_P_nV;
    std::vector<real_type> P;
    std::vector<real_type> P_nV;
    std::vector<real_type> R_new;
    std::vector<real_type> vacc_rate;
  };
  SEIRVModelVarFR(const dust::pars_type<SEIRVModelVarFR>& pars) :
    shared(pars.shared), internal(pars.internal) {
  }
  size_t size() const {
    return shared->dim_C + shared->dim_E + shared->dim_I + shared->dim_R + shared->dim_S + shared->dim_V + 3;
  }
  std::vector<real_type> initial(size_t step, rng_state_type& rng_state) {
    std::vector<real_type> state(shared->dim_C + shared->dim_E + shared->dim_I + shared->dim_R + shared->dim_S + shared->dim_V + 3);
    state[0] = shared->initial_time;
    state[1] = shared->initial_year;
    state[2] = shared->initial_FOI_total;
    std::copy(shared->initial_S.begin(), shared->initial_S.end(), state.begin() + 3);
    std::copy(shared->initial_E.begin(), shared->initial_E.end(), state.begin() + shared->offset_variable_E);
    std::copy(shared->initial_I.begin(), shared->initial_I.end(), state.begin() + shared->offset_variable_I);
    std::copy(shared->initial_R.begin(), shared->initial_R.end(), state.begin() + shared->offset_variable_R);
    std::copy(shared->initial_V.begin(), shared->initial_V.end(), state.begin() + shared->offset_variable_V);
    std::copy(shared->initial_C.begin(), shared->initial_C.end(), state.begin() + shared->offset_variable_C);
    return state;
  }
  void update(size_t step, const real_type * state, rng_state_type& rng_state, real_type * state_next) {
    const real_type time = state[0];
    const real_type * S = state + 3;
    const real_type * E = state + shared->offset_variable_E;
    const real_type * I = state + shared->offset_variable_I;
    const real_type * R = state + shared->offset_variable_R;
    const real_type * V = state + shared->offset_variable_V;
    state_next[0] = time + shared->dt;
    real_type year_i = dust::math::floor(((step + 1) * shared->dt) / (real_type) 365) + 1;
    state_next[1] = year_i + shared->year0 - 1;
    real_type beta = (shared->R0[step + 1 - 1] * shared->dt) / (real_type) shared->t_infectious;
    for (int i = 1; i <= shared->N_age; ++i) {
      internal.I_new[i - 1] = E[i - 1] * shared->rate1;
    }
    for (int i = 1; i <= shared->N_age; ++i) {
      internal.P_nV[i - 1] = S[i - 1] + R[i - 1];
    }
    for (int i = 1; i <= shared->N_age; ++i) {
      internal.R_new[i - 1] = I[i - 1] * shared->rate2;
    }
    for (int i = 1; i <= shared->N_age; ++i) {
      internal.dP1[i - 1] = shared->dP1_all[shared->dim_dP1_all_1 * (static_cast<int>(year_i) - 1) + i - 1] * shared->dt;
    }
    for (int i = 1; i <= shared->N_age; ++i) {
      internal.dP2[i - 1] = shared->dP2_all[shared->dim_dP2_all_1 * (static_cast<int>(year_i) - 1) + i - 1] * shared->dt;
    }
    for (int i = 1; i <= shared->N_age; ++i) {
      internal.inv_P_nV[i - 1] = 1 / (real_type) internal.P_nV[i - 1];
    }
    for (int i = 1; i <= shared->N_age; ++i) {
      internal.P[i - 1] = internal.P_nV[i - 1] + V[i - 1];
    }
    for (int i = 1; i <= shared->N_age; ++i) {
      state_next[shared->offset_variable_C + i - 1] = internal.I_new[i - 1];
    }
    for (int i = 1; i <= shared->N_age; ++i) {
      state_next[shared->offset_variable_I + i - 1] = dust::math::max(shared->Pmin, I[i - 1] + internal.I_new[i - 1] - internal.R_new[i - 1]);
    }
    for (int i = 1; i <= shared->N_age; ++i) {
      internal.inv_P[i - 1] = 1 / (real_type) internal.P[i - 1];
    }
    real_type P_tot = odin_sum1<real_type>(internal.P.data(), 0, shared->dim_P);
    for (int i = 1; i <= shared->N_age; ++i) {
      internal.vacc_rate[i - 1] = shared->vacc_rate_daily[shared->dim_vacc_rate_daily_1 * (static_cast<int>(year_i) - 1) + i - 1] * shared->vaccine_efficacy * shared->dt * internal.P[i - 1];
    }
    real_type FOI_sum = dust::math::min(shared->FOI_max, beta * (odin_sum1<real_type>(I, 0, shared->dim_I) / (real_type) P_tot) + (shared->FOI_spillover[step + 1 - 1] * shared->dt));
    {
       int i = 1;
       state_next[shared->offset_variable_R + i - 1] = dust::math::max(shared->Pmin, R[0] + internal.R_new[0] - internal.vacc_rate[0] * R[0] * internal.inv_P_nV[0] - (internal.dP2[0] * R[0] * internal.inv_P[0]));
    }
    for (int i = 2; i <= shared->N_age; ++i) {
      state_next[shared->offset_variable_R + i - 1] = dust::math::max(shared->Pmin, R[i - 1] + internal.R_new[i - 1] - internal.vacc_rate[i - 1] * R[i - 1] * internal.inv_P_nV[i - 1] + (internal.dP1[i - 1] * R[i - 1 - 1] * internal.inv_P[i - 1 - 1]) - (internal.dP2[i - 1] * R[i - 1] * internal.inv_P[i - 1]));
    }
    {
       int i = 1;
       state_next[shared->offset_variable_V + i - 1] = dust::math::max(shared->Pmin, V[0] + internal.vacc_rate[0] - (internal.dP2[0] * V[0] * internal.inv_P[0]));
    }
    for (int i = 2; i <= shared->N_age; ++i) {
      state_next[shared->offset_variable_V + i - 1] = dust::math::max(shared->Pmin, V[i - 1] + internal.vacc_rate[i - 1] + (internal.dP1[i - 1] * V[i - 1 - 1] * internal.inv_P[i - 1 - 1]) - (internal.dP2[i - 1] * V[i - 1] * internal.inv_P[i - 1]));
    }
    for (int i = 1; i <= shared->N_age; ++i) {
      internal.E_new[i - 1] = dust::random::binomial<real_type>(rng_state, static_cast<int>(S[i - 1]), FOI_sum);
    }
    state_next[2] = FOI_sum;
    for (int i = 1; i <= shared->N_age; ++i) {
      state_next[shared->offset_variable_E + i - 1] = dust::math::max(shared->Pmin, E[i - 1] + internal.E_new[i - 1] - internal.I_new[i - 1]);
    }
    {
       int i = 1;
       state_next[3 + i - 1] = dust::math::max(shared->Pmin, S[0] - internal.E_new[0] - internal.vacc_rate[0] * S[0] * internal.inv_P_nV[0] + internal.dP1[0] - (internal.dP2[0] * S[0] * internal.inv_P[0]));
    }
    for (int i = 2; i <= shared->N_age; ++i) {
      state_next[3 + i - 1] = dust::math::max(shared->Pmin, S[i - 1] - internal.E_new[i - 1] - internal.vacc_rate[i - 1] * S[i - 1] * internal.inv_P_nV[i - 1] + (internal.dP1[i - 1] * S[i - 1 - 1] * internal.inv_P[i - 1 - 1]) - (internal.dP2[i - 1] * S[i - 1] * internal.inv_P[i - 1]));
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
      Rf_error("Incorrect size of dimension %zu of '%s' (expected %d)",
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
dust::pars_type<SEIRVModelVarFR> dust_pars<SEIRVModelVarFR>(cpp11::list user) {
  using real_type = typename SEIRVModelVarFR::real_type;
  auto shared = std::make_shared<SEIRVModelVarFR::shared_type>();
  SEIRVModelVarFR::internal_type internal;
  shared->FOI_max = 1;
  shared->initial_time = 0;
  shared->Pmin = static_cast<real_type>(1e-99);
  shared->dt = NA_REAL;
  shared->N_age = NA_INTEGER;
  shared->n_t_pts = NA_INTEGER;
  shared->n_years = NA_INTEGER;
  shared->t_incubation = NA_REAL;
  shared->t_infectious = NA_REAL;
  shared->t_latent = NA_REAL;
  shared->vaccine_efficacy = NA_REAL;
  shared->year0 = NA_REAL;
  shared->dt = user_get_scalar<real_type>(user, "dt", shared->dt, NA_REAL, NA_REAL);
  shared->N_age = user_get_scalar<int>(user, "N_age", shared->N_age, NA_INTEGER, NA_INTEGER);
  shared->n_t_pts = user_get_scalar<int>(user, "n_t_pts", shared->n_t_pts, NA_INTEGER, NA_INTEGER);
  shared->n_years = user_get_scalar<int>(user, "n_years", shared->n_years, NA_INTEGER, NA_INTEGER);
  shared->t_incubation = user_get_scalar<real_type>(user, "t_incubation", shared->t_incubation, NA_REAL, NA_REAL);
  shared->t_infectious = user_get_scalar<real_type>(user, "t_infectious", shared->t_infectious, NA_REAL, NA_REAL);
  shared->t_latent = user_get_scalar<real_type>(user, "t_latent", shared->t_latent, NA_REAL, NA_REAL);
  shared->vaccine_efficacy = user_get_scalar<real_type>(user, "vaccine_efficacy", shared->vaccine_efficacy, NA_REAL, NA_REAL);
  shared->year0 = user_get_scalar<real_type>(user, "year0", shared->year0, NA_REAL, NA_REAL);
  shared->dim_C = shared->N_age;
  shared->dim_dP1 = shared->N_age;
  shared->dim_dP1_all_1 = shared->N_age;
  shared->dim_dP1_all_2 = shared->n_years;
  shared->dim_dP2 = shared->N_age;
  shared->dim_dP2_all_1 = shared->N_age;
  shared->dim_dP2_all_2 = shared->n_years;
  shared->dim_E = shared->N_age;
  shared->dim_E_0 = shared->N_age;
  shared->dim_E_new = shared->N_age;
  shared->dim_FOI_spillover = shared->n_t_pts;
  shared->dim_I = shared->N_age;
  shared->dim_I_0 = shared->N_age;
  shared->dim_I_new = shared->N_age;
  shared->dim_inv_P = shared->N_age;
  shared->dim_inv_P_nV = shared->N_age;
  shared->dim_P = shared->N_age;
  shared->dim_P_nV = shared->N_age;
  shared->dim_R = shared->N_age;
  shared->dim_R_0 = shared->N_age;
  shared->dim_R_new = shared->N_age;
  shared->dim_R0 = shared->n_t_pts;
  shared->dim_S = shared->N_age;
  shared->dim_S_0 = shared->N_age;
  shared->dim_V = shared->N_age;
  shared->dim_V_0 = shared->N_age;
  shared->dim_vacc_rate = shared->N_age;
  shared->dim_vacc_rate_daily_1 = shared->N_age;
  shared->dim_vacc_rate_daily_2 = shared->n_years;
  shared->initial_year = shared->year0 - 1;
  shared->rate1 = shared->dt / (real_type) (shared->t_incubation + shared->t_latent);
  shared->rate2 = shared->dt / (real_type) shared->t_infectious;
  internal.dP1 = std::vector<real_type>(shared->dim_dP1);
  internal.dP2 = std::vector<real_type>(shared->dim_dP2);
  internal.E_new = std::vector<real_type>(shared->dim_E_new);
  internal.I_new = std::vector<real_type>(shared->dim_I_new);
  shared->initial_C = std::vector<real_type>(shared->dim_C);
  shared->initial_E = std::vector<real_type>(shared->dim_E);
  shared->initial_I = std::vector<real_type>(shared->dim_I);
  shared->initial_R = std::vector<real_type>(shared->dim_R);
  shared->initial_S = std::vector<real_type>(shared->dim_S);
  shared->initial_V = std::vector<real_type>(shared->dim_V);
  internal.inv_P = std::vector<real_type>(shared->dim_inv_P);
  internal.inv_P_nV = std::vector<real_type>(shared->dim_inv_P_nV);
  internal.P = std::vector<real_type>(shared->dim_P);
  internal.P_nV = std::vector<real_type>(shared->dim_P_nV);
  internal.R_new = std::vector<real_type>(shared->dim_R_new);
  internal.vacc_rate = std::vector<real_type>(shared->dim_vacc_rate);
  shared->dim_dP1_all = shared->dim_dP1_all_1 * shared->dim_dP1_all_2;
  shared->dim_dP2_all = shared->dim_dP2_all_1 * shared->dim_dP2_all_2;
  shared->dim_vacc_rate_daily = shared->dim_vacc_rate_daily_1 * shared->dim_vacc_rate_daily_2;
  shared->E_0 = user_get_array_fixed<real_type, 1>(user, "E_0", shared->E_0, {shared->dim_E_0}, NA_REAL, NA_REAL);
  shared->FOI_spillover = user_get_array_fixed<real_type, 1>(user, "FOI_spillover", shared->FOI_spillover, {shared->dim_FOI_spillover}, NA_REAL, NA_REAL);
  shared->I_0 = user_get_array_fixed<real_type, 1>(user, "I_0", shared->I_0, {shared->dim_I_0}, NA_REAL, NA_REAL);
  for (int i = 1; i <= shared->N_age; ++i) {
    shared->initial_C[i - 1] = 0;
  }
  shared->offset_variable_C = shared->dim_E + shared->dim_I + shared->dim_R + shared->dim_S + shared->dim_V + 3;
  shared->offset_variable_E = shared->dim_S + 3;
  shared->offset_variable_I = shared->dim_E + shared->dim_S + 3;
  shared->offset_variable_R = shared->dim_E + shared->dim_I + shared->dim_S + 3;
  shared->offset_variable_V = shared->dim_E + shared->dim_I + shared->dim_R + shared->dim_S + 3;
  shared->R_0 = user_get_array_fixed<real_type, 1>(user, "R_0", shared->R_0, {shared->dim_R_0}, NA_REAL, NA_REAL);
  shared->R0 = user_get_array_fixed<real_type, 1>(user, "R0", shared->R0, {shared->dim_R0}, NA_REAL, NA_REAL);
  shared->S_0 = user_get_array_fixed<real_type, 1>(user, "S_0", shared->S_0, {shared->dim_S_0}, NA_REAL, NA_REAL);
  shared->V_0 = user_get_array_fixed<real_type, 1>(user, "V_0", shared->V_0, {shared->dim_V_0}, NA_REAL, NA_REAL);
  shared->dP1_all = user_get_array_fixed<real_type, 2>(user, "dP1_all", shared->dP1_all, {shared->dim_dP1_all_1, shared->dim_dP1_all_2}, NA_REAL, NA_REAL);
  shared->dP2_all = user_get_array_fixed<real_type, 2>(user, "dP2_all", shared->dP2_all, {shared->dim_dP2_all_1, shared->dim_dP2_all_2}, NA_REAL, NA_REAL);
  for (int i = 1; i <= shared->N_age; ++i) {
    shared->initial_E[i - 1] = shared->E_0[i - 1];
  }
  shared->initial_FOI_total = shared->FOI_spillover[0];
  for (int i = 1; i <= shared->N_age; ++i) {
    shared->initial_I[i - 1] = shared->I_0[i - 1];
  }
  for (int i = 1; i <= shared->N_age; ++i) {
    shared->initial_R[i - 1] = shared->R_0[i - 1];
  }
  for (int i = 1; i <= shared->N_age; ++i) {
    shared->initial_S[i - 1] = shared->S_0[i - 1];
  }
  for (int i = 1; i <= shared->N_age; ++i) {
    shared->initial_V[i - 1] = shared->V_0[i - 1];
  }
  shared->vacc_rate_daily = user_get_array_fixed<real_type, 2>(user, "vacc_rate_daily", shared->vacc_rate_daily, {shared->dim_vacc_rate_daily_1, shared->dim_vacc_rate_daily_2}, NA_REAL, NA_REAL);
  return dust::pars_type<SEIRVModelVarFR>(shared, internal);
}
template <>
cpp11::sexp dust_info<SEIRVModelVarFR>(const dust::pars_type<SEIRVModelVarFR>& pars) {
  const std::shared_ptr<const SEIRVModelVarFR::shared_type> shared = pars.shared;
  cpp11::writable::strings nms({"time", "year", "FOI_total", "S", "E", "I", "R", "V", "C"});
  cpp11::writable::list dim(9);
  dim[0] = cpp11::writable::integers({1});
  dim[1] = cpp11::writable::integers({1});
  dim[2] = cpp11::writable::integers({1});
  dim[3] = cpp11::writable::integers({shared->dim_S});
  dim[4] = cpp11::writable::integers({shared->dim_E});
  dim[5] = cpp11::writable::integers({shared->dim_I});
  dim[6] = cpp11::writable::integers({shared->dim_R});
  dim[7] = cpp11::writable::integers({shared->dim_V});
  dim[8] = cpp11::writable::integers({shared->dim_C});
  dim.names() = nms;
  cpp11::writable::list index(9);
  index[0] = cpp11::writable::integers({1});
  index[1] = cpp11::writable::integers({2});
  index[2] = cpp11::writable::integers({3});
  index[3] = integer_sequence(4, shared->dim_S);
  index[4] = integer_sequence(shared->offset_variable_E + 1, shared->dim_E);
  index[5] = integer_sequence(shared->offset_variable_I + 1, shared->dim_I);
  index[6] = integer_sequence(shared->offset_variable_R + 1, shared->dim_R);
  index[7] = integer_sequence(shared->offset_variable_V + 1, shared->dim_V);
  index[8] = integer_sequence(shared->offset_variable_C + 1, shared->dim_C);
  index.names() = nms;
  size_t len = shared->offset_variable_C + shared->dim_C;
  using namespace cpp11::literals;
  return cpp11::writable::list({
           "dim"_nm = dim,
           "len"_nm = len,
           "index"_nm = index});
}
}

cpp11::sexp dust_SEIRVModelVarFR_gpu_info() {
  return dust::gpu::r::gpu_info();
}
using model_cpu = dust::dust_cpu<SEIRVModelVarFR>;

cpp11::sexp dust_cpu_SEIRVModelVarFR_capabilities() {
  return dust::r::dust_capabilities<model_cpu>();
}

SEXP dust_cpu_SEIRVModelVarFR_alloc(cpp11::list r_pars, bool pars_multi, cpp11::sexp r_time,
                             cpp11::sexp r_n_particles, int n_threads,
                             cpp11::sexp r_seed, bool deterministic,
                             cpp11::sexp gpu_config, cpp11::sexp ode_control) {
  return dust::r::dust_cpu_alloc<SEIRVModelVarFR>(r_pars, pars_multi, r_time, r_n_particles,
                                        n_threads, r_seed, deterministic,
                                        gpu_config, ode_control);
}

SEXP dust_cpu_SEIRVModelVarFR_run(SEXP ptr, cpp11::sexp r_time_end) {
  return dust::r::dust_run<model_cpu>(ptr, r_time_end);
}

SEXP dust_cpu_SEIRVModelVarFR_simulate(SEXP ptr, cpp11::sexp r_time_end) {
  return dust::r::dust_simulate<model_cpu>(ptr, r_time_end);
}

SEXP dust_cpu_SEIRVModelVarFR_run_adjoint(SEXP ptr) {
  return dust::r::dust_run_adjoint<model_cpu>(ptr);
}

SEXP dust_cpu_SEIRVModelVarFR_set_index(SEXP ptr, cpp11::sexp r_index) {
  dust::r::dust_set_index<model_cpu>(ptr, r_index);
  return R_NilValue;
}

SEXP dust_cpu_SEIRVModelVarFR_update_state(SEXP ptr, SEXP r_pars, SEXP r_state,
                                           SEXP r_time, SEXP r_set_initial_state, SEXP index, SEXP reset_step_size) {
  return dust::r::dust_update_state<model_cpu>(ptr, r_pars, r_state, r_time,
                                                      r_set_initial_state, index, reset_step_size);
}

SEXP dust_cpu_SEIRVModelVarFR_state(SEXP ptr, SEXP r_index) {
  return dust::r::dust_state<model_cpu>(ptr, r_index);
}

SEXP dust_cpu_SEIRVModelVarFR_time(SEXP ptr) {
  return dust::r::dust_time<model_cpu>(ptr);
}

void dust_cpu_SEIRVModelVarFR_reorder(SEXP ptr, cpp11::sexp r_index) {
  return dust::r::dust_reorder<model_cpu>(ptr, r_index);
}

SEXP dust_cpu_SEIRVModelVarFR_resample(SEXP ptr, cpp11::doubles r_weights) {
  return dust::r::dust_resample<model_cpu>(ptr, r_weights);
}

SEXP dust_cpu_SEIRVModelVarFR_rng_state(SEXP ptr, bool first_only, bool last_only) {
  return dust::r::dust_rng_state<model_cpu>(ptr, first_only, last_only);
}

SEXP dust_cpu_SEIRVModelVarFR_set_rng_state(SEXP ptr, cpp11::raws rng_state) {
  dust::r::dust_set_rng_state<model_cpu>(ptr, rng_state);
  return R_NilValue;
}

SEXP dust_cpu_SEIRVModelVarFR_set_data(SEXP ptr, cpp11::list data,
                                       bool shared) {
  dust::r::dust_set_data<model_cpu>(ptr, data, shared);
  return R_NilValue;
}

SEXP dust_cpu_SEIRVModelVarFR_compare_data(SEXP ptr) {
  return dust::r::dust_compare_data<model_cpu>(ptr);
}

SEXP dust_cpu_SEIRVModelVarFR_filter(SEXP ptr, SEXP time_end,
                                     bool save_trajectories,
                                     cpp11::sexp time_snapshot,
                                     cpp11::sexp min_log_likelihood) {
  return dust::r::dust_filter<model_cpu>(ptr, time_end,
                                                save_trajectories,
                                                time_snapshot,
                                                min_log_likelihood);
}

void dust_cpu_SEIRVModelVarFR_set_n_threads(SEXP ptr, int n_threads) {
  return dust::r::dust_set_n_threads<model_cpu>(ptr, n_threads);
}

int dust_cpu_SEIRVModelVarFR_n_state(SEXP ptr) {
  return dust::r::dust_n_state<model_cpu>(ptr);
}

void dust_cpu_SEIRVModelVarFR_set_stochastic_schedule(SEXP ptr, SEXP time) {
  dust::r::dust_set_stochastic_schedule<model_cpu>(ptr, time);
}

SEXP dust_cpu_SEIRVModelVarFR_ode_statistics(SEXP ptr) {
  return dust::r::dust_ode_statistics<model_cpu>(ptr);
}
