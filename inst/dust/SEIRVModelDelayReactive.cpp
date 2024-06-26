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
// [[dust::class(SEIRVModelDelayReactive)]]
// [[dust::time_type(discrete)]]
// [[dust::param(case_threshold, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(cluster_threshold, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(dP1_all, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(dP2_all, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(dt, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(E_0, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(E_delay0, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(FOI_spillover, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(I_0, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(I_delay0, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(N_age, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(n_years, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(p_rep, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(R_0, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(R0, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(response_delay, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(S_0, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(t_cam, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(t_incubation, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(t_infectious, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(t_latent, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(V_0, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(vacc_cov_cam, has_default = FALSE, default_value = NULL, rank = 1, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(vacc_rate_daily, has_default = FALSE, default_value = NULL, rank = 2, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(vaccine_efficacy, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
// [[dust::param(year0, has_default = FALSE, default_value = NULL, rank = 0, min = -Inf, max = Inf, integer = FALSE)]]
class SEIRVModelDelayReactive {
public:
  using real_type = double;
  using rng_state_type = dust::random::generator<real_type>;
  using data_type = dust::no_data;
  struct shared_type {
    real_type beta;
    real_type case_threshold;
    real_type cluster_threshold;
    real_type di1;
    real_type di2;
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
    int dim_E_delay;
    int dim_E_delay0;
    int dim_E_new;
    int dim_I;
    int dim_I_0;
    int dim_I_delay;
    int dim_I_delay0;
    int dim_I_new;
    int dim_inv_P;
    int dim_inv_P_nV;
    int dim_P;
    int dim_P_nV;
    int dim_p_rep;
    int dim_R;
    int dim_R_0;
    int dim_R_new;
    int dim_S;
    int dim_S_0;
    int dim_V;
    int dim_V_0;
    int dim_vacc_cov_cam;
    int dim_vacc_rate;
    int dim_vacc_rate_cam;
    int dim_vacc_rate_daily;
    int dim_vacc_rate_daily_1;
    int dim_vacc_rate_daily_2;
    std::vector<real_type> dP1_all;
    std::vector<real_type> dP2_all;
    real_type dt;
    std::vector<real_type> E_0;
    std::vector<real_type> E_delay0;
    real_type FOI_max;
    real_type FOI_spillover;
    std::vector<real_type> I_0;
    std::vector<real_type> I_delay0;
    std::vector<real_type> initial_C;
    real_type initial_C_rep_total;
    std::vector<real_type> initial_E;
    std::vector<real_type> initial_E_delay;
    real_type initial_flag1;
    real_type initial_flag2;
    real_type initial_flag3;
    real_type initial_flag4;
    real_type initial_FOI_total;
    std::vector<real_type> initial_I;
    std::vector<real_type> initial_I_delay;
    std::vector<real_type> initial_R;
    real_type initial_report_rate;
    std::vector<real_type> initial_S;
    real_type initial_time;
    std::vector<real_type> initial_V;
    real_type initial_year;
    int N_age;
    int n_years;
    int np_E_delay;
    int np_I_delay;
    int offset_variable_C;
    int offset_variable_E;
    int offset_variable_E_delay;
    int offset_variable_I;
    int offset_variable_I_delay;
    int offset_variable_R;
    int offset_variable_V;
    real_type one;
    std::vector<real_type> p_rep;
    real_type Pmin;
    std::vector<real_type> R_0;
    real_type R0;
    real_type rate3;
    real_type rate4;
    real_type response_delay;
    std::vector<real_type> S_0;
    real_type t_cam;
    real_type t_incubation;
    real_type t_infectious;
    real_type t_latent;
    std::vector<real_type> V_0;
    std::vector<real_type> vacc_cov_cam;
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
    std::vector<real_type> vacc_rate_cam;
  };
  SEIRVModelDelayReactive(const dust::pars_type<SEIRVModelDelayReactive>& pars) :
    shared(pars.shared), internal(pars.internal) {
  }
  size_t size() const {
    return shared->dim_C + shared->dim_E + shared->dim_E_delay + shared->dim_I + shared->dim_I_delay + shared->dim_R + shared->dim_S + shared->dim_V + 9;
  }
  std::vector<real_type> initial(size_t step, rng_state_type& rng_state) {
    std::vector<real_type> state(shared->dim_C + shared->dim_E + shared->dim_E_delay + shared->dim_I + shared->dim_I_delay + shared->dim_R + shared->dim_S + shared->dim_V + 9);
    state[0] = shared->initial_time;
    state[1] = shared->initial_year;
    state[2] = shared->initial_FOI_total;
    state[3] = shared->initial_C_rep_total;
    state[4] = shared->initial_flag1;
    state[5] = shared->initial_flag2;
    state[6] = shared->initial_flag3;
    state[7] = shared->initial_flag4;
    state[8] = shared->initial_report_rate;
    std::copy(shared->initial_S.begin(), shared->initial_S.end(), state.begin() + 9);
    std::copy(shared->initial_E.begin(), shared->initial_E.end(), state.begin() + shared->offset_variable_E);
    std::copy(shared->initial_E_delay.begin(), shared->initial_E_delay.end(), state.begin() + shared->offset_variable_E_delay);
    std::copy(shared->initial_I.begin(), shared->initial_I.end(), state.begin() + shared->offset_variable_I);
    std::copy(shared->initial_I_delay.begin(), shared->initial_I_delay.end(), state.begin() + shared->offset_variable_I_delay);
    std::copy(shared->initial_R.begin(), shared->initial_R.end(), state.begin() + shared->offset_variable_R);
    std::copy(shared->initial_V.begin(), shared->initial_V.end(), state.begin() + shared->offset_variable_V);
    std::copy(shared->initial_C.begin(), shared->initial_C.end(), state.begin() + shared->offset_variable_C);
    return state;
  }
  void update(size_t step, const real_type * state, rng_state_type& rng_state, real_type * state_next) {
    const real_type time = state[0];
    const real_type C_rep_total = state[3];
    const real_type flag1 = state[4];
    const real_type flag2 = state[5];
    const real_type flag3 = state[6];
    const real_type flag4 = state[7];
    const real_type * S = state + 9;
    const real_type * E = state + shared->offset_variable_E;
    const real_type * E_delay = state + shared->offset_variable_E_delay;
    const real_type * I = state + shared->offset_variable_I;
    const real_type * I_delay = state + shared->offset_variable_I_delay;
    const real_type * R = state + shared->offset_variable_R;
    const real_type * V = state + shared->offset_variable_V;
    real_type case_flag = (C_rep_total >= shared->case_threshold ? 1 : 0);
    state_next[6] = (flag1 == 1 ? 1 : ((flag2 == 1 ? 1 : 0)));
    state_next[0] = time + shared->dt;
    real_type year_i = dust::math::floor(((step + 1) * shared->dt) / (real_type) 365) + 1;
    real_type p_rep_cur = (flag3 == 1 ? shared->p_rep[1] : shared->p_rep[0]);
    state_next[4] = dust::math::min(shared->one, flag1 + (case_flag * shared->rate3));
    state_next[7] = (flag3 == 0 ? 0 : dust::math::min(shared->one, flag4 + shared->rate4));
    state_next[1] = year_i + shared->year0 - 1;
    for (int i = 1; i <= shared->N_age; ++i) {
      internal.I_new[i - 1] = E_delay[static_cast<int>(i + shared->di1) - 1];
    }
    for (int i = 1; i <= shared->N_age; ++i) {
      internal.P_nV[i - 1] = S[i - 1] + R[i - 1];
    }
    for (int i = 1; i <= shared->N_age; ++i) {
      internal.R_new[i - 1] = I_delay[static_cast<int>(i + shared->di2) - 1];
    }
    state_next[8] = p_rep_cur;
    real_type C_rep_new = dust::random::binomial<real_type>(rng_state, static_cast<int>(odin_sum1<real_type>(internal.I_new.data(), 0, shared->dim_I_new)), p_rep_cur);
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
    for (int i = (shared->N_age + 1); i <= shared->np_I_delay; ++i) {
      state_next[shared->offset_variable_I_delay + i - 1] = I_delay[i - shared->N_age - 1];
    }
    for (int i = 1; i <= shared->N_age; ++i) {
      state_next[shared->offset_variable_I_delay + i - 1] = internal.I_new[i - 1];
    }
    for (int i = 1; i <= shared->N_age; ++i) {
      internal.inv_P[i - 1] = 1 / (real_type) internal.P[i - 1];
    }
    real_type P_tot = odin_sum1<real_type>(internal.P.data(), 0, shared->dim_P);
    state_next[3] = C_rep_total + C_rep_new;
    real_type F_I_total = odin_sum1<real_type>(I, 0, shared->dim_I) / (real_type) P_tot;
    real_type FOI_sum = dust::math::min(shared->FOI_max, shared->beta * (odin_sum1<real_type>(I, 0, shared->dim_I) / (real_type) P_tot) + (shared->FOI_spillover * shared->dt));
    for (int i = 1; i <= shared->N_age; ++i) {
      internal.vacc_rate_cam[i - 1] = (flag3 == 0 ? 0 : (shared->vacc_cov_cam[i - 1] * (1 - (V[i - 1] * internal.inv_P[i - 1]))) / (real_type) shared->t_cam);
    }
    real_type cluster_flag = (F_I_total >= shared->cluster_threshold ? 1 : 0);
    for (int i = 1; i <= shared->N_age; ++i) {
      internal.E_new[i - 1] = dust::random::binomial<real_type>(rng_state, static_cast<int>(S[i - 1]), FOI_sum);
    }
    state_next[2] = FOI_sum;
    for (int i = 1; i <= shared->N_age; ++i) {
      internal.vacc_rate[i - 1] = (shared->vacc_rate_daily[shared->dim_vacc_rate_daily_1 * (static_cast<int>(year_i) - 1) + i - 1] + ((flag3 == 0 ? 0 : internal.vacc_rate_cam[i - 1] * dust::math::ceil(1 - flag4)))) * shared->vaccine_efficacy * shared->dt * internal.P[i - 1];
    }
    for (int i = 1; i <= shared->N_age; ++i) {
      state_next[shared->offset_variable_E + i - 1] = dust::math::max(shared->Pmin, E[i - 1] + internal.E_new[i - 1] - internal.I_new[i - 1]);
    }
    for (int i = (shared->N_age + 1); i <= shared->np_E_delay; ++i) {
      state_next[shared->offset_variable_E_delay + i - 1] = E_delay[i - shared->N_age - 1];
    }
    for (int i = 1; i <= shared->N_age; ++i) {
      state_next[shared->offset_variable_E_delay + i - 1] = internal.E_new[i - 1];
    }
    state_next[5] = dust::math::min(shared->one, flag2 + (cluster_flag * shared->rate3));
    {
       int i = 1;
       state_next[shared->offset_variable_R + i - 1] = dust::math::max(shared->Pmin, R[0] + internal.R_new[0] - internal.vacc_rate[0] * R[0] * internal.inv_P_nV[0] - (internal.dP2[0] * R[0] * internal.inv_P[0]));
    }
    for (int i = 2; i <= shared->N_age; ++i) {
      state_next[shared->offset_variable_R + i - 1] = dust::math::max(shared->Pmin, R[i - 1] + internal.R_new[i - 1] - internal.vacc_rate[i - 1] * R[i - 1] * internal.inv_P_nV[i - 1] + (internal.dP1[i - 1] * R[i - 1 - 1] * internal.inv_P[i - 1 - 1]) - (internal.dP2[i - 1] * R[i - 1] * internal.inv_P[i - 1]));
    }
    {
       int i = 1;
       state_next[9 + i - 1] = dust::math::max(shared->Pmin, S[0] - internal.E_new[0] - internal.vacc_rate[0] * S[0] * internal.inv_P_nV[0] + internal.dP1[0] - (internal.dP2[0] * S[0] * internal.inv_P[0]));
    }
    for (int i = 2; i <= shared->N_age; ++i) {
      state_next[9 + i - 1] = dust::math::max(shared->Pmin, S[i - 1] - internal.E_new[i - 1] - internal.vacc_rate[i - 1] * S[i - 1] * internal.inv_P_nV[i - 1] + (internal.dP1[i - 1] * S[i - 1 - 1] * internal.inv_P[i - 1 - 1]) - (internal.dP2[i - 1] * S[i - 1] * internal.inv_P[i - 1]));
    }
    {
       int i = 1;
       state_next[shared->offset_variable_V + i - 1] = dust::math::max(shared->Pmin, V[0] + internal.vacc_rate[0] - (internal.dP2[0] * V[0] * internal.inv_P[0]));
    }
    for (int i = 2; i <= shared->N_age; ++i) {
      state_next[shared->offset_variable_V + i - 1] = dust::math::max(shared->Pmin, V[i - 1] + internal.vacc_rate[i - 1] + (internal.dP1[i - 1] * V[i - 1 - 1] * internal.inv_P[i - 1 - 1]) - (internal.dP2[i - 1] * V[i - 1] * internal.inv_P[i - 1]));
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
dust::pars_type<SEIRVModelDelayReactive> dust_pars<SEIRVModelDelayReactive>(cpp11::list user) {
  using real_type = typename SEIRVModelDelayReactive::real_type;
  auto shared = std::make_shared<SEIRVModelDelayReactive::shared_type>();
  SEIRVModelDelayReactive::internal_type internal;
  shared->dim_p_rep = 2;
  shared->FOI_max = 1;
  shared->initial_C_rep_total = 0;
  shared->initial_flag1 = 0;
  shared->initial_flag2 = 0;
  shared->initial_flag3 = 0;
  shared->initial_flag4 = 0;
  shared->initial_time = 0;
  shared->one = 1;
  shared->Pmin = static_cast<real_type>(1e-99);
  shared->case_threshold = NA_REAL;
  shared->cluster_threshold = NA_REAL;
  shared->dt = NA_REAL;
  shared->FOI_spillover = NA_REAL;
  shared->N_age = NA_INTEGER;
  shared->n_years = NA_INTEGER;
  shared->R0 = NA_REAL;
  shared->response_delay = NA_REAL;
  shared->t_cam = NA_REAL;
  shared->t_incubation = NA_REAL;
  shared->t_infectious = NA_REAL;
  shared->t_latent = NA_REAL;
  shared->vaccine_efficacy = NA_REAL;
  shared->year0 = NA_REAL;
  shared->case_threshold = user_get_scalar<real_type>(user, "case_threshold", shared->case_threshold, NA_REAL, NA_REAL);
  shared->cluster_threshold = user_get_scalar<real_type>(user, "cluster_threshold", shared->cluster_threshold, NA_REAL, NA_REAL);
  shared->dt = user_get_scalar<real_type>(user, "dt", shared->dt, NA_REAL, NA_REAL);
  shared->FOI_spillover = user_get_scalar<real_type>(user, "FOI_spillover", shared->FOI_spillover, NA_REAL, NA_REAL);
  shared->N_age = user_get_scalar<int>(user, "N_age", shared->N_age, NA_INTEGER, NA_INTEGER);
  shared->n_years = user_get_scalar<int>(user, "n_years", shared->n_years, NA_INTEGER, NA_INTEGER);
  shared->R0 = user_get_scalar<real_type>(user, "R0", shared->R0, NA_REAL, NA_REAL);
  shared->response_delay = user_get_scalar<real_type>(user, "response_delay", shared->response_delay, NA_REAL, NA_REAL);
  shared->t_cam = user_get_scalar<real_type>(user, "t_cam", shared->t_cam, NA_REAL, NA_REAL);
  shared->t_incubation = user_get_scalar<real_type>(user, "t_incubation", shared->t_incubation, NA_REAL, NA_REAL);
  shared->t_infectious = user_get_scalar<real_type>(user, "t_infectious", shared->t_infectious, NA_REAL, NA_REAL);
  shared->t_latent = user_get_scalar<real_type>(user, "t_latent", shared->t_latent, NA_REAL, NA_REAL);
  shared->vaccine_efficacy = user_get_scalar<real_type>(user, "vaccine_efficacy", shared->vaccine_efficacy, NA_REAL, NA_REAL);
  shared->year0 = user_get_scalar<real_type>(user, "year0", shared->year0, NA_REAL, NA_REAL);
  shared->beta = (shared->R0 * shared->dt) / (real_type) shared->t_infectious;
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
  shared->dim_S = shared->N_age;
  shared->dim_S_0 = shared->N_age;
  shared->dim_V = shared->N_age;
  shared->dim_V_0 = shared->N_age;
  shared->dim_vacc_cov_cam = shared->N_age;
  shared->dim_vacc_rate = shared->N_age;
  shared->dim_vacc_rate_cam = shared->N_age;
  shared->dim_vacc_rate_daily_1 = shared->N_age;
  shared->dim_vacc_rate_daily_2 = shared->n_years;
  shared->initial_FOI_total = shared->FOI_spillover;
  shared->initial_year = shared->year0 - 1;
  shared->np_E_delay = ((shared->t_incubation + shared->t_latent) / (real_type) shared->dt) * shared->N_age;
  shared->np_I_delay = (shared->t_infectious / (real_type) shared->dt) * shared->N_age;
  shared->p_rep = user_get_array_fixed<real_type, 1>(user, "p_rep", shared->p_rep, {shared->dim_p_rep}, NA_REAL, NA_REAL);
  shared->rate3 = shared->dt / (real_type) shared->response_delay;
  shared->rate4 = shared->dt / (real_type) shared->t_cam;
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
  internal.vacc_rate_cam = std::vector<real_type>(shared->dim_vacc_rate_cam);
  shared->di1 = shared->np_E_delay - shared->N_age;
  shared->di2 = shared->np_I_delay - shared->N_age;
  shared->dim_dP1_all = shared->dim_dP1_all_1 * shared->dim_dP1_all_2;
  shared->dim_dP2_all = shared->dim_dP2_all_1 * shared->dim_dP2_all_2;
  shared->dim_E_delay = shared->np_E_delay;
  shared->dim_E_delay0 = shared->np_E_delay;
  shared->dim_I_delay = shared->np_I_delay;
  shared->dim_I_delay0 = shared->np_I_delay;
  shared->dim_vacc_rate_daily = shared->dim_vacc_rate_daily_1 * shared->dim_vacc_rate_daily_2;
  shared->E_0 = user_get_array_fixed<real_type, 1>(user, "E_0", shared->E_0, {shared->dim_E_0}, NA_REAL, NA_REAL);
  shared->I_0 = user_get_array_fixed<real_type, 1>(user, "I_0", shared->I_0, {shared->dim_I_0}, NA_REAL, NA_REAL);
  for (int i = 1; i <= shared->N_age; ++i) {
    shared->initial_C[i - 1] = 0;
  }
  shared->initial_report_rate = shared->p_rep[0];
  shared->offset_variable_E = shared->dim_S + 9;
  shared->offset_variable_E_delay = shared->dim_E + shared->dim_S + 9;
  shared->R_0 = user_get_array_fixed<real_type, 1>(user, "R_0", shared->R_0, {shared->dim_R_0}, NA_REAL, NA_REAL);
  shared->S_0 = user_get_array_fixed<real_type, 1>(user, "S_0", shared->S_0, {shared->dim_S_0}, NA_REAL, NA_REAL);
  shared->V_0 = user_get_array_fixed<real_type, 1>(user, "V_0", shared->V_0, {shared->dim_V_0}, NA_REAL, NA_REAL);
  shared->vacc_cov_cam = user_get_array_fixed<real_type, 1>(user, "vacc_cov_cam", shared->vacc_cov_cam, {shared->dim_vacc_cov_cam}, NA_REAL, NA_REAL);
  shared->initial_E_delay = std::vector<real_type>(shared->dim_E_delay);
  shared->initial_I_delay = std::vector<real_type>(shared->dim_I_delay);
  shared->dP1_all = user_get_array_fixed<real_type, 2>(user, "dP1_all", shared->dP1_all, {shared->dim_dP1_all_1, shared->dim_dP1_all_2}, NA_REAL, NA_REAL);
  shared->dP2_all = user_get_array_fixed<real_type, 2>(user, "dP2_all", shared->dP2_all, {shared->dim_dP2_all_1, shared->dim_dP2_all_2}, NA_REAL, NA_REAL);
  shared->E_delay0 = user_get_array_fixed<real_type, 1>(user, "E_delay0", shared->E_delay0, {shared->dim_E_delay0}, NA_REAL, NA_REAL);
  shared->I_delay0 = user_get_array_fixed<real_type, 1>(user, "I_delay0", shared->I_delay0, {shared->dim_I_delay0}, NA_REAL, NA_REAL);
  for (int i = 1; i <= shared->N_age; ++i) {
    shared->initial_E[i - 1] = shared->E_0[i - 1];
  }
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
  shared->offset_variable_C = shared->dim_E + shared->dim_E_delay + shared->dim_I + shared->dim_I_delay + shared->dim_R + shared->dim_S + shared->dim_V + 9;
  shared->offset_variable_I = shared->dim_E + shared->dim_E_delay + shared->dim_S + 9;
  shared->offset_variable_I_delay = shared->dim_E + shared->dim_E_delay + shared->dim_I + shared->dim_S + 9;
  shared->offset_variable_R = shared->dim_E + shared->dim_E_delay + shared->dim_I + shared->dim_I_delay + shared->dim_S + 9;
  shared->offset_variable_V = shared->dim_E + shared->dim_E_delay + shared->dim_I + shared->dim_I_delay + shared->dim_R + shared->dim_S + 9;
  shared->vacc_rate_daily = user_get_array_fixed<real_type, 2>(user, "vacc_rate_daily", shared->vacc_rate_daily, {shared->dim_vacc_rate_daily_1, shared->dim_vacc_rate_daily_2}, NA_REAL, NA_REAL);
  for (int i = 1; i <= shared->np_E_delay; ++i) {
    shared->initial_E_delay[i - 1] = shared->E_delay0[i - 1];
  }
  for (int i = 1; i <= shared->np_I_delay; ++i) {
    shared->initial_I_delay[i - 1] = shared->I_delay0[i - 1];
  }
  return dust::pars_type<SEIRVModelDelayReactive>(shared, internal);
}
template <>
cpp11::sexp dust_info<SEIRVModelDelayReactive>(const dust::pars_type<SEIRVModelDelayReactive>& pars) {
  const std::shared_ptr<const SEIRVModelDelayReactive::shared_type> shared = pars.shared;
  cpp11::writable::strings nms({"time", "year", "FOI_total", "C_rep_total", "flag1", "flag2", "flag3", "flag4", "report_rate", "S", "E", "E_delay", "I", "I_delay", "R", "V", "C"});
  cpp11::writable::list dim(17);
  dim[0] = cpp11::writable::integers({1});
  dim[1] = cpp11::writable::integers({1});
  dim[2] = cpp11::writable::integers({1});
  dim[3] = cpp11::writable::integers({1});
  dim[4] = cpp11::writable::integers({1});
  dim[5] = cpp11::writable::integers({1});
  dim[6] = cpp11::writable::integers({1});
  dim[7] = cpp11::writable::integers({1});
  dim[8] = cpp11::writable::integers({1});
  dim[9] = cpp11::writable::integers({shared->dim_S});
  dim[10] = cpp11::writable::integers({shared->dim_E});
  dim[11] = cpp11::writable::integers({shared->dim_E_delay});
  dim[12] = cpp11::writable::integers({shared->dim_I});
  dim[13] = cpp11::writable::integers({shared->dim_I_delay});
  dim[14] = cpp11::writable::integers({shared->dim_R});
  dim[15] = cpp11::writable::integers({shared->dim_V});
  dim[16] = cpp11::writable::integers({shared->dim_C});
  dim.names() = nms;
  cpp11::writable::list index(17);
  index[0] = cpp11::writable::integers({1});
  index[1] = cpp11::writable::integers({2});
  index[2] = cpp11::writable::integers({3});
  index[3] = cpp11::writable::integers({4});
  index[4] = cpp11::writable::integers({5});
  index[5] = cpp11::writable::integers({6});
  index[6] = cpp11::writable::integers({7});
  index[7] = cpp11::writable::integers({8});
  index[8] = cpp11::writable::integers({9});
  index[9] = integer_sequence(10, shared->dim_S);
  index[10] = integer_sequence(shared->offset_variable_E + 1, shared->dim_E);
  index[11] = integer_sequence(shared->offset_variable_E_delay + 1, shared->dim_E_delay);
  index[12] = integer_sequence(shared->offset_variable_I + 1, shared->dim_I);
  index[13] = integer_sequence(shared->offset_variable_I_delay + 1, shared->dim_I_delay);
  index[14] = integer_sequence(shared->offset_variable_R + 1, shared->dim_R);
  index[15] = integer_sequence(shared->offset_variable_V + 1, shared->dim_V);
  index[16] = integer_sequence(shared->offset_variable_C + 1, shared->dim_C);
  index.names() = nms;
  size_t len = shared->offset_variable_C + shared->dim_C;
  using namespace cpp11::literals;
  return cpp11::writable::list({
           "dim"_nm = dim,
           "len"_nm = len,
           "index"_nm = index});
}
}
