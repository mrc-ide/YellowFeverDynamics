// Generated by dust2 (version 0.3.10) - do not edit

// Generated by odin2 (version 0.3.8) - do not edit
#include <dust2/common.hpp>
// [[dust2::class(SEIRVModelDelay)]]
// [[dust2::time_type(discrete)]]
// [[dust2::parameter(time_inc, type = "real_type", rank = 0, required = TRUE, constant = FALSE)]]
// [[dust2::parameter(t_infectious, type = "real_type", rank = 0, required = TRUE, constant = FALSE)]]
// [[dust2::parameter(FOI_spillover, type = "real_type", rank = 0, required = TRUE, constant = FALSE)]]
// [[dust2::parameter(R0, type = "real_type", rank = 0, required = TRUE, constant = FALSE)]]
// [[dust2::parameter(N_age, type = "int", rank = 0, required = TRUE, constant = TRUE)]]
// [[dust2::parameter(vacc_rate_daily, type = "real_type", rank = 2, required = TRUE, constant = FALSE)]]
// [[dust2::parameter(vaccine_efficacy, type = "real_type", rank = 0, required = TRUE, constant = FALSE)]]
// [[dust2::parameter(year0, type = "real_type", rank = 0, required = TRUE, constant = FALSE)]]
// [[dust2::parameter(S_0, type = "real_type", rank = 1, required = TRUE, constant = FALSE)]]
// [[dust2::parameter(E_0, type = "real_type", rank = 1, required = TRUE, constant = FALSE)]]
// [[dust2::parameter(E_delay0, type = "real_type", rank = 1, required = TRUE, constant = FALSE)]]
// [[dust2::parameter(I_0, type = "real_type", rank = 1, required = TRUE, constant = FALSE)]]
// [[dust2::parameter(I_delay0, type = "real_type", rank = 1, required = TRUE, constant = FALSE)]]
// [[dust2::parameter(R_0, type = "real_type", rank = 1, required = TRUE, constant = FALSE)]]
// [[dust2::parameter(V_0, type = "real_type", rank = 1, required = TRUE, constant = FALSE)]]
// [[dust2::parameter(dP1_all, type = "real_type", rank = 2, required = TRUE, constant = FALSE)]]
// [[dust2::parameter(dP2_all, type = "real_type", rank = 2, required = TRUE, constant = FALSE)]]
// [[dust2::parameter(n_years, type = "int", rank = 0, required = TRUE, constant = TRUE)]]
// [[dust2::parameter(np_E_delay, type = "int", rank = 0, required = TRUE, constant = TRUE)]]
// [[dust2::parameter(np_I_delay, type = "int", rank = 0, required = TRUE, constant = TRUE)]]
class SEIRVModelDelay {
public:
  SEIRVModelDelay() = delete;
  using real_type = double;
  using rng_state_type = monty::random::generator<real_type>;
  struct shared_state {
    struct odin_internals_type {
      struct {
        dust2::packing state;
      } packing;
      struct {
        std::array<size_t, 11> state;
      } offset;
    } odin;
    struct dim_type {
      dust2::array::dimensions<1> S;
      dust2::array::dimensions<1> E;
      dust2::array::dimensions<1> E_delay;
      dust2::array::dimensions<1> I;
      dust2::array::dimensions<1> I_delay;
      dust2::array::dimensions<1> R;
      dust2::array::dimensions<1> V;
      dust2::array::dimensions<1> C;
      dust2::array::dimensions<1> dP1;
      dust2::array::dimensions<1> dP2;
      dust2::array::dimensions<1> E_new;
      dust2::array::dimensions<1> I_new;
      dust2::array::dimensions<1> R_new;
      dust2::array::dimensions<1> P_nV;
      dust2::array::dimensions<1> inv_P_nV;
      dust2::array::dimensions<1> P;
      dust2::array::dimensions<1> inv_P;
      dust2::array::dimensions<1> vacc_rate;
      dust2::array::dimensions<1> S_0;
      dust2::array::dimensions<1> E_0;
      dust2::array::dimensions<1> E_delay0;
      dust2::array::dimensions<1> I_0;
      dust2::array::dimensions<1> I_delay0;
      dust2::array::dimensions<1> R_0;
      dust2::array::dimensions<1> V_0;
      dust2::array::dimensions<2> dP1_all;
      dust2::array::dimensions<2> dP2_all;
      dust2::array::dimensions<2> vacc_rate_daily;
    } dim;
    real_type time_inc;
    real_type t_infectious;
    real_type FOI_spillover;
    real_type R0;
    int N_age;
    real_type vaccine_efficacy;
    real_type year0;
    int n_years;
    real_type Pmin;
    real_type FOI_max;
    int np_E_delay;
    int np_I_delay;
    real_type di1;
    real_type di2;
    real_type beta;
    std::vector<real_type> vacc_rate_daily;
    std::vector<real_type> S_0;
    std::vector<real_type> E_0;
    std::vector<real_type> E_delay0;
    std::vector<real_type> I_0;
    std::vector<real_type> I_delay0;
    std::vector<real_type> R_0;
    std::vector<real_type> V_0;
    std::vector<real_type> dP1_all;
    std::vector<real_type> dP2_all;
  };
  struct internal_state {
    std::vector<real_type> I_new;
    std::vector<real_type> R_new;
    std::vector<real_type> P_nV;
    std::vector<real_type> dP1;
    std::vector<real_type> dP2;
    std::vector<real_type> inv_P_nV;
    std::vector<real_type> P;
    std::vector<real_type> inv_P;
    std::vector<real_type> vacc_rate;
    std::vector<real_type> E_new;
  };
  using data_type = dust2::no_data;
  static dust2::packing packing_state(const shared_state& shared) {
    return shared.odin.packing.state;
  }
  static shared_state build_shared(cpp11::list parameters) {
    shared_state::dim_type dim;
    const real_type time_inc = dust2::r::read_real(parameters, "time_inc");
    const real_type t_infectious = dust2::r::read_real(parameters, "t_infectious");
    const real_type FOI_spillover = dust2::r::read_real(parameters, "FOI_spillover");
    const real_type R0 = dust2::r::read_real(parameters, "R0");
    const int N_age = dust2::r::read_int(parameters, "N_age");
    const real_type vaccine_efficacy = dust2::r::read_real(parameters, "vaccine_efficacy");
    const real_type year0 = dust2::r::read_real(parameters, "year0");
    const int n_years = dust2::r::read_int(parameters, "n_years");
    const real_type Pmin = static_cast<real_type>(1e-99);
    const real_type FOI_max = 1;
    const int np_E_delay = dust2::r::read_int(parameters, "np_E_delay");
    const int np_I_delay = dust2::r::read_int(parameters, "np_I_delay");
    const real_type di1 = np_E_delay - N_age;
    const real_type di2 = np_I_delay - N_age;
    const real_type beta = (R0 * time_inc) / t_infectious;
    dim.S.set({static_cast<size_t>(N_age)});
    dim.E.set({static_cast<size_t>(N_age)});
    dim.E_delay.set({static_cast<size_t>(np_E_delay)});
    dim.I.set({static_cast<size_t>(N_age)});
    dim.I_delay.set({static_cast<size_t>(np_I_delay)});
    dim.R.set({static_cast<size_t>(N_age)});
    dim.V.set({static_cast<size_t>(N_age)});
    dim.C.set({static_cast<size_t>(N_age)});
    dim.dP1.set({static_cast<size_t>(N_age)});
    dim.dP2.set({static_cast<size_t>(N_age)});
    dim.E_new.set({static_cast<size_t>(N_age)});
    dim.I_new.set({static_cast<size_t>(N_age)});
    dim.R_new.set({static_cast<size_t>(N_age)});
    dim.P_nV.set({static_cast<size_t>(N_age)});
    dim.inv_P_nV.set({static_cast<size_t>(N_age)});
    dim.P.set({static_cast<size_t>(N_age)});
    dim.inv_P.set({static_cast<size_t>(N_age)});
    dim.vacc_rate.set({static_cast<size_t>(N_age)});
    dim.S_0.set({static_cast<size_t>(N_age)});
    dim.E_0.set({static_cast<size_t>(N_age)});
    dim.E_delay0.set({static_cast<size_t>(np_E_delay)});
    dim.I_0.set({static_cast<size_t>(N_age)});
    dim.I_delay0.set({static_cast<size_t>(np_I_delay)});
    dim.R_0.set({static_cast<size_t>(N_age)});
    dim.V_0.set({static_cast<size_t>(N_age)});
    dim.dP1_all.set({static_cast<size_t>(N_age), static_cast<size_t>(n_years)});
    dim.dP2_all.set({static_cast<size_t>(N_age), static_cast<size_t>(n_years)});
    dim.vacc_rate_daily.set({static_cast<size_t>(N_age), static_cast<size_t>(n_years)});
    std::vector<real_type> vacc_rate_daily(dim.vacc_rate_daily.size);
    dust2::r::read_real_array(parameters, dim.vacc_rate_daily, vacc_rate_daily.data(), "vacc_rate_daily", true);
    std::vector<real_type> S_0(dim.S_0.size);
    dust2::r::read_real_array(parameters, dim.S_0, S_0.data(), "S_0", true);
    std::vector<real_type> E_0(dim.E_0.size);
    dust2::r::read_real_array(parameters, dim.E_0, E_0.data(), "E_0", true);
    std::vector<real_type> E_delay0(dim.E_delay0.size);
    dust2::r::read_real_array(parameters, dim.E_delay0, E_delay0.data(), "E_delay0", true);
    std::vector<real_type> I_0(dim.I_0.size);
    dust2::r::read_real_array(parameters, dim.I_0, I_0.data(), "I_0", true);
    std::vector<real_type> I_delay0(dim.I_delay0.size);
    dust2::r::read_real_array(parameters, dim.I_delay0, I_delay0.data(), "I_delay0", true);
    std::vector<real_type> R_0(dim.R_0.size);
    dust2::r::read_real_array(parameters, dim.R_0, R_0.data(), "R_0", true);
    std::vector<real_type> V_0(dim.V_0.size);
    dust2::r::read_real_array(parameters, dim.V_0, V_0.data(), "V_0", true);
    std::vector<real_type> dP1_all(dim.dP1_all.size);
    dust2::r::read_real_array(parameters, dim.dP1_all, dP1_all.data(), "dP1_all", true);
    std::vector<real_type> dP2_all(dim.dP2_all.size);
    dust2::r::read_real_array(parameters, dim.dP2_all, dP2_all.data(), "dP2_all", true);
    shared_state::odin_internals_type odin;
    odin.packing.state = dust2::packing{
      {"day", {}},
      {"year", {}},
      {"FOI_total", {}},
      {"S", std::vector<size_t>(dim.S.dim.begin(), dim.S.dim.end())},
      {"E", std::vector<size_t>(dim.E.dim.begin(), dim.E.dim.end())},
      {"E_delay", std::vector<size_t>(dim.E_delay.dim.begin(), dim.E_delay.dim.end())},
      {"I", std::vector<size_t>(dim.I.dim.begin(), dim.I.dim.end())},
      {"I_delay", std::vector<size_t>(dim.I_delay.dim.begin(), dim.I_delay.dim.end())},
      {"R", std::vector<size_t>(dim.R.dim.begin(), dim.R.dim.end())},
      {"V", std::vector<size_t>(dim.V.dim.begin(), dim.V.dim.end())},
      {"C", std::vector<size_t>(dim.C.dim.begin(), dim.C.dim.end())}
    };
    odin.packing.state.copy_offset(odin.offset.state.begin());
    return shared_state{odin, dim, time_inc, t_infectious, FOI_spillover, R0, N_age, vaccine_efficacy, year0, n_years, Pmin, FOI_max, np_E_delay, np_I_delay, di1, di2, beta, vacc_rate_daily, S_0, E_0, E_delay0, I_0, I_delay0, R_0, V_0, dP1_all, dP2_all};
  }
  static internal_state build_internal(const shared_state& shared) {
    std::vector<real_type> I_new(shared.dim.I_new.size);
    std::vector<real_type> R_new(shared.dim.R_new.size);
    std::vector<real_type> P_nV(shared.dim.P_nV.size);
    std::vector<real_type> dP1(shared.dim.dP1.size);
    std::vector<real_type> dP2(shared.dim.dP2.size);
    std::vector<real_type> inv_P_nV(shared.dim.inv_P_nV.size);
    std::vector<real_type> P(shared.dim.P.size);
    std::vector<real_type> inv_P(shared.dim.inv_P.size);
    std::vector<real_type> vacc_rate(shared.dim.vacc_rate.size);
    std::vector<real_type> E_new(shared.dim.E_new.size);
    return internal_state{I_new, R_new, P_nV, dP1, dP2, inv_P_nV, P, inv_P, vacc_rate, E_new};
  }
  static void update_shared(cpp11::list parameters, shared_state& shared) {
    shared.time_inc = dust2::r::read_real(parameters, "time_inc", shared.time_inc);
    shared.t_infectious = dust2::r::read_real(parameters, "t_infectious", shared.t_infectious);
    shared.FOI_spillover = dust2::r::read_real(parameters, "FOI_spillover", shared.FOI_spillover);
    shared.R0 = dust2::r::read_real(parameters, "R0", shared.R0);
    shared.vaccine_efficacy = dust2::r::read_real(parameters, "vaccine_efficacy", shared.vaccine_efficacy);
    shared.year0 = dust2::r::read_real(parameters, "year0", shared.year0);
    shared.beta = (shared.R0 * shared.time_inc) / shared.t_infectious;
    dust2::r::read_real_array(parameters, shared.dim.vacc_rate_daily, shared.vacc_rate_daily.data(), "vacc_rate_daily", false);
    dust2::r::read_real_array(parameters, shared.dim.S_0, shared.S_0.data(), "S_0", false);
    dust2::r::read_real_array(parameters, shared.dim.E_0, shared.E_0.data(), "E_0", false);
    dust2::r::read_real_array(parameters, shared.dim.E_delay0, shared.E_delay0.data(), "E_delay0", false);
    dust2::r::read_real_array(parameters, shared.dim.I_0, shared.I_0.data(), "I_0", false);
    dust2::r::read_real_array(parameters, shared.dim.I_delay0, shared.I_delay0.data(), "I_delay0", false);
    dust2::r::read_real_array(parameters, shared.dim.R_0, shared.R_0.data(), "R_0", false);
    dust2::r::read_real_array(parameters, shared.dim.V_0, shared.V_0.data(), "V_0", false);
    dust2::r::read_real_array(parameters, shared.dim.dP1_all, shared.dP1_all.data(), "dP1_all", false);
    dust2::r::read_real_array(parameters, shared.dim.dP2_all, shared.dP2_all.data(), "dP2_all", false);
  }
  static void update_internal(const shared_state& shared, internal_state& internal) {
  }
  static void initial(real_type time, const shared_state& shared, internal_state& internal, rng_state_type& rng_state, real_type* state) {
    state[0] = shared.time_inc;
    state[1] = shared.year0 - 1;
    state[2] = shared.FOI_spillover;
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      state[i - 1 + 3] = shared.S_0[i - 1];
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      state[i - 1 + shared.odin.offset.state[4]] = shared.E_0[i - 1];
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.np_E_delay); ++i) {
      state[i - 1 + shared.odin.offset.state[5]] = shared.E_delay0[i - 1];
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      state[i - 1 + shared.odin.offset.state[6]] = shared.I_0[i - 1];
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.np_I_delay); ++i) {
      state[i - 1 + shared.odin.offset.state[7]] = shared.I_delay0[i - 1];
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      state[i - 1 + shared.odin.offset.state[8]] = shared.R_0[i - 1];
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      state[i - 1 + shared.odin.offset.state[9]] = shared.V_0[i - 1];
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      state[i - 1 + shared.odin.offset.state[10]] = 0;
    }
  }
  static void update(real_type time, real_type dt, const real_type* state, const shared_state& shared, internal_state& internal, rng_state_type& rng_state, real_type* state_next) {
    const auto day = state[0];
    const auto * S = state + 3;
    const auto * E = state + shared.odin.offset.state[4];
    const auto * E_delay = state + shared.odin.offset.state[5];
    const auto * I = state + shared.odin.offset.state[6];
    const auto * I_delay = state + shared.odin.offset.state[7];
    const auto * R = state + shared.odin.offset.state[8];
    const auto * V = state + shared.odin.offset.state[9];
    const real_type year_i = monty::math::floor(day / 365) + 1;
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      internal.I_new[i - 1] = E_delay[static_cast<int>(i + shared.di1) - 1];
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      internal.R_new[i - 1] = I_delay[static_cast<int>(i + shared.di2) - 1];
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      internal.P_nV[i - 1] = S[i - 1] + R[i - 1];
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      internal.dP1[i - 1] = shared.dP1_all[i - 1 + (year_i - 1) * shared.dim.dP1_all.mult[1]] * shared.time_inc;
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      internal.dP2[i - 1] = shared.dP2_all[i - 1 + (year_i - 1) * shared.dim.dP2_all.mult[1]] * shared.time_inc;
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      internal.inv_P_nV[i - 1] = 1 / internal.P_nV[i - 1];
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      internal.P[i - 1] = internal.P_nV[i - 1] + V[i - 1];
    }
    const real_type P_tot = dust2::array::sum<real_type>(internal.P.data(), shared.dim.P);
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      internal.inv_P[i - 1] = 1 / internal.P[i - 1];
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      internal.vacc_rate[i - 1] = shared.vacc_rate_daily[i - 1 + (year_i - 1) * shared.dim.vacc_rate_daily.mult[1]] * shared.vaccine_efficacy * shared.time_inc * internal.P[i - 1];
    }
    const real_type FOI_sum = monty::math::min(shared.FOI_max, shared.beta * (dust2::array::sum<real_type>(I, shared.dim.I) / P_tot) + (shared.FOI_spillover * shared.time_inc));
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      internal.E_new[i - 1] = monty::random::binomial<real_type>(rng_state, static_cast<int>(S[i - 1]), FOI_sum);
    }
    state_next[0] = day + shared.time_inc;
    state_next[1] = year_i + shared.year0 - 1;
    state_next[2] = FOI_sum;
    state_next[3] = monty::math::max(shared.Pmin, S[0] - internal.E_new[0] - internal.vacc_rate[0] * S[0] * internal.inv_P_nV[0] + internal.dP1[0] - (internal.dP2[0] * S[0] * internal.inv_P[0]));
    for (size_t i = 2; i <= static_cast<size_t>(shared.N_age); ++i) {
      state_next[i - 1 + 3] = monty::math::max(shared.Pmin, S[i - 1] - internal.E_new[i - 1] - internal.vacc_rate[i - 1] * S[i - 1] * internal.inv_P_nV[i - 1] + (internal.dP1[i - 1] * S[i - 1 - 1] * internal.inv_P[i - 1 - 1]) - (internal.dP2[i - 1] * S[i - 1] * internal.inv_P[i - 1]));
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      state_next[i - 1 + shared.odin.offset.state[4]] = monty::math::max(shared.Pmin, E[i - 1] + internal.E_new[i - 1] - internal.I_new[i - 1]);
    }
    for (size_t i = (shared.N_age + 1); i <= static_cast<size_t>(shared.np_E_delay); ++i) {
      state_next[i - 1 + shared.odin.offset.state[5]] = E_delay[static_cast<int>(i) - shared.N_age - 1];
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      state_next[i - 1 + shared.odin.offset.state[5]] = internal.E_new[i - 1];
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      state_next[i - 1 + shared.odin.offset.state[6]] = monty::math::max(shared.Pmin, I[i - 1] + internal.I_new[i - 1] - internal.R_new[i - 1]);
    }
    for (size_t i = (shared.N_age + 1); i <= static_cast<size_t>(shared.np_I_delay); ++i) {
      state_next[i - 1 + shared.odin.offset.state[7]] = I_delay[static_cast<int>(i) - shared.N_age - 1];
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      state_next[i - 1 + shared.odin.offset.state[7]] = internal.I_new[i - 1];
    }
    state_next[shared.odin.offset.state[8]] = monty::math::max(shared.Pmin, R[0] + internal.R_new[0] - internal.vacc_rate[0] * R[0] * internal.inv_P_nV[0] - (internal.dP2[0] * R[0] * internal.inv_P[0]));
    for (size_t i = 2; i <= static_cast<size_t>(shared.N_age); ++i) {
      state_next[i - 1 + shared.odin.offset.state[8]] = monty::math::max(shared.Pmin, R[i - 1] + internal.R_new[i - 1] - internal.vacc_rate[i - 1] * R[i - 1] * internal.inv_P_nV[i - 1] + (internal.dP1[i - 1] * R[i - 1 - 1] * internal.inv_P[i - 1 - 1]) - (internal.dP2[i - 1] * R[i - 1] * internal.inv_P[i - 1]));
    }
    state_next[shared.odin.offset.state[9]] = monty::math::max(shared.Pmin, V[0] + internal.vacc_rate[0] - (internal.dP2[0] * V[0] * internal.inv_P[0]));
    for (size_t i = 2; i <= static_cast<size_t>(shared.N_age); ++i) {
      state_next[i - 1 + shared.odin.offset.state[9]] = monty::math::max(shared.Pmin, V[i - 1] + internal.vacc_rate[i - 1] + (internal.dP1[i - 1] * V[i - 1 - 1] * internal.inv_P[i - 1 - 1]) - (internal.dP2[i - 1] * V[i - 1] * internal.inv_P[i - 1]));
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      state_next[i - 1 + shared.odin.offset.state[10]] = internal.I_new[i - 1];
    }
  }
};

#include <cpp11.hpp>
#include <dust2/r/discrete/system.hpp>

[[cpp11::register]]
SEXP dust2_system_SEIRVModelDelay_alloc(cpp11::list r_pars, cpp11::sexp r_time, cpp11::list r_time_control, cpp11::sexp r_n_particles, cpp11::sexp r_n_groups, cpp11::sexp r_seed, cpp11::sexp r_deterministic, cpp11::sexp r_n_threads) {
  return dust2::r::dust2_discrete_alloc<SEIRVModelDelay>(r_pars, r_time, r_time_control, r_n_particles, r_n_groups, r_seed, r_deterministic, r_n_threads);
}
[[cpp11::register]]
SEXP dust2_system_SEIRVModelDelay_run_to_time(cpp11::sexp ptr, cpp11::sexp r_time) {
  return dust2::r::dust2_system_run_to_time<dust2::dust_discrete<SEIRVModelDelay>>(ptr, r_time);
}

[[cpp11::register]]
SEXP dust2_system_SEIRVModelDelay_state(cpp11::sexp ptr, cpp11::sexp r_index_state, cpp11::sexp r_index_particle, cpp11::sexp r_index_group, bool preserve_particle_dimension, bool preserve_group_dimension) {
  return dust2::r::dust2_system_state<dust2::dust_discrete<SEIRVModelDelay>>(ptr, r_index_state, r_index_particle, r_index_group, preserve_particle_dimension, preserve_group_dimension);
}

[[cpp11::register]]
SEXP dust2_system_SEIRVModelDelay_time(cpp11::sexp ptr) {
  return dust2::r::dust2_system_time<dust2::dust_discrete<SEIRVModelDelay>>(ptr);
}

[[cpp11::register]]
SEXP dust2_system_SEIRVModelDelay_set_state_initial(cpp11::sexp ptr) {
  return dust2::r::dust2_system_set_state_initial<dust2::dust_discrete<SEIRVModelDelay>>(ptr);
}

[[cpp11::register]]
SEXP dust2_system_SEIRVModelDelay_set_state(cpp11::sexp ptr, cpp11::list r_state) {
  return dust2::r::dust2_system_set_state<dust2::dust_discrete<SEIRVModelDelay>>(ptr, r_state);
}

[[cpp11::register]]
SEXP dust2_system_SEIRVModelDelay_reorder(cpp11::sexp ptr, cpp11::integers r_index) {
  return dust2::r::dust2_system_reorder<dust2::dust_discrete<SEIRVModelDelay>>(ptr, r_index);
}

[[cpp11::register]]
SEXP dust2_system_SEIRVModelDelay_rng_state(cpp11::sexp ptr) {
  return dust2::r::dust2_system_rng_state<dust2::dust_discrete<SEIRVModelDelay>>(ptr);
}

[[cpp11::register]]
SEXP dust2_system_SEIRVModelDelay_set_rng_state(cpp11::sexp ptr, cpp11::sexp r_rng_state) {
  return dust2::r::dust2_system_set_rng_state<dust2::dust_discrete<SEIRVModelDelay>>(ptr, r_rng_state);
}

[[cpp11::register]]
SEXP dust2_system_SEIRVModelDelay_set_time(cpp11::sexp ptr, cpp11::sexp r_time) {
  return dust2::r::dust2_system_set_time<dust2::dust_discrete<SEIRVModelDelay>>(ptr, r_time);
}

[[cpp11::register]]
SEXP dust2_system_SEIRVModelDelay_update_pars(cpp11::sexp ptr, cpp11::list pars) {
  return dust2::r::dust2_system_update_pars<dust2::dust_discrete<SEIRVModelDelay>>(ptr, pars);
}

[[cpp11::register]]
SEXP dust2_system_SEIRVModelDelay_simulate(cpp11::sexp ptr, cpp11::sexp r_times, cpp11::sexp r_index_state, bool preserve_particle_dimension, bool preserve_group_dimension) {
  return dust2::r::dust2_system_simulate<dust2::dust_discrete<SEIRVModelDelay>>(ptr, r_times, r_index_state, preserve_particle_dimension, preserve_group_dimension);
}
