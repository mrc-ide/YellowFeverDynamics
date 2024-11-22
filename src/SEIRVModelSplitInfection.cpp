// Generated by dust2 (version 0.3.5) - do not edit

// Generated by odin2 (version 0.3.5) - do not edit
#include <dust2/common.hpp>
// [[dust2::class(SEIRVModelSplitInfection)]]
// [[dust2::time_type(discrete)]]
// [[dust2::parameter(time_inc, type = "real_type", rank = 0, required = TRUE, constant = FALSE)]]
// [[dust2::parameter(t_incubation, type = "real_type", rank = 0, required = TRUE, constant = FALSE)]]
// [[dust2::parameter(t_latent, type = "real_type", rank = 0, required = TRUE, constant = FALSE)]]
// [[dust2::parameter(t_infectious, type = "real_type", rank = 0, required = TRUE, constant = FALSE)]]
// [[dust2::parameter(FOI_spillover, type = "real_type", rank = 0, required = TRUE, constant = FALSE)]]
// [[dust2::parameter(R0, type = "real_type", rank = 0, required = TRUE, constant = FALSE)]]
// [[dust2::parameter(N_age, type = "int", rank = 0, required = TRUE, constant = TRUE)]]
// [[dust2::parameter(vacc_rate_daily, type = "real_type", rank = 2, required = TRUE, constant = FALSE)]]
// [[dust2::parameter(vaccine_efficacy, type = "real_type", rank = 0, required = TRUE, constant = FALSE)]]
// [[dust2::parameter(year0, type = "real_type", rank = 0, required = TRUE, constant = FALSE)]]
// [[dust2::parameter(S_0, type = "real_type", rank = 1, required = TRUE, constant = FALSE)]]
// [[dust2::parameter(E_0, type = "real_type", rank = 1, required = TRUE, constant = FALSE)]]
// [[dust2::parameter(I_0, type = "real_type", rank = 1, required = TRUE, constant = FALSE)]]
// [[dust2::parameter(R_0, type = "real_type", rank = 1, required = TRUE, constant = FALSE)]]
// [[dust2::parameter(V_0, type = "real_type", rank = 1, required = TRUE, constant = FALSE)]]
// [[dust2::parameter(dP1_all, type = "real_type", rank = 2, required = TRUE, constant = FALSE)]]
// [[dust2::parameter(dP2_all, type = "real_type", rank = 2, required = TRUE, constant = FALSE)]]
// [[dust2::parameter(n_years, type = "int", rank = 0, required = TRUE, constant = TRUE)]]
class SEIRVModelSplitInfection {
public:
  SEIRVModelSplitInfection() = delete;
  using real_type = double;
  using rng_state_type = monty::random::generator<real_type>;
  struct shared_state {
    struct dim_type {
      dust2::array::dimensions<1> S;
      dust2::array::dimensions<1> E_sylv;
      dust2::array::dimensions<1> E_urb;
      dust2::array::dimensions<1> I_sylv;
      dust2::array::dimensions<1> I_urb;
      dust2::array::dimensions<1> R;
      dust2::array::dimensions<1> V;
      dust2::array::dimensions<1> C_sylv;
      dust2::array::dimensions<1> C_urb;
      dust2::array::dimensions<1> dP1;
      dust2::array::dimensions<1> dP2;
      dust2::array::dimensions<1> E_new_sylv;
      dust2::array::dimensions<1> E_new_urb;
      dust2::array::dimensions<1> I_new_sylv;
      dust2::array::dimensions<1> I_new_urb;
      dust2::array::dimensions<1> R_new_sylv;
      dust2::array::dimensions<1> R_new_urb;
      dust2::array::dimensions<1> P_nV;
      dust2::array::dimensions<1> inv_P_nV;
      dust2::array::dimensions<1> P;
      dust2::array::dimensions<1> inv_P;
      dust2::array::dimensions<1> vacc_rate;
      dust2::array::dimensions<1> S_0;
      dust2::array::dimensions<1> E_0;
      dust2::array::dimensions<1> I_0;
      dust2::array::dimensions<1> R_0;
      dust2::array::dimensions<1> V_0;
      dust2::array::dimensions<2> dP1_all;
      dust2::array::dimensions<2> dP2_all;
      dust2::array::dimensions<2> vacc_rate_daily;
    } dim;
    struct offset_type {
      struct {
        size_t day;
        size_t year;
        size_t FOI_sylv;
        size_t FOI_urb;
        size_t S;
        size_t E_sylv;
        size_t E_urb;
        size_t I_sylv;
        size_t I_urb;
        size_t R;
        size_t V;
        size_t C_sylv;
        size_t C_urb;
      } state;
    } offset;
    real_type time_inc;
    real_type t_incubation;
    real_type t_latent;
    real_type t_infectious;
    real_type FOI_spillover;
    real_type R0;
    int N_age;
    real_type vaccine_efficacy;
    real_type year0;
    int n_years;
    real_type Pmin;
    real_type FOI_max;
    real_type rate1;
    real_type rate2;
    real_type beta;
    real_type FOI_sylv_cur;
    std::vector<real_type> vacc_rate_daily;
    std::vector<real_type> S_0;
    std::vector<real_type> E_0;
    std::vector<real_type> I_0;
    std::vector<real_type> R_0;
    std::vector<real_type> V_0;
    std::vector<real_type> dP1_all;
    std::vector<real_type> dP2_all;
  };
  struct internal_state {
    std::vector<real_type> E_new_sylv;
    std::vector<real_type> I_new_sylv;
    std::vector<real_type> I_new_urb;
    std::vector<real_type> R_new_sylv;
    std::vector<real_type> R_new_urb;
    std::vector<real_type> P_nV;
    std::vector<real_type> dP1;
    std::vector<real_type> dP2;
    std::vector<real_type> inv_P_nV;
    std::vector<real_type> P;
    std::vector<real_type> inv_P;
    std::vector<real_type> vacc_rate;
    std::vector<real_type> E_new_urb;
  };
  using data_type = dust2::no_data;
  static dust2::packing packing_state(const shared_state& shared) {
    return dust2::packing{
      {"day", {}},
      {"year", {}},
      {"FOI_sylv", {}},
      {"FOI_urb", {}},
      {"S", std::vector<size_t>(shared.dim.S.dim.begin(), shared.dim.S.dim.end())},
      {"E_sylv", std::vector<size_t>(shared.dim.E_sylv.dim.begin(), shared.dim.E_sylv.dim.end())},
      {"E_urb", std::vector<size_t>(shared.dim.E_urb.dim.begin(), shared.dim.E_urb.dim.end())},
      {"I_sylv", std::vector<size_t>(shared.dim.I_sylv.dim.begin(), shared.dim.I_sylv.dim.end())},
      {"I_urb", std::vector<size_t>(shared.dim.I_urb.dim.begin(), shared.dim.I_urb.dim.end())},
      {"R", std::vector<size_t>(shared.dim.R.dim.begin(), shared.dim.R.dim.end())},
      {"V", std::vector<size_t>(shared.dim.V.dim.begin(), shared.dim.V.dim.end())},
      {"C_sylv", std::vector<size_t>(shared.dim.C_sylv.dim.begin(), shared.dim.C_sylv.dim.end())},
      {"C_urb", std::vector<size_t>(shared.dim.C_urb.dim.begin(), shared.dim.C_urb.dim.end())}
    };
  }
  static dust2::packing packing_gradient(const shared_state& shared) {
    return dust2::packing{
    };
  }
  static shared_state build_shared(cpp11::list parameters) {
    shared_state::dim_type dim;
    const real_type time_inc = dust2::r::read_real(parameters, "time_inc");
    const real_type t_incubation = dust2::r::read_real(parameters, "t_incubation");
    const real_type t_latent = dust2::r::read_real(parameters, "t_latent");
    const real_type t_infectious = dust2::r::read_real(parameters, "t_infectious");
    const real_type FOI_spillover = dust2::r::read_real(parameters, "FOI_spillover");
    const real_type R0 = dust2::r::read_real(parameters, "R0");
    const int N_age = dust2::r::read_int(parameters, "N_age");
    const real_type vaccine_efficacy = dust2::r::read_real(parameters, "vaccine_efficacy");
    const real_type year0 = dust2::r::read_real(parameters, "year0");
    const int n_years = dust2::r::read_int(parameters, "n_years");
    const real_type Pmin = static_cast<real_type>(1e-99);
    const real_type FOI_max = 1;
    const real_type rate1 = time_inc / (t_incubation + t_latent);
    const real_type rate2 = time_inc / t_infectious;
    const real_type beta = (R0 * time_inc) / t_infectious;
    const real_type FOI_sylv_cur = monty::math::min(FOI_max, FOI_spillover * time_inc);
    dim.S.set({static_cast<size_t>(N_age)});
    dim.E_sylv.set({static_cast<size_t>(N_age)});
    dim.E_urb.set({static_cast<size_t>(N_age)});
    dim.I_sylv.set({static_cast<size_t>(N_age)});
    dim.I_urb.set({static_cast<size_t>(N_age)});
    dim.R.set({static_cast<size_t>(N_age)});
    dim.V.set({static_cast<size_t>(N_age)});
    dim.C_sylv.set({static_cast<size_t>(N_age)});
    dim.C_urb.set({static_cast<size_t>(N_age)});
    dim.dP1.set({static_cast<size_t>(N_age)});
    dim.dP2.set({static_cast<size_t>(N_age)});
    dim.E_new_sylv.set({static_cast<size_t>(N_age)});
    dim.E_new_urb.set({static_cast<size_t>(N_age)});
    dim.I_new_sylv.set({static_cast<size_t>(N_age)});
    dim.I_new_urb.set({static_cast<size_t>(N_age)});
    dim.R_new_sylv.set({static_cast<size_t>(N_age)});
    dim.R_new_urb.set({static_cast<size_t>(N_age)});
    dim.P_nV.set({static_cast<size_t>(N_age)});
    dim.inv_P_nV.set({static_cast<size_t>(N_age)});
    dim.P.set({static_cast<size_t>(N_age)});
    dim.inv_P.set({static_cast<size_t>(N_age)});
    dim.vacc_rate.set({static_cast<size_t>(N_age)});
    dim.S_0.set({static_cast<size_t>(N_age)});
    dim.E_0.set({static_cast<size_t>(N_age)});
    dim.I_0.set({static_cast<size_t>(N_age)});
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
    std::vector<real_type> I_0(dim.I_0.size);
    dust2::r::read_real_array(parameters, dim.I_0, I_0.data(), "I_0", true);
    std::vector<real_type> R_0(dim.R_0.size);
    dust2::r::read_real_array(parameters, dim.R_0, R_0.data(), "R_0", true);
    std::vector<real_type> V_0(dim.V_0.size);
    dust2::r::read_real_array(parameters, dim.V_0, V_0.data(), "V_0", true);
    std::vector<real_type> dP1_all(dim.dP1_all.size);
    dust2::r::read_real_array(parameters, dim.dP1_all, dP1_all.data(), "dP1_all", true);
    std::vector<real_type> dP2_all(dim.dP2_all.size);
    dust2::r::read_real_array(parameters, dim.dP2_all, dP2_all.data(), "dP2_all", true);
    shared_state::offset_type offset;
    offset.state.day = 0;
    offset.state.year = 1;
    offset.state.FOI_sylv = 2;
    offset.state.FOI_urb = 3;
    offset.state.S = 4;
    offset.state.E_sylv = 4 + dim.S.size;
    offset.state.E_urb = offset.state.E_sylv + dim.E_sylv.size;
    offset.state.I_sylv = offset.state.E_urb + dim.E_urb.size;
    offset.state.I_urb = offset.state.I_sylv + dim.I_sylv.size;
    offset.state.R = offset.state.I_urb + dim.I_urb.size;
    offset.state.V = offset.state.R + dim.R.size;
    offset.state.C_sylv = offset.state.V + dim.V.size;
    offset.state.C_urb = offset.state.C_sylv + dim.C_sylv.size;
    return shared_state{dim, offset, time_inc, t_incubation, t_latent, t_infectious, FOI_spillover, R0, N_age, vaccine_efficacy, year0, n_years, Pmin, FOI_max, rate1, rate2, beta, FOI_sylv_cur, vacc_rate_daily, S_0, E_0, I_0, R_0, V_0, dP1_all, dP2_all};
  }
  static internal_state build_internal(const shared_state& shared) {
    std::vector<real_type> E_new_sylv(shared.dim.E_new_sylv.size);
    std::vector<real_type> I_new_sylv(shared.dim.I_new_sylv.size);
    std::vector<real_type> I_new_urb(shared.dim.I_new_urb.size);
    std::vector<real_type> R_new_sylv(shared.dim.R_new_sylv.size);
    std::vector<real_type> R_new_urb(shared.dim.R_new_urb.size);
    std::vector<real_type> P_nV(shared.dim.P_nV.size);
    std::vector<real_type> dP1(shared.dim.dP1.size);
    std::vector<real_type> dP2(shared.dim.dP2.size);
    std::vector<real_type> inv_P_nV(shared.dim.inv_P_nV.size);
    std::vector<real_type> P(shared.dim.P.size);
    std::vector<real_type> inv_P(shared.dim.inv_P.size);
    std::vector<real_type> vacc_rate(shared.dim.vacc_rate.size);
    std::vector<real_type> E_new_urb(shared.dim.E_new_urb.size);
    return internal_state{E_new_sylv, I_new_sylv, I_new_urb, R_new_sylv, R_new_urb, P_nV, dP1, dP2, inv_P_nV, P, inv_P, vacc_rate, E_new_urb};
  }
  static void update_shared(cpp11::list parameters, shared_state& shared) {
    shared.time_inc = dust2::r::read_real(parameters, "time_inc", shared.time_inc);
    shared.t_incubation = dust2::r::read_real(parameters, "t_incubation", shared.t_incubation);
    shared.t_latent = dust2::r::read_real(parameters, "t_latent", shared.t_latent);
    shared.t_infectious = dust2::r::read_real(parameters, "t_infectious", shared.t_infectious);
    shared.FOI_spillover = dust2::r::read_real(parameters, "FOI_spillover", shared.FOI_spillover);
    shared.R0 = dust2::r::read_real(parameters, "R0", shared.R0);
    shared.vaccine_efficacy = dust2::r::read_real(parameters, "vaccine_efficacy", shared.vaccine_efficacy);
    shared.year0 = dust2::r::read_real(parameters, "year0", shared.year0);
    shared.rate1 = shared.time_inc / (shared.t_incubation + shared.t_latent);
    shared.rate2 = shared.time_inc / shared.t_infectious;
    shared.beta = (shared.R0 * shared.time_inc) / shared.t_infectious;
    shared.FOI_sylv_cur = monty::math::min(shared.FOI_max, shared.FOI_spillover * shared.time_inc);
    dust2::r::read_real_array(parameters, shared.dim.vacc_rate_daily, shared.vacc_rate_daily.data(), "vacc_rate_daily", false);
    dust2::r::read_real_array(parameters, shared.dim.S_0, shared.S_0.data(), "S_0", false);
    dust2::r::read_real_array(parameters, shared.dim.E_0, shared.E_0.data(), "E_0", false);
    dust2::r::read_real_array(parameters, shared.dim.I_0, shared.I_0.data(), "I_0", false);
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
    state[3] = 0;
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      state[i - 1 + 4] = shared.S_0[i - 1];
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      state[i - 1 + shared.offset.state.E_sylv] = shared.E_0[i - 1];
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      state[i - 1 + shared.offset.state.E_urb] = 0;
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      state[i - 1 + shared.offset.state.I_sylv] = shared.I_0[i - 1];
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      state[i - 1 + shared.offset.state.I_urb] = 0;
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      state[i - 1 + shared.offset.state.R] = shared.R_0[i - 1];
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      state[i - 1 + shared.offset.state.V] = shared.V_0[i - 1];
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      state[i - 1 + shared.offset.state.C_sylv] = 0;
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      state[i - 1 + shared.offset.state.C_urb] = 0;
    }
  }
  static void update(real_type time, real_type dt, const real_type* state, const shared_state& shared, internal_state& internal, rng_state_type& rng_state, real_type* state_next) {
    const auto day = state[0];
    const auto * S = state + 4;
    const auto * E_sylv = state + shared.offset.state.E_sylv;
    const auto * E_urb = state + shared.offset.state.E_urb;
    const auto * I_sylv = state + shared.offset.state.I_sylv;
    const auto * I_urb = state + shared.offset.state.I_urb;
    const auto * R = state + shared.offset.state.R;
    const auto * V = state + shared.offset.state.V;
    const real_type year_i = monty::math::floor(day / 365) + 1;
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      internal.E_new_sylv[i - 1] = monty::random::binomial<real_type>(rng_state, static_cast<int>(S[i - 1]), shared.FOI_sylv_cur);
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      internal.I_new_sylv[i - 1] = E_sylv[i - 1] * shared.rate1;
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      internal.I_new_urb[i - 1] = E_urb[i - 1] * shared.rate1;
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      internal.R_new_sylv[i - 1] = I_sylv[i - 1] * shared.rate2;
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      internal.R_new_urb[i - 1] = I_urb[i - 1] * shared.rate2;
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
    const real_type FOI_urb_cur = monty::math::min(shared.FOI_max, shared.beta * ((dust2::array::sum<real_type>(I_sylv, shared.dim.I_sylv) + dust2::array::sum<real_type>(I_urb, shared.dim.I_urb)) / P_tot));
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      internal.E_new_urb[i - 1] = monty::random::binomial<real_type>(rng_state, static_cast<int>(S[i - 1]), FOI_urb_cur);
    }
    state_next[0] = day + shared.time_inc;
    state_next[1] = year_i + shared.year0 - 1;
    state_next[2] = shared.FOI_sylv_cur;
    state_next[3] = FOI_urb_cur;
    state_next[4] = monty::math::max(shared.Pmin, S[0] - internal.E_new_sylv[0] - internal.E_new_urb[0] - internal.vacc_rate[0] * S[0] * internal.inv_P_nV[0] + internal.dP1[0] - (internal.dP2[0] * S[0] * internal.inv_P[0]));
    for (size_t i = 2; i <= static_cast<size_t>(shared.N_age); ++i) {
      state_next[i - 1 + 4] = monty::math::max(shared.Pmin, S[i - 1] - internal.E_new_sylv[i - 1] - internal.E_new_urb[i - 1] - internal.vacc_rate[i - 1] * S[i - 1] * internal.inv_P_nV[i - 1] + (internal.dP1[i - 1] * S[i - 1 - 1] * internal.inv_P[i - 1 - 1]) - (internal.dP2[i - 1] * S[i - 1] * internal.inv_P[i - 1]));
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      state_next[i - 1 + shared.offset.state.E_sylv] = monty::math::max(shared.Pmin, E_sylv[i - 1] + internal.E_new_sylv[i - 1] - internal.I_new_sylv[i - 1]);
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      state_next[i - 1 + shared.offset.state.E_urb] = monty::math::max(shared.Pmin, E_urb[i - 1] + internal.E_new_urb[i - 1] - internal.I_new_urb[i - 1]);
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      state_next[i - 1 + shared.offset.state.I_sylv] = monty::math::max(shared.Pmin, I_sylv[i - 1] + internal.I_new_sylv[i - 1] - internal.R_new_sylv[i - 1]);
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      state_next[i - 1 + shared.offset.state.I_urb] = monty::math::max(shared.Pmin, I_urb[i - 1] + internal.I_new_urb[i - 1] - internal.R_new_urb[i - 1]);
    }
    state_next[shared.offset.state.R] = monty::math::max(shared.Pmin, R[0] + internal.R_new_sylv[0] + internal.R_new_urb[0] - internal.vacc_rate[0] * R[0] * internal.inv_P_nV[0] - (internal.dP2[0] * R[0] * internal.inv_P[0]));
    for (size_t i = 2; i <= static_cast<size_t>(shared.N_age); ++i) {
      state_next[i - 1 + shared.offset.state.R] = monty::math::max(shared.Pmin, R[i - 1] + internal.R_new_sylv[i - 1] + internal.R_new_urb[i - 1] - internal.vacc_rate[i - 1] * R[i - 1] * internal.inv_P_nV[i - 1] + (internal.dP1[i - 1] * R[i - 1 - 1] * internal.inv_P[i - 1 - 1]) - (internal.dP2[i - 1] * R[i - 1] * internal.inv_P[i - 1]));
    }
    state_next[shared.offset.state.V] = monty::math::max(shared.Pmin, V[0] + internal.vacc_rate[0] - (internal.dP2[0] * V[0] * internal.inv_P[0]));
    for (size_t i = 2; i <= static_cast<size_t>(shared.N_age); ++i) {
      state_next[i - 1 + shared.offset.state.V] = monty::math::max(shared.Pmin, V[i - 1] + internal.vacc_rate[i - 1] + (internal.dP1[i - 1] * V[i - 1 - 1] * internal.inv_P[i - 1 - 1]) - (internal.dP2[i - 1] * V[i - 1] * internal.inv_P[i - 1]));
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      state_next[i - 1 + shared.offset.state.C_sylv] = internal.I_new_sylv[i - 1];
    }
    for (size_t i = 1; i <= static_cast<size_t>(shared.N_age); ++i) {
      state_next[i - 1 + shared.offset.state.C_urb] = internal.I_new_urb[i - 1];
    }
  }
  static auto zero_every(const shared_state& shared) {
    return dust2::zero_every_type<real_type>();
  }
};

#include <cpp11.hpp>
#include <dust2/r/discrete/system.hpp>

[[cpp11::register]]
SEXP dust2_system_SEIRVModelSplitInfection_alloc(cpp11::list r_pars, cpp11::sexp r_time, cpp11::list r_time_control, cpp11::sexp r_n_particles, cpp11::sexp r_n_groups, cpp11::sexp r_seed, cpp11::sexp r_deterministic, cpp11::sexp r_n_threads) {
  return dust2::r::dust2_discrete_alloc<SEIRVModelSplitInfection>(r_pars, r_time, r_time_control, r_n_particles, r_n_groups, r_seed, r_deterministic, r_n_threads);
}
[[cpp11::register]]
SEXP dust2_system_SEIRVModelSplitInfection_run_to_time(cpp11::sexp ptr, cpp11::sexp r_time) {
  return dust2::r::dust2_system_run_to_time<dust2::dust_discrete<SEIRVModelSplitInfection>>(ptr, r_time);
}

[[cpp11::register]]
SEXP dust2_system_SEIRVModelSplitInfection_state(cpp11::sexp ptr, cpp11::sexp r_index_state, cpp11::sexp r_index_particle, cpp11::sexp r_index_group, bool preserve_particle_dimension, bool preserve_group_dimension) {
  return dust2::r::dust2_system_state<dust2::dust_discrete<SEIRVModelSplitInfection>>(ptr, r_index_state, r_index_particle, r_index_group, preserve_particle_dimension, preserve_group_dimension);
}

[[cpp11::register]]
SEXP dust2_system_SEIRVModelSplitInfection_time(cpp11::sexp ptr) {
  return dust2::r::dust2_system_time<dust2::dust_discrete<SEIRVModelSplitInfection>>(ptr);
}

[[cpp11::register]]
SEXP dust2_system_SEIRVModelSplitInfection_set_state_initial(cpp11::sexp ptr) {
  return dust2::r::dust2_system_set_state_initial<dust2::dust_discrete<SEIRVModelSplitInfection>>(ptr);
}

[[cpp11::register]]
SEXP dust2_system_SEIRVModelSplitInfection_set_state(cpp11::sexp ptr, cpp11::list r_state) {
  return dust2::r::dust2_system_set_state<dust2::dust_discrete<SEIRVModelSplitInfection>>(ptr, r_state);
}

[[cpp11::register]]
SEXP dust2_system_SEIRVModelSplitInfection_reorder(cpp11::sexp ptr, cpp11::integers r_index) {
  return dust2::r::dust2_system_reorder<dust2::dust_discrete<SEIRVModelSplitInfection>>(ptr, r_index);
}

[[cpp11::register]]
SEXP dust2_system_SEIRVModelSplitInfection_rng_state(cpp11::sexp ptr) {
  return dust2::r::dust2_system_rng_state<dust2::dust_discrete<SEIRVModelSplitInfection>>(ptr);
}

[[cpp11::register]]
SEXP dust2_system_SEIRVModelSplitInfection_set_rng_state(cpp11::sexp ptr, cpp11::sexp r_rng_state) {
  return dust2::r::dust2_system_set_rng_state<dust2::dust_discrete<SEIRVModelSplitInfection>>(ptr, r_rng_state);
}

[[cpp11::register]]
SEXP dust2_system_SEIRVModelSplitInfection_set_time(cpp11::sexp ptr, cpp11::sexp r_time) {
  return dust2::r::dust2_system_set_time<dust2::dust_discrete<SEIRVModelSplitInfection>>(ptr, r_time);
}

[[cpp11::register]]
SEXP dust2_system_SEIRVModelSplitInfection_update_pars(cpp11::sexp ptr, cpp11::list pars) {
  return dust2::r::dust2_system_update_pars<dust2::dust_discrete<SEIRVModelSplitInfection>>(ptr, pars);
}

[[cpp11::register]]
SEXP dust2_system_SEIRVModelSplitInfection_simulate(cpp11::sexp ptr, cpp11::sexp r_times, cpp11::sexp r_index_state, bool preserve_particle_dimension, bool preserve_group_dimension) {
  return dust2::r::dust2_system_simulate<dust2::dust_discrete<SEIRVModelSplitInfection>>(ptr, r_times, r_index_state, preserve_particle_dimension, preserve_group_dimension);
}
