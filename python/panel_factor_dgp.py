# Lint as: python3
"""Generate simulated data for testing causal methods."""

"""
Copyright 2020 Google LLC

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import itertools
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
import seaborn as sns
import tqdm


class SimulateDGP(object):
  """Formulate and sample from data generating process (DGP) for testing.

  The generated data for the counterfactual follows an autoregressive process
  of order one (AR(1)) meaning the outcome of interest y(t) at time 't' is
  generated from the past history y(t-1) + some factor structure and an
  independent identically distributed (iid) noise component. The factor
  structure is described in more detail in 'generate_factors' but adds further
  time series structure to the data.

  The object of interest is the treatment effect which is applied additively to
  the counterfactual at each relevant period. Samples from this dgp are drawn
  from 'get_new_data' which adds the iid noise component to the outcomes.
  """

  def __init__(self,
               num_periods,
               num_entries,
               num_factors,
               treatment_start,
               treat_proportion=0.5,
               rho=0.2,
               rho_scale=1,
               loading_scale=1,
               intercept_scale=1,
               treatment_impact=0.1,
               treatment_decay=0.7,
               treatment_decay_scale=1,
               treatment_impact_scale=1,
               number_treatments=1,
               conditional_impact=0):
    """Generate simulated data for testing causal methods.

    Args:
      num_periods: The number of periods for the generated series as a positive
        integer.
      num_entries: The number of entries for the generated series as a positive
        integer.
      num_factors: The number of common factors to simulate as a positive
        integer greater than 2.
      treatment_start: The period at which treatment begins as a positive
        integer in 0.2*num_periods < treatment_start < 0.8*num_periods.
      treat_proportion: The proportion of the entire sample that gets exposed to
        treatment in (0.1,0.9).
      rho: Autoregressive parameter for the outcome in the treated group as a
        positive real in [0,1).
      rho_scale: The scaling associated with the autoregressive parameter for
        the control group. If control accounts are suspected to grow slower than
        treated accounts this scaling should be in positive reals [0,1].
      loading_scale: Shifts the distribution of factor loadings taking values as
        positive reals in (0,2). A value of 1 draws loadings from the same
        distribution. Values closer to 0 shift the distribution of loadings for
        the control group towards 0 and towards 1 for the treated group. Values
        closer to two have the opposite effect.
      intercept_scale: Shifts the distribution by a positive real of entry
        intercepts to be smaller or larger for the control group relative to the
        treatment group. A value of 1 will have equal scaling for both groups.
        When the scaling goes to 0 all intercepts for the control group
        degenerate to 1.
      treatment_impact: The initial treatment impact to be applied across all
        treated units at the first period of treatment in positive reals [0,1).
      treatment_decay: The exponential decay parameter to be applied to the
        treatment effect across all units in positive reals [0,1). The decay
        rate is applied after the first period of treatment.
      treatment_decay_scale: Treatment effect decay heterogeneity as a positive
        real in  [0,1]. This will allow for the decay to vary across treated
        units if set to a non-zero value.
      treatment_impact_scale: Treatment effect heterogeneity as a positive real
        in [0,1]. This will allow for the impact to vary across treated units if
        set to a non-zero value.
      number_treatments: The number of treatments to consider as a positive
        integer less than or equal to 4.
      conditional_impact: Adds additional impact to accounts if they are larger
        (>75th percentile) as a non-negative real in [0,1).
    """

    self._num_periods = num_periods
    self._num_entries = num_entries
    self._num_factors = num_factors
    self._treatment_start = treatment_start
    self._treat_proportion = treat_proportion
    self._rho = rho
    self._rho_scale = rho_scale
    self._beta_param = loading_scale
    self._intercept_scale = intercept_scale
    self._treat_impact = treatment_impact
    self._treat_impact_scale = treatment_impact_scale
    self._treat_decay = treatment_decay
    self._treat_decay_scale = treatment_decay_scale
    self._num_treatments = number_treatments
    self._conditional_impact = conditional_impact

    assert self._treatment_start < self._num_periods * 0.8, (
        'Please allow at least 20% of the time interval '
        'to be post-treatment periods.')
    assert self._treatment_start > self._num_periods * 0.2, (
        'Please allow at least 20% of the time interval to be pre-treatment '
        'periods')
    assert self._treat_proportion > 0.1, (
        'Please allow at least 10% of the sample to be treated.')
    assert self._treat_proportion < 0.9, (
        'Please allow at least 10% of the sample to be in the control group.')
    assert np.isin(self._num_treatments, np.arange(
        1, 5)), 'Please limit the number of treatments to be between [1,4]'
    self.treatment_groups = [
        'treated_%i' % i for i in range(self._num_treatments)
    ]
    self.waiting = pd.DataFrame(
        0, index=range(self._num_entries), columns=self.treatment_groups)
    self.treated = self.waiting.copy().astype(bool)
    self._num_treated = {}
    self.treatment_combos = list(
        itertools.chain.from_iterable([
            itertools.combinations(self.treatment_groups, i)
            for i in np.arange(2, self._num_treatments + 1)
        ]))
    for treat in self.treatment_groups:
      self.treated.loc[:, treat] = stats.bernoulli.rvs(
          p=self._treat_proportion, size=self._num_entries).astype(bool)
      self._num_treated[treat] = sum(self.treated[treat])
      self.waiting.loc[:,
                       treat] = self._treatment_start + stats.geom(p=0.1).rvs(
                           self._num_entries)
      # Check to make sure "treated" units are actually treated.
      if np.any(self.waiting.loc[:, treat] >= self._num_periods):
        self.treated.loc[self.waiting.loc[:, treat] >= self._num_periods,
                         treat] = 0
    # TODO(b/154617471): Allow correlation of treatment with
    # (un)observables.
    self.counterfactual_outcomes, self.factors, self.loadings, \
    self.rho, self.intercept = (self._generate_counterfactual())
    self._large_accounts = self.counterfactual_outcomes.apply(
        lambda x: x.mean() > self.counterfactual_outcomes.mean().quantile(.75))

  def get_new_draw(self):
    """Draw from the data generating process.

    Returns:
      New realization of the outcomes.
    """
    counterfactual = self.counterfactual_outcomes.add(
        stats.norm(scale=.5).rvs(size=self.counterfactual_outcomes.shape))
    treated_outcomes = counterfactual.copy()
    treatment_status = {}
    treat_impact = {}
    # For each individual treatment add the treatment impact.
    for treat in self.treatment_groups:
      treat_impact[treat] = pd.DataFrame(
          0, index=treated_outcomes.index, columns=treated_outcomes.columns)
      treatment_status[treat] = treat_impact[treat].copy()
      # For each treated entry apply treatment status indicator and impact.
      for i in np.where(self.treated[treat])[0]:
        treatment_status[treat][i] = np.append(
            np.zeros(self.waiting[treat][i]),
            np.ones(self._num_periods - self.waiting[treat][i]))
        # Impact for each individual treatment.
        treat_impact[treat].loc[self.waiting[treat][i]:,
                                i] = self._generate_treatment_impact(
                                    treatment_duration=self._num_periods -
                                    self.waiting[treat][i],
                                    impact=self._treat_impact,
                                    decay=self._treat_decay,
                                    impact_scale=self._treat_impact_scale,
                                    decay_scale=self._treat_decay_scale,
                                    conditional_impact=self._large_accounts[i] *
                                    self._conditional_impact)
        treated_outcomes.loc[self.waiting[treat][i]:, i] = (
            treated_outcomes.loc[self.waiting[treat][i]:,
                                 i].add(treat_impact[treat][i]))
    # For each combination of treatments add the joint effect.
    for treat in self.treatment_combos:
      treat_impact[treat] = pd.DataFrame(
          0, index=treated_outcomes.index, columns=treated_outcomes.columns)
      for i in np.where(self.treated.loc[:, treat].all(1))[0]:
        treated_outcomes.loc[self.waiting.loc[:, treat].loc[i].max():,
                             i].subtract(.05)
        treat_impact[treat].loc[self.waiting.loc[:, treat].loc[i].max():,
                                i] = -.05
    # Combine the information into a single data frame.
    treat_status = pd.concat(
        [treatment_status[key].stack() for key in treatment_status.keys()], 1)
    treat_status.columns = self.treatment_groups
    data = pd.concat((treat_status, np.exp(
        treated_outcomes.stack()), np.exp(counterfactual.stack())), 1)
    data.columns = ['treatperiod_%i' % i for i in range(self._num_treatments)
                   ] + ['target', 'counter-factual']
    data.index.names = ['period', 'entry']
    return data, treated_outcomes.sub(counterfactual), treat_impact

  @staticmethod
  def _draw_trunc_norm(loc, scale, a=1e-5, b=0.9):
    """Draw from a truncated normal distribution with endpoints [a,b].

    This is a modification of the scipy.stats truncated normal distribution
    where the end points are explicitly defined. The draws will be confined to
    (0,0.9] where the lower bound is taken to be '1e-5'.

    Args:
      loc: The location parameter associated with the normal distribution as a
        real number.
      scale: The scale parameter associated with the normal distribution as a
        positive real.
      a: Lower bound for the truncated normal distribution as a real number.
      b: Upper bound for the truncated normal distribution as a real number.

    Returns:
      A frozen instance of stats.truncnorm evaluated at loc and scale with
      endpoints a and b.
    """
    try:
      a, b = (a - loc) / scale, (b - loc) / scale
      return stats.truncnorm(a=a, b=b, loc=loc, scale=scale)

    except ZeroDivisionError:
      return stats.norm(loc=loc, scale=0)

  def _get_random_innovations(self):
    """Generate random period shocks for factor models.

    Returns:
      Series object of length 'num_periods' containing jumps at random
      intervals.
    """
    discrete_unif = stats.randint(low=1, high=13)
    grid = np.arange(1, 52)
    how_many_jumps = discrete_unif.rvs()
    jump_location = np.random.choice(grid, how_many_jumps, replace=False)
    jump_location.sort()
    jump_location = pd.Series(np.append(np.append(0, jump_location), 52))
    jump_length = jump_location.diff(1).dropna()
    c = 1
    innovation = {}
    for j in jump_length:
      innovation[c] = pd.Series(
          np.ones(int(j)),
          index=np.arange(jump_location[c - 1], jump_location[c - 1] + j))
      c += 1
    innovation = pd.concat(innovation, 1).fillna(0)
    innovation = innovation @ np.random.uniform(
        -1, 1, size=[len(innovation.columns), 1])
    return pd.Series(
        np.tile(innovation.values.ravel(),
                int(np.ceil(self._num_periods / 52)))[:self._num_periods],
        index=range(self._num_periods))

  def get_autoregressive_process(self,
                                 name,
                                 noise,
                                 innovation,
                                 rho=0.2,
                                 loadings=None):
    """Construct AR(1) process for simulations.

    The stochastic process for a random variable X is generated
    deterministically from the previous value scaled by 'rho' + an idiosyncratic
    shock governed by 'noise'. The initial state is purely random and we allow
    500 periods for the process to converge prior to taking samples.

    Args:
      name: Label for generated series.
      noise: Frozen random variable to generate noi
      innovation: Output of "get_random_innovations" or vector of length
        'num_periods' containing various jumps.
      rho: Parameter value to scale dependence on past realizations of the
        random variable as a scalar in [0,1).
      loadings: Parameter values to scale the innovations as a vector in all
        reals.

    Returns:
      Series object for an AR(1) process with additional
      random and deterministic innovations.
    """
    series = pd.Series(index=range(500 + self._num_periods), name=str(name))
    for t in series.index:
      if t == 0:
        series.loc[t] = noise.rvs()
      else:
        series.loc[t] = rho * series.loc[t - 1] + noise.rvs()
      if t > 499:
        if loadings is None:
          series.loc[t] += innovation[t - 500]
        else:
          series.loc[t] += (innovation.loc[t - 500] @ loadings)

    return series.loc[500:].reset_index(drop=True)

  def generate_factors(self, sigma=0.1):
    """Generates common factors from a specified data generating process.

    We generate a prespecified number of factors to resemble what occurs in real
    world data when factors are estimated directly. To this end the first factor
    is typically a time trend and the following are various seasonality
    components.
      - The first factor is always generated as t + N(0,sigma)
      - The remaining factors are AR(1) with seasonality effects, e.g.,
          f_t = fixed_effect + rho f_{t-1} + N(0,sigma)

      The second and third factors are always month and quarter effects. Other
      factors are then randomly sampled in the following procedure:
        - How many breaks ~ discrete-Unif(1,12).
        - Where are the breaks sampled without replacement from a grid [1-51].

      The fixed effects are included to generate jumps at differing intervals
      reflecting observed factor behavior in real data.

    Args:
      sigma: Standard deviation for iid noise component in generated factors as
        a positive real.

    Returns:
      Simulated factors for use in monte carlo estimation.
    """
    assert self._num_factors > 2, 'Please specify at least 3 factors.'
    factors = pd.DataFrame(
        index=range(self._num_periods),
        columns=['Factor %s' % (i + 1) for i in range(self._num_factors)])
    rand_normal = stats.norm(0, sigma)
    # Quarter and month effects for factors 2 & 3.
    quarter_effects = np.random.uniform(-1, 1, size=[4, 1])
    month_effects = np.random.uniform(-1, 1, size=[13, 1])
    quarter_effects = pd.Series(
        np.tile(
            np.kron(quarter_effects, np.ones([13, 1])).ravel(),
            int(np.ceil(self._num_periods / 52)))[:self._num_periods],
        index=range(self._num_periods))
    month_effects = pd.Series(
        np.tile(
            np.kron(month_effects, np.ones([4, 1])).ravel(),
            int(np.ceil(self._num_periods / 52)))[:self._num_periods],
        index=range(self._num_periods))
    # First factor is a pure time trend.
    factors['Factor 1'] = np.arange(
        1, self._num_periods + 1) / self._num_periods + rand_normal.rvs(
            size=self._num_periods)
    # Second factor.
    factors['Factor 2'] = self.get_autoregressive_process(
        'Factor 2', rand_normal, quarter_effects)
    # Third factor.
    factors['Factor 3'] = self.get_autoregressive_process(
        'Factor 3', rand_normal, month_effects)
    # All other factors.
    if self._num_factors > 3:
      for i in np.arange(4, self._num_factors + 1):
        factors['Factor %s' % i] = self.get_autoregressive_process(
            'Factor %s' % i, rand_normal, self._get_random_innovations())

    return factors

  # TODO(b/154616170): Determine upper bound by treatment attributes.
  def _generate_treatment_impact(self,
                                 treatment_duration,
                                 impact=0.1,
                                 decay=0.7,
                                 impact_scale=0,
                                 decay_scale=0,
                                 conditional_impact=0):
    """Generate treatment impact from a simple impact decay model.

    Args:
      treatment_duration: The number of periods an entry undergoes treatment as
        a positive integer.
      impact: The impact for the first period effect as a real number.
      decay: What proportion of the effect remains period over period as a real
        in [0,1).
      impact_scale: Positive real in (0,1] to scale the initial impact allowing
        for heterogeneity across treated units.
      decay_scale: Positive real in (0,1] to scale the decay rate allowing for
        heterogeneity across treated units.
      conditional_impact: Non-negative real in [0,1] to allow for heterogeneity
        in the treatment based on some observable, e.g., account size.

    Returns:
      The treatment impact drawn from decaying truncated normal distributions.
    """
    initial_impact = self._draw_trunc_norm(
        loc=impact, scale=impact_scale, b=0.2).rvs()
    decay_rate = self._draw_trunc_norm(loc=decay, scale=decay_scale).rvs()

    return [(conditional_impact + initial_impact) * (decay_rate**i)
            for i in range(treatment_duration)]

  def _generate_counterfactual(self):
    """Generate counterfactual outcomes for all entries.

    These outcomes are generated using an AR(1) factor model.

    Returns:
      Counterfactual outcomes for every entry.
    """
    data = pd.DataFrame(
        0, index=range(self._num_periods), columns=range(self._num_entries))
    # TODO(b/154617921): Allow for common factor structure to be violated.
    factors = self.generate_factors()
    noise = stats.norm(scale=0)
    loadings = pd.DataFrame(
        0, index=range(self._num_entries), columns=factors.columns)
    rho = pd.Series(0, index=range(self._num_entries))
    intercept = rho.copy()
    # TODO(b/154618251): Allow distribution of counterfactuals to vary
    # across treatment combinations.
    # All treated entries have the same counterfactual distribution.
    treat_status = self.treated.any(1)
    for i in tqdm.tqdm(data.columns):
      if treat_status[i]:
        intercept.loc[i] = stats.expon(0, 1).rvs()
        loadings.loc[i] = stats.beta(2 - self._beta_param,
                                     1).rvs(size=[self._num_factors])
        rho.loc[i] = np.random.uniform(0, self._rho)
        data[i] = self.get_autoregressive_process(
            'outcome %s' % i,
            noise,
            factors,
            rho=rho.loc[i],
            loadings=loadings.loc[i]) + intercept[i]
      else:
        intercept.loc[i] = stats.expon(0, self._intercept_scale).rvs()
        loadings.loc[i] = stats.beta(
            1,
            2 - self._beta_param,
        ).rvs(size=[self._num_factors])
        rho.loc[i] = np.random.uniform(0, self._rho * self._rho_scale)
        data[i] = self.get_autoregressive_process(
            'outcome %s' % i,
            noise,
            factors,
            rho=rho.loc[i],
            loadings=loadings.loc[i]) + intercept[i]
    return data, factors, loadings, rho, intercept

  def plot_factors(self):
    """Generate plots for factors.

    Returns:
      Time series plot for each factor.
    """
    plt.figure(figsize=[10, 10])
    self.factors.plot(
        kind='line', figsize=[10, 10], ls='--', lw=2, subplots=True)

  def plot_treatment_effects(self, treatment_effect):
    """Plot average and individual treatment effects.

    Args:
      treatment_effect: True treatment effects generated from simulation.

    Returns:
      Time series of individual and average treatment effects.
    """
    # Sort true impacts.
    treatment_effect = treatment_effect.replace({
        0: np.nan
    }).apply(lambda x: pd.Series(x.dropna().values))
    # Generate plots for true individual level effects and implied average.
    _, ax = plt.subplots(1, 2, sharey=True, figsize=[20, 8])
    ax[0].plot(treatment_effect.mean(1), label='Average Impact')
    ax[0].set_title('Average Impact')
    ax[0].set_xlabel('Periods Since First Treatment')
    ax[0].axhline(0, c='red')
    ax[1].plot(treatment_effect)
    ax[1].set_title('Individual Impacts')
    ax[1].axhline(0, c='red')
    ax[1].set_xlabel('Periods Since First Treatment')
    plt.suptitle('Average Treatment Effect')

  def plot_autoregressive_params(self):
    """Plot the empirical distribution of simulation AR(1) parameters.

    Returns:
      Kernel density plot of AR(1) parameters.
    """
    plt.figure(figsize=[10, 10])
    plt.title('Autoregressive Parameters (Persistence)')
    sns.distplot(
        self.rho[self.treated.any(1).ravel().astype(bool)],
        label='treated',
        color='red')
    sns.distplot(
        self.rho[~self.treated.any(1).ravel().astype(bool)],
        label='control',
        color='black')
    plt.legend()

  def plot_factor_loadings(self):
    """Plot the empirical distribution of factor loadings.

    Returns:
      Kernel density plots of factor loadings.
    """
    plt.figure(figsize=[10, 10])
    plt.suptitle('Factor Loadings')
    for i in range(4):
      plt.subplot(2, 2, i + 1)
      sns.distplot(
          self.loadings.loc[self.treated.any(1).ravel()].iloc[:, i],
          label='Treated Group',
          color='red')
      sns.distplot(
          self.loadings.loc[~self.treated.any(1).ravel()].iloc[:, i],
          label='Control Group',
          color='black')
      plt.legend()
