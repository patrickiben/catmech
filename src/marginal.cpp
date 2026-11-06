#include <RcppArmadillo.h>
#include <map>
#include <string>
#include <vector>

// Original blocked Gibbs implementation. Covariance convention:
// C ~ IW(df, S), density proportional to
// |C|^(-(df+d+1)/2) exp(-tr(S C^-1)/2).
// Huang--Wand: C|a ~ IW(d+1, 4 diag(1/a)),
// a_k ~ IG(1/2, 1/A^2), A=sqrt(prior_scale).
namespace {
arma::mat symmetric(const arma::mat& a) { return (a + a.t()) * 0.5; }

arma::mat inverse_spd(const arma::mat& a) {
  arma::mat out;
  if (!a.is_finite() || !arma::inv_sympd(out, symmetric(a)))
    Rcpp::stop("Marginal sampler covariance or precision is not positive definite.");
  return symmetric(out);
}

arma::mat chol_lower(const arma::mat& a) {
  arma::mat out;
  if (!a.is_finite() || !arma::chol(out, symmetric(a), "lower"))
    Rcpp::stop("Marginal sampler covariance or precision is not positive definite.");
  return out;
}

arma::vec normals(unsigned int n) {
  arma::vec out(n);
  for (unsigned int j = 0; j < n; ++j) out[j] = R::rnorm(0, 1);
  return out;
}

arma::mat inverse_wishart(double df, const arma::mat& scale) {
  unsigned int d = scale.n_rows;
  if (df <= d - 1) Rcpp::stop("Invalid inverse-Wishart degrees of freedom.");
  arma::mat a(d, d, arma::fill::zeros);
  for (unsigned int j = 0; j < d; ++j) {
    a(j, j) = std::sqrt(R::rchisq(df - j));
    for (unsigned int i = j + 1; i < d; ++i) a(i, j) = R::rnorm(0, 1);
  }
  arma::mat b = arma::solve(arma::trimatl(a), chol_lower(scale).t());
  arma::mat out = symmetric(b.t() * b);
  if (!out.is_finite()) Rcpp::stop("Nonfinite inverse-Wishart draw.");
  return out;
}

arma::vec auxiliary_draw(const arma::mat& covariance, double scale) {
  unsigned int d = covariance.n_rows;
  arma::vec precision = inverse_spd(covariance).diag();
  arma::vec out(d);
  for (unsigned int k = 0; k < d; ++k)
    out[k] = 1 / R::rgamma((2.0 + d) / 2, 1 / (2 * precision[k] + 1 / scale));
  if (!out.is_finite() || arma::any(out <= 0))
    Rcpp::stop("Nonfinite Huang-Wand auxiliary draw.");
  return out;
}

struct ConditionalCovariance {
  arma::mat predictor;
  arma::mat regression;  // k predictors by two measurement responses
  arma::mat residual;
};

ConditionalCovariance conditional_covariance(const arma::mat& covariance) {
  const unsigned int d = covariance.n_rows, k = d - 2;
  ConditionalCovariance out;
  out.residual = covariance.submat(0, 0, 1, 1);
  if (k) {
    out.predictor = covariance.submat(2, 2, d - 1, d - 1);
    out.regression = inverse_spd(out.predictor) * covariance.submat(2, 0, d - 1, 1);
    out.residual -= covariance.submat(0, 2, 1, d - 1) * out.regression;
    out.residual = symmetric(out.residual);
  } else {
    out.predictor.set_size(0, 0);
    out.regression.set_size(0, 2);
  }
  return out;
}

// Exact conjugate block conditional for the proper prior:
// S=Q_XX and Omega=Q_RR-Q_RX S^-1 Q_XR have independent HW priors;
// B|Omega ~ MN(0, regression_scale^2 I_k, Omega). R=(T,W).
arma::mat factorized_covariance_draw(const arma::mat& values, const arma::vec& auxiliary,
                                     double regression_scale) {
  const unsigned int n = values.n_rows, d = values.n_cols, k = d - 2;
  arma::mat response = values.cols(0, 1);
  if (!k) {
    arma::mat scale = response.t() * response;
    scale.diag() += 4 / auxiliary;
    return inverse_wishart(3.0 + n, scale);
  }
  const double precision = 1 / (regression_scale * regression_scale);
  arma::mat x = values.cols(2, d - 1), cross = x.t() * x;
  arma::mat predictor_scale = cross;
  predictor_scale.diag() += 4 / auxiliary.tail(k);
  arma::mat s = inverse_wishart(k + 1.0 + n, predictor_scale);
  arma::mat coefficient_precision = cross;
  coefficient_precision.diag() += precision;
  arma::mat upper = chol_lower(coefficient_precision).t();
  arma::mat mean = arma::solve(arma::trimatu(upper),
    arma::solve(arma::trimatl(upper.t()), x.t() * response));
  arma::mat residual = response - x * mean;
  // Stable form of R'R - M'(X'X + precision I)M.
  arma::mat conditional_scale = residual.t() * residual + precision * mean.t() * mean;
  conditional_scale.diag() += 4 / auxiliary.head(2);
  arma::mat omega = inverse_wishart(3.0 + n, conditional_scale);
  arma::mat z(k, 2);
  for (unsigned int j = 0; j < 2; ++j) z.col(j) = normals(k);
  arma::mat b = mean + arma::solve(arma::trimatu(upper), z) * chol_lower(omega).t();
  arma::mat out(d, d);
  out.submat(2, 2, d - 1, d - 1) = s;
  out.submat(2, 0, d - 1, 1) = s * b;
  out.submat(0, 2, 1, d - 1) = b.t() * s;
  out.submat(0, 0, 1, 1) = omega + b.t() * s * b;
  return symmetric(out);
}

arma::vec factorized_auxiliary_draw(const arma::mat& covariance, double scale) {
  ConditionalCovariance block = conditional_covariance(covariance);
  arma::vec out(covariance.n_rows);
  out.head(2) = auxiliary_draw(block.residual, scale);
  if (covariance.n_rows > 2) out.tail(covariance.n_rows - 2) = auxiliary_draw(block.predictor, scale);
  return out;
}

double log_sum_exp(double a, double b) {
  double larger = std::max(a, b);
  return larger + std::log1p(std::exp(std::min(a, b) - larger));
}

// Reversible scale-group move, integrating the Huang--Wand auxiliaries.
// Scaling one U column by c and the corresponding Psi row/column by c gives
// Jacobian c^(G+d+1). The conditional Gaussian U density contributes c^-G.
// The remaining covariance-coordinate Jacobian is c^(d+1).
double scale_log_ratio(const arma::mat& ybar, const arma::vec& sizes,
                       const arma::vec& mu, const arma::mat& u,
                       const arma::mat& si, const arma::mat& pi,
                       const arma::mat& psi,
                       unsigned int j, double eta, double scale,
                       const std::string& family, double regression_scale) {
  const double d = u.n_cols;
  arma::mat residual = ybar - u;
  residual.each_row() -= mu.t();
  arma::vec delta = std::expm1(eta) * u.col(j);
  double likelihood = arma::dot(sizes, delta % (residual * si.col(j)) -
    0.5 * arma::square(delta) * si(j, j));
  double prior_and_jacobian;
  if (family == "conditional_regression") {
    const double k = d - 2;
    if (j < 2) {
      // Omega transforms with the measurement coordinate; B's matching
      // column transforms too. Its normal density cancels that B Jacobian.
      double log_precision = std::log(2 * pi(j, j));
      double log_rate = -std::log(scale);
      double change = log_sum_exp(log_precision - 2 * eta, log_rate) -
        log_sum_exp(log_precision, log_rate);
      prior_and_jacobian = -3 * eta - 2 * change;
    } else {
      ConditionalCovariance block = conditional_covariance(psi);
      arma::mat sinv = inverse_spd(block.predictor);
      arma::rowvec coefficient = block.regression.row(j - 2);
      double quadratic = arma::as_scalar(coefficient * pi.submat(0, 0, 1, 1) * coefficient.t());
      double log_precision = std::log(2 * sinv(j - 2, j - 2));
      double log_rate = -std::log(scale);
      double change = log_sum_exp(log_precision - 2 * eta, log_rate) -
        log_sum_exp(log_precision, log_rate);
      // The Q->(S,B,Omega) prior density contributes -2 log|S|.
      prior_and_jacobian = -(2 + 2 * k) * eta - (2 + k) / 2 * change +
        (k - 1) * eta - 0.5 * quadratic / (regression_scale * regression_scale) * std::expm1(-2 * eta);
    }
  } else if (family == "huang_wand") {
    double log_precision = std::log(2 * pi(j, j));
    double log_rate = -std::log(scale);
    double log_change = log_sum_exp(log_precision - 2 * eta, log_rate) -
      log_sum_exp(log_precision, log_rate);
    prior_and_jacobian = -(2 + 2 * d) * eta - (2 + d) / 2 * log_change + (d + 1) * eta;
  } else {
    const double df = d + 2;
    prior_and_jacobian = -df * eta - .25 * scale * pi(j, j) * std::expm1(-2 * eta);
  }
  double out = likelihood + prior_and_jacobian;
  return R_finite(out) ? out : R_NegInf;
}

void scale_interweave(arma::mat& u, arma::mat& psi, const arma::mat& ybar,
                      const arma::vec& sizes, const arma::vec& mu,
                      const arma::mat& sigma, double scale, const std::string& family,
                      double regression_scale,
                      arma::vec& proposed, arma::vec& accepted) {
  arma::mat si = inverse_spd(sigma), pi = inverse_spd(psi);
  for (unsigned int j = 0; j < u.n_cols; ++j) {
    // Fixed engineering settings, never adapted to interval performance.
    for (unsigned int step = 0; step < 2; ++step) {
      const double eta = R::rnorm(0, 1), c = std::exp(eta);
      const double ratio = scale_log_ratio(ybar, sizes, mu, u, si, pi, psi,
        j, eta, scale, family, regression_scale);
      proposed[j] += 1;
      if (std::log(R::runif(0, 1)) < ratio) {
        u.col(j) *= c;
        psi.row(j) *= c; psi.col(j) *= c;
        pi.row(j) /= c; pi.col(j) /= c;
        accepted[j] += 1;
      }
    }
  }
}

struct Pattern {
  arma::uvec observed;
  arma::uvec missing;
  std::vector<unsigned int> rows;
  arma::mat observed_values;
};

void save_matrix(Rcpp::NumericVector& output, const arma::mat& value,
                 unsigned int index) {
  std::copy(value.begin(), value.end(), output.begin() + index * value.n_elem);
}

std::vector<Pattern> build_prediction_patterns(const arma::mat& y) {
  const unsigned int n = y.n_rows, d = y.n_cols;
  std::map<std::string, unsigned int> index;
  std::vector<Pattern> patterns;
  for (unsigned int i = 0; i < n; ++i) {
    std::string key;
    std::vector<unsigned int> obs;
    for (unsigned int k = 0; k < d; ++k) {
      bool present = R_finite(y(i, k));
      key += present ? '1' : '0';
      if (present) obs.push_back(k);
    }
    if (index.count(key) == 0) {
      Pattern p;
      p.observed = arma::conv_to<arma::uvec>::from(obs);
      patterns.push_back(p);
      index[key] = patterns.size() - 1;
    }
    patterns[index[key]].rows.push_back(i);
  }
  for (Pattern& pattern : patterns) {
    if (pattern.observed.n_elem) {
      pattern.observed_values.set_size(pattern.observed.n_elem, pattern.rows.size());
      for (unsigned int r = 0; r < pattern.rows.size(); ++r) {
        arma::vec response = y.row(pattern.rows[r]).t();
        pattern.observed_values.col(r) = response.elem(pattern.observed);
      }
    }
  }
  return patterns;
}

arma::vec predict_target_patterns(const arma::mat& y, const arma::uvec& group,
                                  unsigned int g, const std::vector<Pattern>& patterns,
                                  const arma::vec& mu, const arma::mat& sigma,
                                  const arma::mat& psi, bool draw) {
  const unsigned int d = y.n_cols;
  arma::mat pi = inverse_spd(psi);
  std::vector<arma::mat> precision(g, pi), observed_inverse;
  arma::mat rhs(g, d, arma::fill::zeros);
  for (const Pattern& pattern : patterns) {
    arma::mat oi;
    if (pattern.observed.n_elem) {
      oi = inverse_spd(sigma.submat(pattern.observed, pattern.observed));
      arma::mat embedded(d, d, arma::fill::zeros);
      embedded.submat(pattern.observed, pattern.observed) = oi;
      arma::vec mu_obs = mu.elem(pattern.observed);
      const unsigned int m_obs = pattern.observed.n_elem;
      const unsigned int* obs_idx = pattern.observed.memptr();
      arma::vec diff(m_obs), contribution(m_obs);
      for (unsigned int r = 0; r < pattern.rows.size(); ++r) {
        unsigned int row = pattern.rows[r];
        unsigned int j = group[row];
        precision[j] += embedded;
        const double* obs_col = pattern.observed_values.colptr(r);
        for (unsigned int k = 0; k < m_obs; ++k)
          diff[k] = obs_col[k] - mu_obs[k];
        contribution = oi * diff;
        for (unsigned int k = 0; k < m_obs; ++k)
          rhs(j, obs_idx[k]) += contribution[k];
      }
    }
    observed_inverse.push_back(oi);
  }
  arma::mat u(g, d);
  std::vector<arma::vec> site_location(g);
  for (unsigned int j = 0; j < g; ++j) {
    arma::mat upper = chol_lower(precision[j]).t();
    arma::vec value = arma::solve(arma::trimatu(upper),
      arma::solve(arma::trimatl(upper.t()), rhs.row(j).t()));
    if (draw) value += arma::solve(arma::trimatu(upper), normals(d));
    u.row(j) = value.t();
    site_location[j] = mu + value;
  }
  arma::vec out = y.col(0);
  for (unsigned int p = 0; p < patterns.size(); ++p) {
    const Pattern& pattern = patterns[p];
    if (pattern.observed.n_elem && pattern.observed[0] == 0) continue;
    arma::rowvec regression;
    double variance = sigma(0, 0);
    if (pattern.observed.n_elem) {
      arma::uvec target(1); target[0] = 0;
      regression = sigma.submat(target, pattern.observed) * observed_inverse[p];
      variance -= arma::as_scalar(regression * sigma.submat(pattern.observed, target));
    }
    if (!R_finite(variance) || variance <= 0)
      Rcpp::stop("Marginal target conditional variance is not positive.");
    const unsigned int m_obs = pattern.observed.n_elem;
    const unsigned int* obs_idx = pattern.observed.memptr();
    arma::vec diff(m_obs);
    for (unsigned int r = 0; r < pattern.rows.size(); ++r) {
      unsigned int row = pattern.rows[r];
      const arma::vec& location = site_location[group[row]];
      double value = location[0];
      if (m_obs) {
        const double* obs_col = pattern.observed_values.colptr(r);
        for (unsigned int k = 0; k < m_obs; ++k)
          diff[k] = obs_col[k] - location[obs_idx[k]];
        value += arma::as_scalar(regression * diff);
      }
      out[row] = draw ? R::rnorm(value, std::sqrt(variance)) : value;
    }
  }
  return out;
}

arma::vec predict_target(const arma::mat& y, const arma::uvec& group,
                         const arma::vec& mu, const arma::mat& sigma,
                         const arma::mat& psi, bool draw) {
  return predict_target_patterns(y, group, group.max() + 1, build_prediction_patterns(y),
    mu, sigma, psi, draw);
}
}

RcppExport SEXP catmech_marginal_gibbs_cpp(SEXP input_sexp, SEXP iterations_sexp,
                                         SEXP start_sexp, SEXP prior_sexp,
                                         SEXP fixed_sexp) {
  BEGIN_RCPP
  Rcpp::RNGScope scope;
  Rcpp::List input(input_sexp), prior(prior_sexp);
  arma::mat observed_y = Rcpp::as<arma::mat>(input["y"]);
  Rcpp::IntegerVector subject(input["subj"]);
  const unsigned int n = observed_y.n_rows, d = observed_y.n_cols;
  const unsigned int iterations = Rcpp::as<unsigned int>(iterations_sexp);
  const double scale = Rcpp::as<double>(prior["scale"]);
  const std::string family = Rcpp::as<std::string>(prior["family"]);
  const double regression_scale = prior.containsElementNamed("regression_scale") ?
    Rcpp::as<double>(prior["regression_scale"]) : 10;
  const bool interweave = prior.containsElementNamed("interweave") ?
    Rcpp::as<bool>(prior["interweave"]) : true;
  if (n < 1 || d < 1 || iterations < 1 || subject.size() != n ||
      !R_finite(scale) || scale <= 0 ||
      !R_finite(regression_scale) || regression_scale <= 0 ||
      (family != "huang_wand" && family != "inverse_wishart" && family != "conditional_regression") ||
      (family == "conditional_regression" && d < 2))
    Rcpp::stop("Malformed marginal sampler input.");
  int largest = 0;
  for (unsigned int i = 0; i < n; ++i) {
    if (subject[i] == NA_INTEGER || subject[i] < 1)
      Rcpp::stop("Marginal sampler site identifiers must be positive integers.");
    largest = std::max(largest, subject[i]);
  }
  const unsigned int g = largest;
  if (g > n) Rcpp::stop("Marginal sampler site identifiers must be contiguous.");
  arma::uvec group(n);
  arma::vec sizes(g, arma::fill::zeros);
  for (unsigned int i = 0; i < n; ++i) {
    group[i] = subject[i] - 1;
    sizes[group[i]] += 1;
  }
  if (arma::any(sizes == 0)) Rcpp::stop("Marginal sampler site identifiers must be contiguous.");

  std::map<std::string, unsigned int> pattern_index;
  std::vector<Pattern> patterns;
  arma::mat y = observed_y;
  for (unsigned int i = 0; i < n; ++i) {
    std::string key;
    std::vector<unsigned int> obs, mis;
    for (unsigned int j = 0; j < d; ++j) {
      bool missing = R_IsNA(y(i, j)) || R_IsNaN(y(i, j));
      if (!missing && !R_finite(y(i, j))) Rcpp::stop("Nonfinite marginal observed data.");
      key += missing ? '0' : '1';
      if (missing) { mis.push_back(j); y(i, j) = 0; }
      else obs.push_back(j);
    }
    if (mis.empty()) continue;
    if (pattern_index.count(key) == 0) {
      Pattern pattern;
      pattern.observed = arma::conv_to<arma::uvec>::from(obs);
      pattern.missing = arma::conv_to<arma::uvec>::from(mis);
      patterns.push_back(pattern);
      pattern_index[key] = patterns.size() - 1;
    }
    patterns[pattern_index[key]].rows.push_back(i);
  }
  for (Pattern& pattern : patterns) {
    if (pattern.observed.n_elem) {
      pattern.observed_values.set_size(pattern.observed.n_elem, pattern.rows.size());
      for (unsigned int r = 0; r < pattern.rows.size(); ++r) {
        arma::vec value = y.row(pattern.rows[r]).t();
        pattern.observed_values.col(r) = value.elem(pattern.observed);
      }
    }
  }

  arma::vec mu(d, arma::fill::zeros), a_sigma(d, arma::fill::ones), a_psi(d, arma::fill::ones);
  arma::mat sigma(d, d, arma::fill::eye), psi(d, d, arma::fill::eye);
  arma::mat u(g, d, arma::fill::zeros);
  arma::vec scale_proposed(d, arma::fill::zeros), scale_accepted(d, arma::fill::zeros);
  if (!Rf_isNull(start_sexp)) {
    Rcpp::List start(start_sexp);
    if (start.containsElementNamed("mu")) mu = Rcpp::as<arma::vec>(start["mu"]);
    if (start.containsElementNamed("sigma")) sigma = Rcpp::as<arma::mat>(start["sigma"]);
    if (start.containsElementNamed("psi")) psi = Rcpp::as<arma::mat>(start["psi"]);
    if (start.containsElementNamed("aux_sigma")) a_sigma = Rcpp::as<arma::vec>(start["aux_sigma"]);
    if (start.containsElementNamed("aux_psi")) a_psi = Rcpp::as<arma::vec>(start["aux_psi"]);
    if (start.containsElementNamed("u")) u = Rcpp::as<arma::mat>(start["u"]);
    if (start.containsElementNamed("y")) {
      arma::mat initial = Rcpp::as<arma::mat>(start["y"]);
      if (initial.n_rows != n || initial.n_cols != d || !initial.is_finite())
        Rcpp::stop("Malformed marginal completed-data start.");
      for (const Pattern& pattern : patterns)
        for (unsigned int row : pattern.rows)
          for (unsigned int col : pattern.missing) y(row, col) = initial(row, col);
    }
  }
  bool fixed = !Rf_isNull(fixed_sexp);
  if (fixed) {
    Rcpp::List values(fixed_sexp);
    sigma = Rcpp::as<arma::mat>(values["sigma"]);
    psi = Rcpp::as<arma::mat>(values["psi"]);
  }
  if (mu.n_elem != d || sigma.n_rows != d || sigma.n_cols != d ||
      psi.n_rows != d || psi.n_cols != d || u.n_rows != g || u.n_cols != d ||
      a_sigma.n_elem != d || a_psi.n_elem != d || !mu.is_finite() || !u.is_finite() ||
      !a_sigma.is_finite() || !a_psi.is_finite() || arma::any(a_sigma <= 0) || arma::any(a_psi <= 0))
    Rcpp::stop("Malformed marginal parameter start.");
  chol_lower(sigma); chol_lower(psi);
  Rcpp::NumericMatrix means(d, iterations), aux_sigma(d, iterations), aux_psi(d, iterations);
  Rcpp::NumericVector sigma_draws(static_cast<R_xlen_t>(d) * d * iterations);
  Rcpp::NumericVector psi_draws(static_cast<R_xlen_t>(d) * d * iterations);
  sigma_draws.attr("dim") = Rcpp::IntegerVector::create(d, d, iterations);
  psi_draws.attr("dim") = Rcpp::IntegerVector::create(d, d, iterations);
  const double covariance_df = family == "huang_wand" ? d + 1.0 : d + 2.0;
  for (unsigned int iteration = 0; iteration < iterations; ++iteration) {
    if (iteration % 128 == 0) Rcpp::checkUserInterrupt();
    // Augment only missing entries. Every observed component remains immutable.
    std::vector<arma::vec> site_location(g);
    for (unsigned int j = 0; j < g; ++j) site_location[j] = mu + u.row(j).t();
    for (const Pattern& pattern : patterns) {
      arma::mat regression, conditional = sigma.submat(pattern.missing, pattern.missing);
      if (pattern.observed.n_elem) {
        regression = sigma.submat(pattern.missing, pattern.observed) *
          inverse_spd(sigma.submat(pattern.observed, pattern.observed));
        conditional -= regression * sigma.submat(pattern.observed, pattern.missing);
      }
      arma::mat chol = chol_lower(conditional);
      const unsigned int m_obs = pattern.observed.n_elem, m_mis = pattern.missing.n_elem;
      const unsigned int* obs_idx = pattern.observed.memptr();
      const unsigned int* mis_idx = pattern.missing.memptr();
      arma::vec diff(m_obs), center(m_mis), draw(m_mis);
      for (unsigned int r = 0; r < pattern.rows.size(); ++r) {
        unsigned int row = pattern.rows[r];
        const arma::vec& location = site_location[group[row]];
        for (unsigned int k = 0; k < m_mis; ++k)
          center[k] = location[mis_idx[k]];
        if (m_obs) {
          const double* obs_col = pattern.observed_values.colptr(r);
          for (unsigned int k = 0; k < m_obs; ++k)
            diff[k] = obs_col[k] - location[obs_idx[k]];
          center += regression * diff;
        }
        draw = center + chol * normals(m_mis);
        for (unsigned int k = 0; k < m_mis; ++k)
          y(row, mis_idx[k]) = draw[k];
      }
    }
    arma::mat ybar(g, d, arma::fill::zeros);
    for (unsigned int i = 0; i < n; ++i) ybar.row(group[i]) += y.row(i);
    ybar.each_col() /= sizes;
    // Draw mu after integrating site effects, then each u conditional on mu.
    // This is one block draw from p(mu,u | completed Y, Sigma,Psi).
    arma::mat si = inverse_spd(sigma), pi = inverse_spd(psi);
    arma::mat precision(d, d, arma::fill::zeros);
    arma::vec rhs(d, arma::fill::zeros);
    std::map<unsigned int, arma::mat> mean_precision, site_covariance, site_chol;
    for (unsigned int j = 0; j < g; ++j) {
      unsigned int count = static_cast<unsigned int>(sizes[j]);
      if (mean_precision.count(count) == 0) {
        mean_precision[count] = inverse_spd(psi + sigma / count);
        site_covariance[count] = inverse_spd(pi + count * si);
        site_chol[count] = chol_lower(site_covariance[count]);
      }
      precision += mean_precision[count];
      rhs += mean_precision[count] * ybar.row(j).t();
    }
    arma::mat upper = chol_lower(precision).t();
    mu = arma::solve(arma::trimatu(upper), arma::solve(arma::trimatl(upper.t()), rhs)) +
      arma::solve(arma::trimatu(upper), normals(d));
    for (unsigned int j = 0; j < g; ++j) {
      unsigned int count = static_cast<unsigned int>(sizes[j]);
      arma::vec location = site_covariance[count] * (count * si * (ybar.row(j).t() - mu));
      u.row(j) = (location + site_chol[count] * normals(d)).t();
    }
    if (!fixed) {
      arma::mat residual = y;
      for (unsigned int i = 0; i < n; ++i) residual.row(i) -= mu.t() + u.row(group[i]);
      if (family == "conditional_regression") {
        sigma = factorized_covariance_draw(residual, a_sigma, regression_scale);
        psi = factorized_covariance_draw(u, a_psi, regression_scale);
      } else {
        arma::mat sigma_scale = residual.t() * residual;
        arma::mat psi_scale = u.t() * u;
        if (family == "huang_wand") {
          sigma_scale.diag() += 4 / a_sigma;
          psi_scale.diag() += 4 / a_psi;
        } else {
          sigma_scale.diag() += 0.5 * scale;
          psi_scale.diag() += 0.5 * scale;
        }
        sigma = inverse_wishart(covariance_df + n, sigma_scale);
        psi = inverse_wishart(covariance_df + g, psi_scale);
      }
      if (family == "huang_wand") a_sigma = auxiliary_draw(sigma, scale);
      if (family == "conditional_regression") a_sigma = factorized_auxiliary_draw(sigma, scale);
      if (interweave)
        scale_interweave(u, psi, ybar, sizes, mu, sigma, scale, family, regression_scale,
          scale_proposed, scale_accepted);
      // Mandatory after the collapsed-auxiliary scale move, including rejection.
      if (family == "huang_wand") a_psi = auxiliary_draw(psi, scale);
      if (family == "conditional_regression") a_psi = factorized_auxiliary_draw(psi, scale);
    }
    if (!mu.is_finite() || !y.is_finite()) Rcpp::stop("Nonfinite marginal Gibbs state.");
    for (unsigned int j = 0; j < d; ++j) {
      means(j, iteration) = mu[j];
      aux_sigma(j, iteration) = a_sigma[j];
      aux_psi(j, iteration) = a_psi[j];
    }
    save_matrix(sigma_draws, sigma, iteration);
    save_matrix(psi_draws, psi, iteration);
  }
  return Rcpp::List::create(Rcpp::Named("mu") = means,
    Rcpp::Named("sigma") = sigma_draws, Rcpp::Named("psi") = psi_draws,
    Rcpp::Named("aux_sigma") = aux_sigma, Rcpp::Named("aux_psi") = aux_psi,
    Rcpp::Named("scale_proposed") = scale_proposed,
    Rcpp::Named("scale_accepted") = scale_accepted,
    Rcpp::Named("last") = Rcpp::List::create(Rcpp::Named("mu") = mu,
      Rcpp::Named("sigma") = sigma, Rcpp::Named("psi") = psi,
      Rcpp::Named("aux_sigma") = a_sigma, Rcpp::Named("aux_psi") = a_psi,
      Rcpp::Named("y") = y, Rcpp::Named("u") = u));
  END_RCPP
}

RcppExport SEXP catmech_marginal_ratio_cpp(SEXP y_sexp, SEXP subject_sexp,
                                         SEXP mu_sexp, SEXP u_sexp,
                                         SEXP sigma_sexp, SEXP psi_sexp,
                                         SEXP j_sexp, SEXP eta_sexp, SEXP prior_sexp) {
  BEGIN_RCPP
  arma::mat y = Rcpp::as<arma::mat>(y_sexp), u = Rcpp::as<arma::mat>(u_sexp);
  arma::uvec group = Rcpp::as<arma::uvec>(subject_sexp); group -= 1;
  Rcpp::List prior(prior_sexp);
  arma::vec sizes(u.n_rows, arma::fill::zeros);
  arma::mat ybar(u.n_rows, y.n_cols, arma::fill::zeros);
  for (unsigned int i = 0; i < y.n_rows; ++i) {
    sizes[group[i]] += 1;
    ybar.row(group[i]) += y.row(i);
  }
  ybar.each_col() /= sizes;
  arma::mat psi = Rcpp::as<arma::mat>(psi_sexp);
  double regression_scale = prior.containsElementNamed("regression_scale") ?
    Rcpp::as<double>(prior["regression_scale"]) : 10;
  return Rcpp::wrap(scale_log_ratio(ybar, sizes, Rcpp::as<arma::vec>(mu_sexp), u,
    inverse_spd(Rcpp::as<arma::mat>(sigma_sexp)), inverse_spd(psi), psi,
    Rcpp::as<unsigned int>(j_sexp) - 1, Rcpp::as<double>(eta_sexp),
    Rcpp::as<double>(prior["scale"]), Rcpp::as<std::string>(prior["family"]), regression_scale));
  END_RCPP
}

RcppExport SEXP catmech_marginal_factor_cpp(SEXP values_sexp, SEXP auxiliary_sexp,
                                          SEXP regression_scale_sexp, SEXP n_sexp) {
  BEGIN_RCPP
  Rcpp::RNGScope scope;
  arma::mat values = Rcpp::as<arma::mat>(values_sexp);
  arma::vec auxiliary = Rcpp::as<arma::vec>(auxiliary_sexp);
  const unsigned int n = Rcpp::as<unsigned int>(n_sexp), d = values.n_cols;
  const double regression_scale = Rcpp::as<double>(regression_scale_sexp);
  Rcpp::NumericVector covariance(static_cast<R_xlen_t>(d) * d * n);
  covariance.attr("dim") = Rcpp::IntegerVector::create(d, d, n);
  for (unsigned int i = 0; i < n; ++i) {
    if (i % 128 == 0) Rcpp::checkUserInterrupt();
    save_matrix(covariance, factorized_covariance_draw(values, auxiliary, regression_scale), i);
  }
  return Rcpp::List::create(Rcpp::Named("covariance") = covariance);
  END_RCPP
}

RcppExport SEXP catmech_marginal_hw_cpp(SEXP n_sexp, SEXP dimension_sexp, SEXP scale_sexp) {
  BEGIN_RCPP
  Rcpp::RNGScope scope;
  const unsigned int n = Rcpp::as<unsigned int>(n_sexp), d = Rcpp::as<unsigned int>(dimension_sexp);
  const double scale = Rcpp::as<double>(scale_sexp);
  if (n < 1 || d < 1 || !R_finite(scale) || scale <= 0) Rcpp::stop("Invalid Huang-Wand prior request.");
  Rcpp::NumericVector covariance(static_cast<R_xlen_t>(d) * d * n);
  covariance.attr("dim") = Rcpp::IntegerVector::create(d, d, n);
  Rcpp::NumericMatrix auxiliary(d, n);
  for (unsigned int i = 0; i < n; ++i) {
    if (i % 128 == 0) Rcpp::checkUserInterrupt();
    arma::vec a(d);
    for (unsigned int k = 0; k < d; ++k) {
      a[k] = 1 / R::rgamma(0.5, scale);
      auxiliary(k, i) = a[k];
    }
    save_matrix(covariance, inverse_wishart(d + 1.0, arma::diagmat(4 / a)), i);
  }
  return Rcpp::List::create(Rcpp::Named("covariance") = covariance,
    Rcpp::Named("auxiliary") = auxiliary, Rcpp::Named("A") = std::sqrt(scale));
  END_RCPP
}

RcppExport SEXP catmech_marginal_predict_cpp(SEXP y_sexp, SEXP group_sexp,
                                           SEXP mu_sexp, SEXP sigma_sexp,
                                           SEXP psi_sexp, SEXP draw_sexp) {
  BEGIN_RCPP
  Rcpp::RNGScope scope;
  arma::mat y = Rcpp::as<arma::mat>(y_sexp);
  arma::uvec group = Rcpp::as<arma::uvec>(group_sexp);
  if (y.n_rows != group.n_elem || group.n_elem == 0 || group.max() >= group.n_elem)
    Rcpp::stop("Malformed marginal prediction sites.");
  return Rcpp::wrap(predict_target(y, group, Rcpp::as<arma::vec>(mu_sexp),
    Rcpp::as<arma::mat>(sigma_sexp), Rcpp::as<arma::mat>(psi_sexp), Rcpp::as<bool>(draw_sexp)));
  END_RCPP
}

RcppExport SEXP catmech_marginal_derived_cpp(SEXP input_sexp, SEXP fit_sexp,
                                           SEXP keep_sexp, SEXP center_sexp,
                                           SEXP scale_sexp) {
  BEGIN_RCPP
  Rcpp::List input(input_sexp), fit(fit_sexp);
  arma::mat y = Rcpp::as<arma::mat>(input["y"]), mu = Rcpp::as<arma::mat>(fit["mu"]);
  arma::uvec group = Rcpp::as<arma::uvec>(input["subj"]); group -= 1;
  arma::cube sigma = Rcpp::as<arma::cube>(fit["sigma"]), psi = Rcpp::as<arma::cube>(fit["psi"]);
  Rcpp::IntegerVector keep(keep_sexp);
  arma::vec center = Rcpp::as<arma::vec>(center_sexp), scale = Rcpp::as<arma::vec>(scale_sexp);
  const unsigned int d = y.n_cols, k = d - 2;
  std::vector<unsigned int> recorded, absent;
  for (unsigned int i = 0; i < y.n_rows; ++i) if (!R_finite(y(i, 0))) {
    if (R_finite(y(i, 1))) recorded.push_back(i); else absent.push_back(i);
  }
  const unsigned int base_count = k + 1 + !recorded.empty() + !absent.empty();
  const unsigned int count = base_count + 2 * (3 + 2 * k);
  const unsigned int g = group.max() + 1;
  const std::vector<Pattern> prediction_patterns = (!recorded.empty() || !absent.empty()) ?
    build_prediction_patterns(y) : std::vector<Pattern>();
  Rcpp::NumericMatrix out(keep.size(), count);
  for (int row = 0; row < keep.size(); ++row) {
    if (row % 256 == 0) Rcpp::checkUserInterrupt();
    const unsigned int it = keep[row] - 1;
    if (it >= mu.n_cols || it >= sigma.n_slices || it >= psi.n_slices)
      Rcpp::stop("Invalid marginal retained draw index.");
    arma::mat total = sigma.slice(it) + psi.slice(it);
    double intercept = center[0] + scale[0] * mu(0, it);
    if (k) {
      arma::vec beta = inverse_spd(total.submat(2, 2, d - 1, d - 1)) * total.submat(2, 0, d - 1, 0);
      for (unsigned int j = 0; j < k; ++j) {
        double coefficient = beta[j] * scale[0] / scale[j + 2];
        out(row, j + 1) = coefficient;
        intercept -= coefficient * (center[j + 2] + scale[j + 2] * mu(j + 2, it));
      }
    }
    out(row, 0) = intercept;
    if (!recorded.empty() || !absent.empty()) {
      arma::vec predicted = predict_target_patterns(y, group, g, prediction_patterns,
        mu.col(it), sigma.slice(it), psi.slice(it), false);
      unsigned int column = k + 1;
      if (!recorded.empty()) {
        double sum = 0;
        for (unsigned int index : recorded) sum += predicted[index];
        out(row, column++) = center[0] + scale[0] * sum / recorded.size();
      }
      if (!absent.empty()) {
        double sum = 0;
        for (unsigned int index : absent) sum += predicted[index];
        out(row, column++) = center[0] + scale[0] * sum / absent.size();
      }
    }
    unsigned int column = base_count;
    for (unsigned int component = 0; component < 2; ++component) {
      ConditionalCovariance block = conditional_covariance(component ? psi.slice(it) : sigma.slice(it));
      out(row, column++) = block.residual(0, 0);
      out(row, column++) = block.residual(1, 0);
      out(row, column++) = block.residual(1, 1);
      for (unsigned int predictor = 0; predictor < k; ++predictor) {
        out(row, column++) = block.regression(predictor, 0);
        out(row, column++) = block.regression(predictor, 1);
      }
    }
  }
  return out;
  END_RCPP
}
