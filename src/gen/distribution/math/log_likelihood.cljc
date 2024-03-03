(ns gen.distribution.math.log-likelihood
  "Log-likelihood implementations for various primitive distributions."
  (:require [gen.distribution.math.gamma :as g]))

(defn bernoulli
  "Returns the log-likelihood of a [Bernoulli
  distribution](https://en.wikipedia.org/wiki/Bernoulli_distribution)
  parameterized by probability `p` at the boolean value `v`."
  [p v]
  {:pre [(<= 0 p 1)]}
  (Math/log (if v p (- 1.0 p))))

(defn gamma
  "Returns the log-likelihood of the [Gamma
  distribution](https://en.wikipedia.org//wiki/Gamma_distribution) parameterized
  by `shape` and `scale` at the value `v`.

  The implementation follows the algorithm described on the Gamma
  distribution's [Wikipedia
  page](https://en.wikipedia.org//wiki/Gamma_distribution#Maximum_likelihood_estimation)."
  [shape scale v]
  (if (pos? v)
    (- (* (dec shape) (Math/log v))
       (/ v scale)
       (g/log-gamma-fn shape)
       (* shape (Math/log scale)))
    ##-Inf))

(defn beta
  "Returns the log-likelihood of the [Beta
  distribution](https://en.wikipedia.org/wiki/Beta_distribution) parameterized by
  `alpha` and `beta` at the value `v`.

  The implementation follows the algorithm described on the Beta
  distribution's [Wikipedia
  page](https://en.wikipedia.org/wiki/Beta_distribution#Probability_density_function)."
  [alpha beta v]
  {:pre [(pos? alpha) (pos? beta)]}
  (if (< 0 v 1)
    (- (+ (* (- alpha 1) (Math/log v))
          (* (- beta 1) (Math/log (- 1 v))))
       (g/log-beta-fn alpha beta))
    ##-Inf))

(defn cauchy
  "Returns the log-likelihood of a [Cauchy
  distribution](https://en.wikipedia.org/wiki/Cauchy_distribution) parameterized
  by `location` and `scale` at the value `v`.

  The implementation follows the algorithm described on the Cauchy
  distribution's [Wikipedia
  page](https://en.wikipedia.org/wiki/Cauchy_distribution#Probability_density_function_(PDF))."
  [location scale v]
  (let [normalized (/ (- v location) scale)
        norm**2    (* normalized normalized)]
    (- (- g/log-pi)
       (Math/log scale)
       (Math/log (+ 1 norm**2)))))

(defn delta
  "Returns the log-likelihood of the [Dirac delta
  distribution](https://en.wikipedia.org/wiki/Dirac_delta_function) centered
  around `center` at the value `v`."
  [center v]
  (if (= center v) 0.0 ##-Inf))

(defn exponential
  "Returns the log-likelihood of the [exponential
  distribution](https://en.wikipedia.org/wiki/Exponential_distribution) with
  rate parameter `rate` at the value `v`."
  [rate v]
  (if (>= v 0)
    (- (Math/log rate) (* rate v))
    ##-Inf))

(defn laplace
  "Returns the log-likelihood of the [Laplace
  distribution](https://en.wikipedia.org/wiki/Laplace_distribution) with
  `location` and `scale` parameters at the value `v`.

  The implementation follows the algorithm described on the Laplace
  distribution's [Wikipedia
  page](https://en.wikipedia.org/wiki/Laplace_distribution#Probability_density_function)."
  [location scale v]
  (- (+ (Math/log (* 2.0 scale))
        (/ (Math/abs ^double (- v location))
           scale))))

(defn gaussian
  "Returns the log-likelihood of the [Gaussian
  distribution](https://en.wikipedia.org/wiki/Gaussian_distribution) with
  mean `mu` and standard deviation `sigma` at the value `v`.

  The implementation follows the algorithm described on the Gaussian
  distribution's [Wikipedia
  page](https://en.wikipedia.org/wiki/Normal_distribution#Operations_on_a_single_normal_variable):

  Given $z = \\left(\\frac{x-\\mu}{\\sigma}\\right)^2$:

  $$
  \\begin{aligned}
  \\ln p(x) &= -\\frac{1}{2}z^2 - \\ln (\\sigma \\sqrt{2\\pi}) \\\\
            &= -\\frac{1}{2}\\left(z^2 + 2 \\ln (\\sigma \\sqrt{2\\pi}) \\right) \\\\
            &= -\\frac{1}{2}\\left(z^2 + \\ln (\\sigma^2) + \\ln (2\\pi) \\right)
  \\end{aligned}
  $$"
  [mu sigma v]
  (let [v-mu:sigma (/ (- v mu) sigma)]
    (* -0.5 (+ g/log-2pi
               (* 2 (Math/log sigma))
               (* v-mu:sigma
                  v-mu:sigma)))))

(defn uniform
  "Returns the log-likelihood of the continuous [uniform
  distribution](https://en.wikipedia.org/wiki/Continuous_uniform_distribution)
  with inclusive lower bound `a` and inclusive upper bound `b` at the value
  `v`."
  [a b v]
  (if (<= a v b)
    (- (Math/log (- b a)))
    ##-Inf))

(defn student-t
  "Returns the log-likelihood of the [non-standardized Student's
  t-distribution](https://en.wikipedia.org/wiki/Student's_t-distribution#Location-scale_transformation)
  parametrized by `location`, `scale` and degrees-of-freedom `nu` at the value
  `v`.

  This distribution is also known as the location-scale t-distribution.

  The implementation follows the algorithm described on the
  distribution's [Wikipedia
  page](https://en.wikipedia.org/wiki/Student's_t-distribution#Location-scale_transformation)."
  ([nu v] (student-t nu 0 1 v))
  ([nu location scale v]
   (let [inc-nu      (inc nu)
         half-inc-nu (* 0.5 inc-nu)
         normalized  (/ (- v location) scale)
         norm**2     (* normalized normalized)]
     (- (g/log-gamma-fn half-inc-nu)
        (g/log-gamma-fn (* 0.5 nu))
        (Math/log scale)
        (* 0.5 (+ g/log-pi (Math/log nu)))
        (* half-inc-nu (Math/log (inc (/ norm**2 nu))))))))
