(ns gen.distribution.math.gamma)

;; These come in handy in the implementations below and are worth caching.

(def ^:no-doc log-pi
  (Math/log Math/PI))

(def ^:no-doc log-2pi
  (Math/log (* 2 Math/PI)))

(def ^:no-doc sqrt-2pi
  (Math/sqrt (* 2 Math/PI)))

;; ## Log-likelihood implementations

(def ^:no-doc gamma-coefficients
  "Coefficients for the Lanczos approximation to the natural log of the Gamma
  function described in [section 6.1 of Numerical
  Recipes](http://phys.uri.edu/nigh/NumRec/bookfpdf/f6-1.pdf)."
  [76.18009172947146
   -86.50532032941677
   24.01409824083091
   -1.231739572450155
   0.1208650973866179e-2
   -0.5395239384953e-5])

(defn ^:no-doc log-gamma-fn
  "Returns the natural log of the value of the [Gamma
  function](https://en.wikipedia.org/wiki/Gamma_function) evaluated at `x`

  This function implements the Lanczos approximation described in [section 6.1
  of Numerical Recipes](http://phys.uri.edu/nigh/NumRec/bookfpdf/f6-1.pdf)."
  [x]
  (let [tmp (+ x 5.5)
        tmp (- (* (+ x 0.5) (Math/log tmp)) tmp)
        n   (dec (count gamma-coefficients))
        ser (loop [i   0
                   x+1 (inc x)
                   acc 1.000000000190015]
              (if (> i n)
                acc
                (let [coef (nth gamma-coefficients i nil)]
                  (recur (inc i)
                         (inc x+1)
                         (+ acc (/ coef x+1))))))]
    (+ tmp (Math/log (* sqrt-2pi (/ ser x))))))

(defn log-beta-fn
  "Returns the natural log of the value of the [Beta
  function](https://en.wikipedia.org/wiki/Beta_function) evaluated at inputs `a`
  and `b`."
  [a b]
  (- (+ (log-gamma-fn a)
        (log-gamma-fn b))
     (log-gamma-fn (+ a b))))
