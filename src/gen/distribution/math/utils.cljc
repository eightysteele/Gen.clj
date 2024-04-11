(ns gen.distribution.math.utils
  "Math utilities and constants.")

(def ^:no-doc log-pi
  (Math/log Math/PI))

(def ^:no-doc log-2pi
  (Math/log (* 2 Math/PI)))

(def ^:no-doc sqrt-2pi
  (Math/sqrt (* 2 Math/PI)))

(def ^:no-doc half-log-2pi
  (* 0.5 (Math/log (* 2.0 Math/PI))))
