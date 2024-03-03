(ns gen.distribution.math.gamma.gamma-test
  (:require [com.gfredericks.test.chuck.clojure-test :refer [checking]]
            [clojure.test :refer [deftest is testing]]
            [gen.distribution.math.gamma :as g]
            [gen.distribution.math.combinatorics.factorial :as f]
            [gen.generators :refer [gen-double within]]
            [same.core :refer [ish? with-comparator]])
  (:import (org.apache.commons.numbers.gamma LogGamma)))

;; 170! is the max factorial possible in a double
(def MAX_N_DOUBLE 170)

(def log-pi
  (Math/log Math/PI))

(deftest log-gamma-fn-tests
  (testing "log-Gamma ~matches log(factorial)"
    (with-comparator (within 1e-12)
      (doseq [n (range MAX_N_DOUBLE)]
        (is (ish? (Math/log (f/factorial (dec n)))
                  (g/log-gamma-fn n))))))

  (testing "log-gamma-fn ~matches commons.numbers.gamma.LogGamma"
    (with-comparator (within 1e-12)
      (doseq [n (range MAX_N_DOUBLE)]
        (is (ish? (LogGamma/value(n)
                  (g/log-gamma-fn n)))))))

  (with-comparator (within 1e-12)
    (checking "Euler's reflection formula"
              [z (gen-double 0.001 0.999)]
              (is (ish? (+ (g/log-gamma-fn (- 1 z))
                           (g/log-gamma-fn z))
                        (- log-pi
                           (Math/log
                            (Math/sin (* Math/PI z)))))))))
