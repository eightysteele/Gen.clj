(ns gen.distribution.math.gamma.gamma-test
  (:require [com.gfredericks.test.chuck.clojure-test :refer [checking]]
            [clojure.test :refer [deftest is testing]]
            [gen.distribution.math.gamma :as g]
            [gen.generators :refer [gen-double within]]
            [same.core :refer [ish? with-comparator]]))

(def log-pi
  (Math/log Math/PI))

(defn factorial
  "Factorial implementation for testing."
  [n]
  (if (zero? n)
    1
    (* n (factorial (dec n)))))

(deftest log-gamma-fn-tests
  (testing "log-Gamma ~matches log(factorial)"
    (with-comparator (within 1e-11)
      (doseq [n (range 1 15)]
        (is (ish? (Math/log (factorial (dec n)))
                  (g/log-gamma-fn n))))))


  (with-comparator (within 1e-12)
    (checking "Euler's reflection formula"
              [z (gen-double 0.001 0.999)]
              (is (ish? (+ (g/log-gamma-fn (- 1 z))
                           (g/log-gamma-fn z))
                        (- log-pi
                           (Math/log
                            (Math/sin (* Math/PI z)))))))))
