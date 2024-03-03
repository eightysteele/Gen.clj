(ns gen.distribution.math.combinatorics.factorial-test
  (:require [gen.distribution.math.combinatorics.factorial :as f]
            [clojure.test :refer [deftest is testing]]))

(def MAX_N_DOUBLE 170)

(deftest factorial-zero-tests
  (testing "factorial of zero"
    (is (= 1.0 (f/factorial 0)))))

(deftest factorial-too-big-tests
  (testing "factorial of n > 170 should be positive infinity"
    (is (= ##Inf (f/factorial (inc MAX_N_DOUBLE))))))

(deftest factorial-all-tests
  (testing "factorial n from 1 to 170"
    (doseq [n (range MAX_N_DOUBLE)]
            (let [expected (double  (reduce *' 1 (range 1 (inc n))))]
        (is (= expected (f/factorial n))
            (str n "! factorial mismatch"))))))
