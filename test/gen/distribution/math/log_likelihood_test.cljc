(ns gen.distribution.math.log-likelihood-test
  (:require [com.gfredericks.test.chuck.clojure-test :refer [checking]]
            [clojure.test :refer [deftest is]]
            [gen.distribution.math.log-likelihood :as ll]
            [gen.distribution :as distribution]
            [gen.distribution-test :as dt]
            [gen.generators :refer [gen-double within]]
            [same.core :refer [ish? with-comparator]]))

(defn ->logpdf [f]
  (fn [& args]
    (reify distribution/LogPDF
      (logpdf [_ v]
        (apply f (concat args [v]))))))

(deftest gamma-tests
  (dt/gamma-tests (->logpdf ll/gamma)))

(deftest beta-tests
  (dt/beta-tests (->logpdf ll/beta)))

(deftest bernoulli-tests
  (dt/bernoulli-tests (->logpdf ll/bernoulli)))

(deftest cauchy-tests
  (dt/cauchy-tests (->logpdf ll/cauchy)))

(deftest delta-tests
  (dt/delta-tests (->logpdf ll/delta)))

(deftest exponential-tests
  (dt/exponential-tests (->logpdf ll/exponential))

  (checking "exponential will never produce negative values"
            [rate (gen-double -100 100)
             v    (gen-double -100 -0.00001)]
            (is (= ##-Inf (ll/exponential rate v))))

  (checking "rate 0.0 produces #-Inf"
            [v (gen-double -100 100)]
            (is (= ##-Inf (ll/exponential 0.0 v)))))

(deftest laplace-test
  (dt/laplace-tests (->logpdf ll/laplace)))

(deftest gaussian-tests
  (dt/normal-tests (->logpdf ll/gaussian)))

(deftest uniform-tests
  (dt/uniform-tests (->logpdf ll/uniform)))

(deftest student-t-tests
  (dt/student-t-tests (->logpdf ll/student-t)))
