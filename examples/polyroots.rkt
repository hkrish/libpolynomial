#lang typed/racket

(require typed/racket/unsafe
         math/flonum)

(unsafe-require/typed
 "defs.rkt"
 (poly_quadroots (-> Float Float Float (Listof Float)))
 (poly_cubicroots (-> Float Float Float Float (Listof Float)))
 (poly_quartroots (-> Float Float Float Float Float (Listof Float)))
 (poly_quintroots (-> Float Float Float Float Float Float (Listof Float))))

(unsafe-require/typed
 "utils.rkt"
 (random-polys/quad (-> Number (Listof (List Float Float Float))))
 (random-polys/cubic (-> Number (Listof (List Float Float Float Float))))
 (random-polys/quart (-> Number (Listof (List Float Float Float Float Float))))
 (random-polys/quint (-> Number (Listof (List Float Float Float Float Float Float))))
 ;; (testdata/cubic (Listof (List (List Float Float Float Float) (Listof Float))))
 (testdata/cubic (Listof (List Float Float Float Float)))
 (displayln-digit-diff (->* (Any Any) (#:prefix-a Any #:prefix-b Any) Void)))



(: polyroots/quad (-> Float Float Float (Listof Float)))
(define (polyroots/quad a B c)
  (cond
    [(< (abs a) epsilon.0)
     (if (< (abs B) epsilon.0)
         '()
         (list (/ (- c) B)))]
    [else
     (let-values ([(a b c D) (let* ([b (/ B -2.0)]
                                    [d (- (* b b) (* a c))])
                               (if (or (<= (abs d) 1e-8) (>= (abs d) 1e8))
                                   (let* ([infnorm0 (max (abs a) (abs b) (abs c))]
                                          [infnorm (if (= 0. infnorm0) epsilon.0 infnorm0)]
                                          [scale (flexp2 (- (flround (fllog2 infnorm))))]
                                          [a (* a scale)]
                                          [b (* b scale)]
                                          [c (* c scale)])
                                     (values a b c (- (* b b) (* a c))))
                                   (values a b c d)))])
       (if (>= D (- epsilon.0))
           (let* ([Q (if (< D 0.) 0. (flsqrt D))]
                  [R (+ b (* (if (< b 0.) -1. 1.) Q))])
             (define-values (x1 x2) (if (= 0. R)
                                        (values (/ c a) (/ (- c) a))
                                        (values (/ R a) (/ c R))))
             (cond
               [(<= (abs (- x1 x2)) epsilon.0) (list x1)]
               [else (if (< x1 x2) (list x1 x2) (list x2 x1))]))
           '()))]))

(: -NR-rec/cubic (-> Float Float Float Float Float Float (Values Float Float Float)))
(define (-NR-rec/cubic a b c d s x)
  (define ec (+ 1. epsilon.0))
  (define ax (* a x))
  (define b1 (+ ax b))
  (define c2 (+ (* b1 x) c))
  (define qd (+ (* (+ ax b1) x) c2))
  (define q (+ (* c2 x) d))
  (define xn (if (= qd 0.) x (- x (/ (/ q qd) ec))))
  (if (> (* xn s) (* x s))
      (-NR-rec/cubic a b c d s xn)
      (values x b1 c2)))

(: polyroots/cubic (-> Float Float Float Float (Listof Float)))
(define (polyroots/cubic a0 b0 c0 d0)
  (define infnorm (max (abs a0) (max (abs b0) (max (abs c0)  (abs d0)))))
  (define-values (a b c d)
    (if (or (<= (abs infnorm) 1e-7) (>= (abs infnorm) 1e7))
        (let* ([scale (flexp2 (- (flround (fllog2 (if (= 0. infnorm) epsilon.0 infnorm)))))]
               [a (* a0 scale)]
               [b (* b0 scale)]
               [c (* c0 scale)]
               [d (* d0 scale)])
          (values a b c d))
        (values a0 b0 c0 d0)))
  (define-values (rs qa qb qc)
    (cond
      [(= 0. a) (values '() b c d)]
      [(= 0. d) (values '(0.) a b c)]
      [else
       (let* ([x (/ (/ (- b) a) 3.)]
              [ax (* a x)]
              [b1 (+ ax b)]
              [c2 (+ (* b1 x) c)]
              [qd (+ (* (+ ax b1) x) c2)]
              [q (+ (* c2 x) d)]
              [t0 (/ q a)]
              [r0 (flexpt (flabs t0) (/ 1. 3.))]
              [s (if (< t0 0.)  -1.0  1.0)]
              [t (/ (- qd) a)]
              ; λ ~= 1.324717957244746025960909 in λ^3 = λ + 1
              ; See Kahan's notes on why 1.324718*... works.
              [r (if (> t 0.) (* 1.324717957244746025960909 (flmax r0 (flsqrt t))) r0)]
              [x0 (- x (* s r))])
         (if (not (= x0 x))
             (let ()
               (define-values (x b1 c2) (-NR-rec/cubic a b c d s x0))
               (if (> (* (* (abs a) x) x) (abs (/ d x)))
                   (let ([c2 (/ (- d) x)]
                         [b1 (/ (- c2 c) x)])
                     (values (list x) a b1 c2))
                   (values (list x) a b1 c2)))
             (values (list x) a b1 c2)))]))
  (remove-duplicates
   (append rs (polyroots/quad qa qb qc))
   (lambda ((a : Float) (b : Float)) (<= (abs (- a b)) epsilon.0))))



(: test-perf/ffi (case->
                  (-> (-> Float Float Float (Listof Float))
                      (Listof (List Float Float Float))
                      Real Void)
                  (-> (-> Float Float Float Float (Listof Float))
                      (Listof (List Float Float Float Float))
                      Real Void)
                  (-> (-> Float Float Float Float Float (Listof Float))
                      (Listof (List Float Float Float Float Float))
                      Real Void)
                  (-> (-> Float Float Float Float Float Float (Listof Float))
                      (Listof (List Float Float Float Float Float Float))
                      Real Void)))
(define (test-perf/ffi fn polys N)
  (define count (length polys))
  (define times  (inexact->exact (ceiling (/ N count))))
  (time (for* ([i (in-range 0 times)]
               [c (in-list polys)])
          (apply fn c))))


(: test-poly-results/ffi (case->
                          (-> (-> Float Float Float (Listof Float))
                              (-> Float Float Float (Listof Float))
                              (Listof (List Float Float Float))
                              Void)
                          (-> (-> Float Float Float Float (Listof Float))
                              (-> Float Float Float Float (Listof Float))
                              (Listof (List Float Float Float Float))
                              Void)
                          (-> (-> Float Float Float Float Float (Listof Float))
                              (-> Float Float Float Float Float (Listof Float))
                              (Listof (List Float Float Float Float Float))
                              Void)
                          (-> (-> Float Float Float Float Float Float (Listof Float))
                              (-> Float Float Float Float Float Float (Listof Float))
                              (Listof (List Float Float Float Float Float Float))
                              Void)))
(define (test-poly-results/ffi fn ffifn polys)
  (define accuracy (* 10. epsilon.0))
  (for ([ca (in-list polys)])
    (define rs (sort (apply fn ca) fl<))
    (define as (sort (apply ffifn ca) fl<))
    (if (= (length rs) (length as))
        (when (not (for/and : Boolean ([r (in-list rs)]
                                       [a (in-list as)])
                     (<= (abs (- r a)) accuracy)))
          (define ulp-errors (for/list : (Listof Float) ([r (in-list rs)]
                                                         [a (in-list as)])
                               (flulp-error r a)))
          (displayln-digit-diff as rs #:prefix-a "Accuracy ")
          (displayln (~a "Ulp-err  " ulp-errors)))
        (displayln-digit-diff as rs #:prefix-a "Mismatch: Number of Roots: "))))




(let* ([times 100000]
       [count 10000]
       [quads/test (random-polys/quad 10)]
       [cubics/test (random-polys/cubic 10)]
       [quads (random-polys/quad count)]
       [cubics (random-polys/cubic count)]
       [quarts (random-polys/quart count)]
       [quints (random-polys/quint count)]
       )
  (displayln "--------------------------------------------------")
  (displayln "Accuracy Tests:")

  (displayln "Quadratic FFI vs Racket:")
  (test-poly-results/ffi polyroots/quad poly_quadroots quads/test)

  (displayln "Cubic FFI vs Racket/1:")
  (test-poly-results/ffi polyroots/cubic poly_cubicroots testdata/cubic)
  (displayln "Cubic FFI vs Racket/2:")
  (test-poly-results/ffi polyroots/cubic poly_cubicroots cubics/test)

  (displayln "--------------------------------------------------")
  (displayln "Performace tests")

  (collect-garbage)
  (display "FFI/Quad            \t")
  (test-perf/ffi poly_quadroots quads times)
  (collect-garbage)
  (display "Typed/Racket/Quad   \t")
  (test-perf/ffi polyroots/quad quads times)

  (collect-garbage)
  (display "FFI/Cubic           \t")
  (test-perf/ffi poly_cubicroots cubics times)
  (collect-garbage)
  (display "Typed/Racket/Cubic  \t")
  (test-perf/ffi polyroots/cubic cubics times)

  (collect-garbage)
  (display "FFI/Quart           \t")
  (test-perf/ffi poly_quartroots quarts times)

  (collect-garbage)
  (display "FFI/Quint           \t")
  (test-perf/ffi poly_quintroots quints times))
