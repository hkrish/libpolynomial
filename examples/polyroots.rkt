#lang typed/racket

(require math/flonum
         typed/racket/unsafe
         racket/performance-hint)

(require/typed profile
  (profile-thunk (-> (-> Any) Any)))

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
 (testdata/quart (Listof (List Float Float Float Float Float)))
 (testdata/quint (Listof (List Float Float Float Float Float Float)))
 (displayln-digit-diff (->* (Any Any) (#:prefix-a Any #:prefix-b Any) Void)))


(: polyroots/quad (-> Float Float Float (Listof Float)))
(define (polyroots/quad a B c)
  (define-values (an bn cn D)
    (let* ([b (/ B -2.0)]
           [d (- (* b b) (* a c))])
      (if (or (<= (abs d) 1e-8) (>= (abs d) 1e8))
          ;; Normalise coefficients à la *Jenkins & Traub's RPOLY*
          ;; so that all the bits in the mantissa are preserved.
          (let* ([infnorm0 (max (abs a) (abs b) (abs c))]
                 [infnorm (if (= 0. infnorm0) epsilon.0 infnorm0)]
                 [scale (flexp2 (- (flround (fllog2 infnorm))))]
                 [a (* a scale)]
                 [b (* b scale)]
                 [c (* c scale)])
            (values a b c (- (* b b) (* a c))))
          (values a b c d))))
  (cond
    [(< (abs an) epsilon.0)
     (if (< (abs bn) epsilon.0)
         '()
         (list (/ (- cn) B)))]
    [else
     (if (>= D (- epsilon.0))
         (let* ([Q (if (< D 0.) 0. (flsqrt D))]
                [R (+ bn (* (if (< bn 0.) -1. 1.) Q))])
           (define-values (x1 x2) (if (= 0. R)
                                      (values (/ cn an) (/ (- cn) an))
                                      (values (/ R an) (/ cn R))))
           (cond
             [(<= (abs (- x1 x2)) epsilon.0) (list x1)]
             [else (if (< x1 x2) (list x1 x2) (list x2 x1))]))
         '())]))


(: polyroots/cubic (-> Float Float Float Float (Listof Float)))
(define (polyroots/cubic a0 b0 c0 d0)
  (define infnorm (max (abs a0) (abs b0) (abs c0) (abs d0)))
  (define-values (a b c d)
    (if (or (<= (abs infnorm) 1e-7) (>= (abs infnorm) 1e7))
        ;; Normalise coefficients à la *Jenkins & Traub's RPOLY*
        ;; so that all the bits in the mantissa are preserved.
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
               (define-values (x b1 c2)
                 (let loop : (Values Float Float Float)
                      ([x x0])
                      (define ec (+ 1. epsilon.0))
                      (define ax (* a x))
                      (define b1 (+ ax b))
                      (define c2 (+ (* b1 x) c))
                      (define qd (+ (* (+ ax b1) x) c2))
                      (define q (+ (* c2 x) d))
                      (define xn (if (= qd 0.) x (- x (/ (/ q qd) ec))))
                      (if (> (* xn s) (* x s))
                          (loop xn)
                          (values x b1 c2))))
               (if (> (* (* (abs a) x) x) (abs (/ d x)))
                   (let ([c2 (/ (- d) x)]
                         [b1 (/ (- c2 c) x)])
                     (values (list x) a b1 c2))
                   (values (list x) a b1 c2)))
             (values (list x) a b1 c2)))]))
  (remove-duplicates
   (append rs (polyroots/quad qa qb qc))
   (lambda ((a : Float) (b : Float)) (<= (abs (- a b)) epsilon.0))))


(: -root-polish/quart (-> Float Float Float Float Float Float Float))
(define (-root-polish/quart a b c d e x0)
  (define da (* a 4))
  (define db (* b 3))
  (define dc (* c 2))
  (define z (+ e (* x0 (+ d (* x0 (+ c (* x0 (+ b (* x0 a)))))))))
  (let loop ([x x0] [z z] [n-rec 10])
    (define dz (+ d (* x (+ dc (* x (+ db (* x da)))))))
    (define nx (if (= 0. dz) x (- x (/ z dz))))
    (define nz (+ e (* nx (+ d (* nx (+ c (* nx (+ b (* nx a)))))))))
    (if (and (not (= nx x)) (< (abs nz) (abs z)) (> (abs nz) epsilon.0) (> n-rec 0))
        (loop nx nz (sub1 n-rec))
        x)))


(: polyroots/quart (->* (Float Float Float Float Float)
                        (#:polish-roots? Boolean)
                        (Listof Float)))
(define (polyroots/quart a0 b0 c0 d0 e0 #:polish-roots? (polish? #t))
  (define infnorm (max (abs a0) (abs b0) (abs c0) (abs d0)))
  (define-values (an bn cn dn en)
    (if (or (<= (abs infnorm) 1e-7) (>= (abs infnorm) 1e7))
        ;; Normalise coefficients à la *Jenkins & Traub's RPOLY*
        ;; so that all the bits in the mantissa are preserved.
        (let* ([scale (flexp2 (- (flround (fllog2 (if (= 0. infnorm) epsilon.0 infnorm)))))]
               [a (* a0 scale)]
               [b (* b0 scale)]
               [c (* c0 scale)]
               [d (* d0 scale)]
               [e (* e0 scale)])
          (values a b c d e))
        (values a0 b0 c0 d0 e0)))
  (define rs (cond
               [(< (abs an) epsilon.0) (polyroots/cubic bn cn dn en)]
               [(< (abs en) epsilon.0) (polyroots/cubic an bn cn dn)]
               [else
                ;; Original quartic `P(x) = ac x^4 + bc x^3 + cc x^2 + dc x + ec = 0`
                ;; Convert to monic form `P(x) → x^4 + a x^3 + b x^2 + c x + d = 0`
                (define a (/ bn an))
                (define b (/ cn an))
                (define c (/ dn an))
                (define d (/ en an))
                ;; Powers of a
                (define a2 (* a a))
                (define a3 (* a2 a))
                (define a4 (* a3 a))
                (define a_4d (/ a 4.0))
                ;; Substitute `x → y-a/4`  =>  `Q(y) → y^4 + e y^2 + f y + g = 0`
                (define e (-  b (/ (* 3.0 a2) 8.0)))
                (define f (+ c (/ a3 8.0) (/ (* a b) -2.0)))
                (define g (+ d (/ (* 3.0 a4) -256.0) (/ (* a2 b) 16.0) (/ (* a c) -4.0)))
                (map (lambda ((r : Float)) (- r a_4d))
                     (cond
                       [(<= (abs g) epsilon.0)
                        ;; One root @ 0
                        (cons 0.
                              (if (<= (abs f) epsilon.0)
                                  ;; `g~=0, f~=0`  ⇒  solve `y^2 + e = 0`
                                  (if (and (> (abs e) epsilon.0) (< e 0.))
                                      (let ([qr1 (sqrt (- e))])
                                        (list qr1 (- qr1)))
                                      '())
                                  ;; `g ~= 0`  ⇒  solve `y^3 + ey + f = 0`
                                  (polyroots/cubic 1.0 0.0 e f)))]
                       [(<= (abs f) epsilon.0)
                        ;; `f ~= 0`  ⇒  `y^4 + e y^2 + g`  is a quadratic in y^2
                        (for/fold ([l : (Listof Float) '()])
                                  ([r (in-list (polyroots/quad 1. e g))])
                          (cons r (cons (- r) l)))]
                       [else
                        ;; Coefficients of cubic in h^2
                        (define b1 (* 2. e))
                        (define c1 (- (* e e) (* 4. g)))
                        (define d1 (- (* f f)))
                        (define h2 (apply max (polyroots/cubic 1. b1 c1 d1)))
                        (define h (flsqrt h2))
                        (define j (/ (- (+ e h2) (/ f h)) 2.))
                        (append (polyroots/quad 1. h j)
                                (polyroots/quad 1. (- h) (/ g j)))]))]))
  (define rsu (remove-duplicates rs (lambda ((a : Float) (b : Float))
                                      (<= (abs (- a b)) epsilon.0))))
  (if polish?
      (for/list ([r (in-list rsu)])
        (-root-polish/quart an bn cn dn en r))
      rsu))


(begin-encourage-inline

  (: -eval-quint (-> Float Float Float Float Float Float Float (Values Float Float)))
  (define (-eval-quint a b c d e f t)
    (define val (+ f (* t (+ e (* t (+ d (* t (+ c (* t (+ b (* t a)))))))))))
    (define e1 (* (+ (abs b) (* t (abs a))) epsilon.0))
    (define e2 (+ (* (+ (abs c) (* t (abs b))) epsilon.0) (* t e1)))
    (define e3 (+ (* (+ (abs d) (* t (abs c))) epsilon.0) (* t e2)))
    (define e4 (+ (* (+ (abs e) (* t (abs d))) epsilon.0) (* t e3)))
    (define error (+ (* (+ (abs f) (* t (abs e))) epsilon.0) (* t e4)))
    (values val error))


  (: -root-polish/quint (-> Float Float Float Float Float Float Float Float))
  (define (-root-polish/quint a b c d e f x0)
    (define da (* a 5))
    (define db (* b 4))
    (define dc (* c 3))
    (define dd (* d 2))
    (define z (+ f (* x0 (+ e (* x0 (+ d (* x0 (+ c (* x0 (+ b (* x0 a)))))))))))
    (define max-rec 10)
    (let loop ([x x0] [z z] [n-rec 10])
      (define dz (+ e (* x (+ dd (* x (+ dc (* x (+ db (* x da)))))))))
      (define nx (if (= 0. dz) x (- x (/ z dz))))
      (define nz (+ f (* nx (+ e (* nx (+ d (* nx (+ c (* nx (+ b (* nx a)))))))))))
      (if (and (not (= nx x)) (< (abs nz) (abs z)) (> (abs nz) epsilon.0) (> n-rec 0))
          (loop nx nz (sub1 n-rec))
          (if (< (abs nz) (abs z)) nx x)))))


;; Refines a root (root search) in an interval known to contain exactly one root
;; of a quintic polynomial.
;;
;; cf. "Bisected Direct Quadratic Regula Falsi"
;; - Robert G. Gottlieb and Blair F. Thompson
;;
;; I have stolen the tolerance bound and stopping criterion from (flbracketed-root)
;;
(: -bracketed-root/quint (-> Float Float Float Float Float Float Float Float Float Float Float))
(define (-bracketed-root/quint a b c d e f x1 y1 x2 y2)
  (define x-last0 : Float (if (< y1 y2) x1 x2))
  (define-values (xup0 xdn0 yup0 ydn0) (if (< y1 0.)
                                           (values x2 x1 y2 y1)
                                           (values x1 x2 y1 y2)))
  (define-values (root y xup xdn yup ydn)
    (for/fold : (Values Float Float Float Float Float Float)
        ([x-last x-last0]
         [y y1]
         [xup xup0]
         [xdn xdn0]
         [yup yup0]
         [ydn ydn0])
        ([i (in-range 0 50)])
      (define min-abs-xupxdn (min (abs xup) (abs xdn)))
      (define ε (if (min-abs-xupxdn . <= . +max-subnormal.0) +min.0 (* min-abs-xupxdn epsilon.0)))
      #:break (or (= y 0.) (<= (abs (- xdn xup)) ε))
      (define delta  (/ (- xup xdn) 2.))
      (define xm     (/ (+ xup xdn) 2.))
      (define ym     (+ f (* xm (+ e (* xm (+ d (* xm (+ c (* xm (+ b (* xm a)))))))))))
      (define a_d    (* 2. delta delta))
      (define a_     (/ (+ yup ydn (* -2. ym)) a_d))
      (define b_     (/ (- yup ydn) (* 2. delta)))
      (define xden   (flsqrt (- 1. (* (/ (* 4. a_) b_) (/ ym b_)))))
      (define x      (- xm (/ (* 2. ym) (* b_ (+ 1. xden)))))
      #:break (= x x-last)
      (define yn     (+ f (* x (+ e (* x (+ d (* x (+ c (* x (+ b (* x a)))))))))))
      (define-values (xupn xdnn yupn ydnn)
        (cond
          [(and (> yn 0.) (< ym 0.)) (values x xm yn ym)]
          [(> yn 0.)                 (values x xdn yn ydn)]
          [(> ym 0.)                (values xm x ym yn)]
          [else                     (values xup x yup yn)]
          ))
      (values x yn xupn xdnn yupn ydnn)))
  root)


(: polyroots/quint (-> Float Float Float Float Float Float (Listof Float)))
(define (polyroots/quint a0 b0 c0 d0 e0 f0)
  (define infnorm (max (abs a0) (abs b0) (abs c0) (abs d0) (abs f0)))
  (define-values (a0l b0l c0l d0l e0l f0l)
    (if (or (<= (abs infnorm) 1e-7) (>= (abs infnorm) 1e7))
        ;; Normalise coefficients à la *Jenkins & Traub's RPOLY*
        ;; so that all the bits in the mantissa are preserved.
        (let* ([scale (flexp2 (- (flround (fllog2 (if (= 0. infnorm) epsilon.0 infnorm)))))]
               [a (* a0 scale)]
               [b (* b0 scale)]
               [c (* c0 scale)]
               [d (* d0 scale)]
               [e (* e0 scale)]
               [f (* f0 scale)])
          (values a b c d e f))
        (values a0 b0 c0 d0 e0 f0)))
  (define rs
    (cond
      [(< (abs a0l) epsilon.0)
       (polyroots/quart b0l c0l d0l e0l f0l #:polish-roots? #f)]
      [(= 0. f0l)
       (cons 0. (polyroots/quart a0l b0l c0l d0l e0l #:polish-roots? #f))]
      [else
       ;; Convert to monic form `P(x) → x^5 + ac x^4 + bc x^3 + cc x^2 + dc x + ec = 0`
       (define ac (/ b0l a0l))
       (define bc (/ c0l a0l))
       (define cc (/ d0l a0l))
       (define dc (/ e0l a0l))
       (define ec (/ f0l a0l))
       ;; Convert to a *depressed* quintic
       ;; `x → y-ac/5` ⇒ `Q(y) → y^5 + a y^3 + b y^2 + c y + d`
       (define a_5d (/ ac 5.0))
       (define ac2 (* ac ac))
       (define ac3 (* ac2 ac))
       (define cc25 (* 25.0 cc))
       (define acbc (* ac bc))
       (define c1 (+ (* 3.0 ac3) (* -15.0 acbc) (* 2.0 cc25)))
       (define a (- bc (* 0.4 ac2)))          ; y^3
       (define b (* (+ c1 ac3 (- cc25)) 0.04)) ; y^2
       (define c (- dc (* c1 ac 0.008)))       ; y^1
       (define d (+ (* ac (- (* ac (+ (* 4.0 ac3) (* -25.0 acbc) (* 5.0 cc25)))
                             (* 625.0 dc)) 0.00032) ec)) ; y^0
       (cond
         [(= 0. (abs d))
          (map (lambda ((r : Float)) (- r a_5d))
               (cons 0. (polyroots/quart 1. 0. a b c #:polish-roots? #f)))]
         [else
          ;; Calculate the *isolator polynomials* (a pair of quadratics)
          ;; cf. "Isolator Polynomials" - Sederberg T., Chang G.
          (define a2 (* a a))
          (define b2 (* b b))
          ;; First isolator polynomial A(x)
          (define siAx2 (+ (* 12. a2 a) (* 45. b2) (* -40. a c)))
          (define siAx (+ (* 8. a2 b) (* 60. b c) (* -50. a d)))
          (define siA1 (+ (* 4. a2 c) (* 75. b d)))
          ;; second isolator polynomial B(x)
          (define siBx2 (* 10. a))
          (define siBx (* -15. b))
          (define siB1 (* 4. a2))
          ;; Calculcate bounds of real roots of the original quintic
          ;; Lagrange's root bounds
          ;; Given P(x) = x^5 + ac x^4 + bc x^3 + cc x^2 + dc x + ec
          ;;    Upper bound on positive roots = 1 + (-B^(1/k))
          ;;        where, B is the largest negative coeff. and
          ;;               k is the number of leading positive coeff.
          ;; Minimum bound = -1 * upper-bound of P(-x)
          ;; Upper bound:
          (define B+ (min ac bc cc dc ec))
          (define bound-right
            (if (> B+ 0.)
                0.
                (let ([k+ (if (< ac 0.) 1.
                              (if (< bc 0.) 2.
                                  (if (< cc 0.) 3.
                                      (if (< dc 0.) 4.
                                          (if (< ec 0.) 5. 0.)))))])
                  (+ 1. (if (= 1. k+) (- B+) (flexpt (- B+) (/ 1. k+)))))))
          ;; Lower bound:
          (define ac- (- ac))
          (define cc- (- cc))
          (define ec- (- ec))
          (define B- (min ac- bc cc- dc ec-))
          (define bound-left
            (if (> B- 0.)
                0.
                (let ([k- (if (< ac- 0.) 1.
                              (if (< bc 0.) 2.
                                  (if (< cc- 0.) 3.
                                      (if (< dc 0.) 4.
                                          (if (< ec- 0.) 5. 0.)))))])
                  (- (+ 1. (if (= 1. k-) (- B-) (flexpt (- B-) (/ 1. k-))))))))
          (define rBounds `(,@(sort
                               (map (lambda ((r : Float)) (- r a_5d))
                                    (append (polyroots/quad siAx2 siAx siA1)
                                            (polyroots/quad siBx2 siBx siB1)))
                               fl<)
                            ,bound-right))
          ;; Find a suitable bracket in which we are certain the polynimial has a root
          (define-values (x1 y1 x2 y2)
            (let loop : (Values Float Float Float Float)
                 ([b bound-left] [bs rBounds])
                 (cond
                   [(empty? bs) (values +inf.0 +inf.0 +inf.0 +inf.0)]
                   [else
                    (let-values ([(bv be) (-eval-quint 1.0 ac bc cc dc ec b)]
                                 [(bsv bse) (-eval-quint 1.0 ac bc cc dc ec (car bs))])
                      (define at-zero? (or (<= (abs bv) (* 4. be))
                                           (<= (abs bsv) (* 4. bse))))
                      (if (or at-zero? (<= (* bv bsv) 0.))
                          (if at-zero?
                              (values (car bs) bsv (car bs) bsv)
                              (values b bv (car bs) bsv))
                          (loop (car bs) (cdr bs))))])))
          (if (infinite? x1)
              ;; All Quintic polynomials must have atleast 1 real root.
              ;; Normally this branch should never be taken.
              ;; although...
              empty
              (let ()
                ;; Numerical precision is still an issue in these cases
                (define xraw (if (<= (abs (- x1 x2)) epsilon.0)
                                 (min x1 x2) ; Already at a root
                                 (-bracketed-root/quint a0l b0l c0l d0l e0l f0l
                                                        x1 y1 x2 y2)))
                (define x (-root-polish/quint a0l b0l c0l d0l e0l f0l xraw))
                (define bq (+ b0l (* a0l x)))
                (define cq (+ c0l (* bq x)))
                (define dq (+ d0l (* cq x)))
                (define eq (+ e0l (* dq x)))
                (cons x (polyroots/quart a0l bq cq dq eq))))])]))
  ;; Remove duplicates and "polish the roots
  (define rsu (remove-duplicates rs
                                 (lambda ((a : Float) (b : Float))
                                   (<= (abs (- a b)) epsilon.0))))
  (for/list ([r (in-list rsu)])
    (-root-polish/quint a0l b0l c0l d0l e0l f0l r)))



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


(: test-profile/ffi (case->
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
(define (test-profile/ffi fn polys N)
  (define count (length polys))
  (define times  (inexact->exact (ceiling (/ N count))))
  (for* ([i (in-range 0 times)]
         [c (in-list polys)])
    (apply fn c)))


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
        (unless (for/and : Boolean ([r (in-list rs)]
                                    [a (in-list as)])
                  (<= (abs (- r a)) accuracy))
          (define ulp-errors (for/list : (Listof Float) ([r (in-list rs)]
                                                         [a (in-list as)])
                               (flulp-error r a)))
          (displayln (~a "Poly:    " ca))
          (displayln-digit-diff as rs #:prefix-a "Accuracy ")
          (displayln (~a "Ulp-err  " ulp-errors)))
        (displayln-digit-diff as rs #:prefix-a "Mismatch: Number of Roots: "))))




(let* ([times 100000]
       [count 10000]
       [quads/test (random-polys/quad 10)]
       [cubics/test (random-polys/cubic 10)]
       [quarts/test (random-polys/quart 10)]
       [quints/test (random-polys/quint 10)]
       [quads (random-polys/quad count)]
       [cubics (random-polys/cubic count)]
       [quarts (random-polys/quart count)]
       [quints (random-polys/quint count)]
       )
  (displayln "--------------------------------------------------")
  (displayln "Accuracy Tests:")

  (displayln "Quadratic FFI vs Racket:")
  (test-poly-results/ffi polyroots/quad poly_quadroots quads/test)

  (displayln "\nCubic FFI vs Racket/1:")
  (test-poly-results/ffi polyroots/cubic poly_cubicroots testdata/cubic)
  (displayln "Cubic FFI vs Racket/2:")
  (test-poly-results/ffi polyroots/cubic poly_cubicroots cubics/test)

  (displayln "\nQuartic FFI vs Racket/1:")
  (test-poly-results/ffi polyroots/quart poly_quartroots testdata/quart)
  (displayln "Quartic FFI vs Racket/2:")
  (test-poly-results/ffi polyroots/quart poly_quartroots quarts/test)

  (displayln "\nQuintic FFI vs Racket/1:")
  (test-poly-results/ffi polyroots/quint poly_quintroots testdata/quint)
  (displayln "Quintic FFI vs Racket/2:")
  (test-poly-results/ffi polyroots/quint poly_quintroots quints/test)

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
  (display "Typed/Racket/Quartic \t")
  (test-perf/ffi polyroots/quart quarts times)

  (collect-garbage)
  (display "FFI/Quintic          \t")
  (test-perf/ffi poly_quintroots quints times)
  (collect-garbage)
  (display "Typed/Racket/Quintic \t")
  (test-perf/ffi polyroots/quint quints times))
