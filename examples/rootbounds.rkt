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


(begin-encourage-inline

  (: finite? (-> Real Boolean))
  (define (finite? n) (or (nan? n) (infinite? n)))


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


    (: polish-rec (-> Float Float Integer Float))
    (define (polish-rec x z rec)
      (define dz (+ e (* x (+ dd (* x (+ dc (* x (+ db (* x da)))))))))
      (define nx (if (= 0. dz) x (- x (/ z dz))))
      (define nz (+ f (* nx (+ e (* nx (+ d (* nx (+ c (* nx (+ b (* nx a)))))))))))
      (if (and (not (= nx x)) (< (abs nz) (abs z)) (> (abs nz) epsilon.0) (> rec 0))
          (polish-rec nx nz (sub1 rec))
          (if (< (abs nz) (abs z)) nx x)))
    (polish-rec x0 z max-rec)))


;; Refines a root (root search) in an interval known to contain exactly one root
;; of a quintic polynomial.
;;
;; cf. "Bisected Direct Quadratic Regula Falsi"
;; - Robert G. Gottlieb and Blair F. Thompson
;;
(: -bracket-root-quint (-> Float Float Float Float Float Float Float Float Float Float Float))
(define (-bracket-root-quint a b c d e f x1 x2 y1 y2)
  (define x-last0 : Float (if (< y1 y2) x1 x2))
  (define-values (xup0 xdn0 yup0 ydn0) (if (< y1 0.)
                                           (values x2 x1 y2 y1)
                                           (values x1 x2 y1 y2)))
  ;; (define n-iter 0)
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
      #:break (or (= y 0) (<= (abs (- xdn xup)) ε))
      (define delta  (/ (- xup xdn) 2.))
      (define xm     (/ (+ xup xdn) 2.))
      (define ym     (+ f (* xm (+ e (* xm (+ d (* xm (+ c (* xm (+ b (* xm a)))))))))))
      (define a_d    (* 2. delta delta))
      (define a_     (/ (+ yup ydn (* -2. ym)) a_d))
      (define b_     (/ (- yup ydn) (* 2. delta)))
      (define xden   (flsqrt (- 1. (* (/ (* 4. a_) b_) (/ ym b_)))))
      (define x      (- xm (/ (* 2. ym) (* b_ (+ 1. xden)))))
      #:break (= x x-last)
      ;; (set! n-iter (add1 i))
      (define yn     (+ f (* x (+ e (* x (+ d (* x (+ c (* x (+ b (* x a)))))))))))
      (define-values (xupn xdnn yupn ydnn)
        (cond
          [(and (> yn 0.) (< ym 0.)) (values x xm yn ym)]
          [(> yn 0.)                 (values x xdn yn ydn)]
          [(> ym 0.)                (values xm x ym yn)]
          [else                     (values xup x yup yn)]))
      (values x yn xupn xdnn yupn ydnn)))
  ;; (displayln (format "Interations (~a, ~a): ~a" x1 x2 n-iter))
  root)


(: flbracketed-root2 (-> Float Float Float Float Float Float Flonum Flonum Flonum Flonum Flonum))
(define (flbracketed-root2 an bn cn dn en fn a fa b fb)
  (let loop ([bisected? #t] [a a] [fa fa] [b b] [fb fb] [c a] [fc fa] [d 0.0] [n 0])
    (define min-abs-ab (min (abs a) (abs b)))
    (define ε (if (min-abs-ab . <= . +max-subnormal.0) +min.0 (* min-abs-ab epsilon.0)))
    (cond
      ;; If we got it right, return it
      [(= fb 0.0)  b]
      ;; If a and b are too close, return b
      [((abs (- a b)) . <= . ε) b]
      [(n . < . 500)
       (let*-values
           ([(bisect? s)
             (cond
               ;; Weird rules for forcing bisection to make progress
               [(or (and bisected?       ((abs (- b c)) . < . (* 128.0 ε)))
                    (and (not bisected?) ((abs (- c d)) . < . (* 128.0 ε))))
                (values #t 0.0)]
               ;; Get an interpolated point
               [else
                (define fa-fb (- fa fb))
                (define fb-fc (- fb fc))
                (define fc-fa (- fc fa))
                (define da (* fa-fb (- fc-fa)))
                (define db (* fb-fc (- fa-fb)))
                (define dc (* fc-fa (- fb-fc)))
                (cond [(or (= da 0.0) (= db 0.0) (= dc 0.0))
                       ;; Secant method because quadratic method will fail
                       (values #f (- b (* fb (/ (- b a) (- fa-fb)))))]
                      [else
                       ;; Inverse quadratic interpolation method
                       (values #f (+ (/ (* a fb fc) da)
                                     (/ (* b fc fa) db)
                                     (/ (* c fa fb) dc)))])])]
            ;; Run tests to determine whether it's a bad point
            [(bisected? s)
             (cond
               [(or bisect?
                    (not (let ([s0  (/ (+ (* 3.0 a) b) 4.0)])
                           (if (<= s0 b) (<= s0 s b) (<= b s s0))))
                    (and bisected?       ((abs (- s b)) . >= . (* 0.5 (abs (- b c)))))
                    (and (not bisected?) ((abs (- s b)) . >= . (* 0.5 (abs (- c d))))))
                ;; Bisect this time
                (values #t (* 0.5 (+ a b)))]
               [else
                (values #f s)])]
            [(d)  c]
            [(c fc)  (values b fb)]
            [(fs)  (+ fn (* s (+ en (* s (+ dn (* s (+ cn (* s (+ bn (* s an))))))))))]
            ;; Replace the endpoint with the same sign as s
            [(a fa b fb)  (if ((* fa fs) . < . 0.0)
                              (values a fa s fs)
                              (values s fs b fb))]
            ;; Make sure b is closer
            [(a fa b fb)  (if ((abs fa) . < . (abs fb))
                              (values b fb a fa)
                              (values a fa b fb))])
         (loop bisected? a fa b fb c fc d (+ n 1)))]
      [else b])))


(: poly-normalize (-> Float Float Float Float Float Float
                      (Values Float Float Float Float Float Float)))
(define (poly-normalize a b c d e f)
  (define infnorm (max (abs a) (abs b) (abs c) (abs d) (abs f)))
  (let* ([scale (flexp2 (- (flround (fllog2 (if (= 0. infnorm) epsilon.0 infnorm)))))]
         [an (* a scale)]
         [bn (* b scale)]
         [cn (* c scale)]
         [dn (* d scale)]
         [en (* e scale)]
         [fn (* f scale)])
    (values an bn cn dn en fn)))


(define-type BoundsFn (-> Float Float Float Float Float Float (Pairof Float Float)))

(: bounds/laguerre BoundsFn)
(define (bounds/laguerre a b c d e f)
  (define ac (/ b a))
  (define bc (/ c a))
  (define cc (/ d a))
  (define dc (/ e a))
  (define ec (/ f a))
  (define rmean (/ (- ac) 5.))
  (define rstdd (* (/ 4. 5.) (sqrt (abs (- (* ac ac) (* 2.5 bc))))))
  (define bound-left (- rmean rstdd))
  (define bound-right (+ rmean rstdd))
  (cons bound-left bound-right))

(: bounds/lagrange BoundsFn)
(define (bounds/lagrange a b c d e f)
  (define-values (an bn cn dn en fn) (poly-normalize a b c d e f))
  (define absan (abs an))
  (define absfn (abs fn))
  (define boundu (+ 1. (max (/ (abs bn) absan) (/ (abs cn) absan) (/ (abs dn) absan)
                            (/ (abs en) absan) (/ absfn absan))))
  (cons (- boundu) boundu))


(begin-encourage-inline
  (: bounds/cauchy-helper (-> (List Float Float Float Float Float) Float))
  (define (bounds/cauchy-helper cs)
    (define n- (- (fl (for/sum : Real ([c (in-list cs)]) (if (< c 0.) 1. 0.)))))
    (for/fold : Float
        ([bu 0.])
        ([c (in-list cs)]
         [k (in-naturals 1)])
      (if (< c 0.)
          (max bu (flexpt (* n- c) (/ 1. (fl k))))
          bu))))

(: bounds/cauchy BoundsFn)
(define (bounds/cauchy a b c d e f)
  (define-values (an bn cn dn en fn) (poly-normalize a b c d e f))
  (define cs (list (/ bn an) (/ cn an) (/ dn an) (/ en an) (/ fn an)))
  (define cs- (list (/ (- bn) an) (/ cn an) (/ (- dn) an) (/ en an) (/ (- fn) an)))
  (cons (- (bounds/cauchy-helper cs-)) (bounds/cauchy-helper cs)))


(: bounds/lagrange2 BoundsFn)
(define (bounds/lagrange2 an bn cn dn en fn)
  (define cs+ (list (/ bn an) (/ cn an) (/ dn an) (/ en an) (/ fn an)))
  (define cs- (list (/ (- bn) an) (/ cn an) (/ (- dn) an) (/ en an) (/ (- fn) an)))
  (define-values (k+ B+ k- B-)
    (for/fold : (Values Float Float Float Float)
        ([k+ 0.]
         [B+ 0.]
         [k- 0.]
         [B- 0.])
        ([i (in-range 1. 6.)]
         [c+ (in-list cs+)]
         [c- (in-list cs-)])
      (define-values (kn+ Bn+)
        (if (< c+ 0.)
            (values (if (= 0. k+) i k+) (min B+ c+))
            (values k+ B+)))
      (define-values (kn- Bn-)
        (if (< c- 0.)
            (values (if (= 0. k-) i k-) (min B- c-))
            (values k- B-)))
      (values kn+ Bn+ kn- Bn-)))
  (define bu (+ 1. (if (= 1. k+) (- B+) (flexpt (- B+) (/ 1. k+)))))
  (define bl (- (+ 1. (if (= 1. k-) (- B-) (flexpt (- B-) (/ 1. k-))))))
  (cons bl bu))


(: bounds/lagrange2-impl2 BoundsFn)
(define (bounds/lagrange2-impl2 an bn cn dn en fn)
  (define ac (/ bn an))
  (define bc (/ cn an))
  (define cc (/ dn an))
  (define dc (/ en an))
  (define ec (/ fn an))
  (define ac- (- ac))
  (define cc- (- cc))
  (define ec- (- ec))
  (define B+ (min ac bc cc dc ec))
  (define bu
    (if (> B+ 0.)
        0.
        (let ([k+ (if (< ac 0.) 1.
                      (if (< bc 0.) 2.
                          (if (< cc 0.) 3.
                              (if (< dc 0.) 4.
                                  (if (< ec 0.) 5. 0.)))))])
          (+ 1. (if (= 1. k+) (- B+) (flexpt (- B+) (/ 1. k+)))))))
  (define B- (min ac- bc cc- dc ec-))
  (define bl
    (if (> B- 0.)
        0.
        (let ([k- (if (< ac- 0.) 1.
                      (if (< bc 0.) 2.
                          (if (< cc- 0.) 3.
                              (if (< dc 0.) 4.
                                  (if (< ec- 0.) 5. 0.)))))])
          (- (+ 1. (if (= 1. k-) (- B-) (flexpt (- B-) (/ 1. k-))))))))
  (cons bl bu))


(: bracket-root (-> (List Float Float Float Float Float Float) Float Float Float))
(define (bracket-root ps lb ub)
  (match-define (list a b c d e f) ps)
  (define-values (ylb ylbe) (-eval-quint a b c d e f lb))
  (define-values (yub yube) (-eval-quint a b c d e f ub))
  (-bracket-root-quint a b c d e f lb ub ylb yub))

(: bracket-root2 (-> (List Float Float Float Float Float Float) Float Float Float))
(define (bracket-root2 ps lb ub)
  (match-define (list a b c d e f) ps)
  (define-values (ylb ylbe) (-eval-quint a b c d e f lb))
  (define-values (yub yube) (-eval-quint a b c d e f ub))
  (flbracketed-root2 a b c d e f lb ylb ub yub))

;; (define quints (random-polys/quint 100))

;; (displayln quints)

;; (test-perf/ffi polyroots/quint quints 1000)

;; (test-perf/ffi poly_quintroots quints 1000)




(define dqs '((3198.3652844233393 -5263.891082877705 3473.96514058461 -983.2723253353927 118.93285630171899 3198.3645221700654)
              ;; (-0.6787702423448035)
              (6167.571999797359 2234.5984792014215 1584.3439161922431 -129.9479962383686 3.459235351545105 6167.57199538209)
              ;; (-1.018512639063058)
              (1898.2588450512478 1119.9762460009108 467.0698372317855 -12.212462634641785 1.1154481782204775 1898.2588411445997)
              ;; (-1.083224927899759)
              (22228992000000.0 -55572480000000.0 57117707539200.0 -30104081308800.0 7961050818000.0 -815594524200.0)
              ;; (0.269731766859053 0.7302682331409409)
              (1. -10. 15. 4. -16. 400.)
              ;; (-2.2123376379624875742 7.9945186684323459459)
              (1. -0.0037841796875 4.61935997009277e-6 -2.25554686039686e-9 4.40536496171262e-13 -2.77555756156289e-17)
              ))


(displayln "bounds/laguerre  (NOTE: Only works if all roots are real!)")
(for/list : (Listof (Pairof Float Float))
    ([q (in-list dqs)])
  (apply bounds/laguerre q))

(displayln "bounds/lagrange")
(for/list : (Listof (Pairof Float Float))
    ([q (in-list dqs)])
  (apply bounds/lagrange q))

(displayln "bounds/cauchy")
(for/list : (Listof (Pairof Float Float))
    ([q (in-list dqs)])
  (apply bounds/cauchy q))


(displayln "bounds/lagrange2")
(for/list : (Listof (Pairof Float Float))
    ([q (in-list dqs)])
  (apply bounds/lagrange2 q))

(displayln "bounds/lagrange2-impl2")
(for/list : (Listof (Pairof Float Float))
    ([q (in-list dqs)])
  (apply bounds/lagrange2-impl2 q))


;; (bracket-root (list-ref dqs 3) -3. 0.3)
(define psrg (current-pseudo-random-generator))

(: time-bracket (-> (List Float Float Float Float Float Float) Integer Void))
(define (time-bracket cs n)
  (match-define (list a b c d e f) cs)
  (time (for ([i (in-range 0 n)])
          ;; (define x1 (- -3.0 (* 400 (random psrg))))
          (define x1 -7.3)
          (define y1 (+ f (* x1 (+ e (* x1 (+ d (* x1 (+ c (* x1 (+ b (* x1 a)))))))))))
          (define x2 1.)
          (define y2 (+ f (* x2 (+ e (* x2 (+ d (* x2 (+ c (* x2 (+ b (* x2 a)))))))))))
          (-bracket-root-quint a b c d e f x1 x2 y1 y2))))

(: time-bounds-bracket (->  BoundsFn (List Float Float Float Float Float Float) Integer Void))
(define (time-bounds-bracket bfn cs n)
  (match-define (list a b c d e f) cs)
  (time (for ([i (in-range 0 n)])
          (define bs (bfn a b c d e f))
          (define x1 (car bs))
          (define y1 (+ f (* x1 (+ e (* x1 (+ d (* x1 (+ c (* x1 (+ b (* x1 a)))))))))))
          (define x2 0.0002)
          (define y2 (+ f (* x2 (+ e (* x2 (+ d (* x2 (+ c (* x2 (+ b (* x2 a)))))))))))
          (-bracket-root-quint a b c d e f x1 x2 y1 y2))))

(: time-bounds-flbracketed (->  BoundsFn (List Float Float Float Float Float Float) Integer Void))
(define (time-bounds-flbracketed bfn cs n)
  (match-define (list a b c d e f) cs)
  (time (for ([i (in-range 0 n)])
          (define bs (bfn a b c d e f))
          (define x1 (car bs))
          (define y1 (+ f (* x1 (+ e (* x1 (+ d (* x1 (+ c (* x1 (+ b (* x1 a)))))))))))
          (define x2 0.0002)
          (define y2 (+ f (* x2 (+ e (* x2 (+ d (* x2 (+ c (* x2 (+ b (* x2 a)))))))))))
          (flbracketed-root2 a b c d e f x1 y1 x2 y2))))

(: profile-bounds-bracket (->  BoundsFn (List Float Float Float Float Float Float) Integer Any))
(define (profile-bounds-bracket bfn cs n)
  (match-define (list a b c d e f) cs)
  (profile-thunk
   (thunk
    (for ([i (in-range 0 n)])
      (define bs (bfn a b c d e f))
      (define x1 (car bs))
      (define y1 (+ f (* x1 (+ e (* x1 (+ d (* x1 (+ c (* x1 (+ b (* x1 a)))))))))))
      (define x2 1.)
      (define y2 (+ f (* x2 (+ e (* x2 (+ d (* x2 (+ c (* x2 (+ b (* x2 a)))))))))))
      (-bracket-root-quint a b c d e f x1 x2 y1 y2)))))



(collect-garbage)
(display "Perf - bounds/laguerre        \t:")
(time-bounds-bracket bounds/laguerre (list-ref dqs 5) 100000)

(collect-garbage)
(display "Perf - bounds/lagrange        \t:")
(time-bounds-bracket bounds/lagrange (list-ref dqs 5) 100000)

(collect-garbage)
(display "Perf - bounds/cauchy          \t:")
(time-bounds-bracket bounds/cauchy (list-ref dqs 5) 100000)

(collect-garbage)
(display "Perf - bounds/lagrange2       \t:")
(time-bounds-bracket bounds/lagrange (list-ref dqs 5) 100000)

(collect-garbage)
(display "Perf - bounds/lagrange2-impl2 \t:")
(time-bounds-bracket bounds/lagrange2-impl2 (list-ref dqs 5) 100000)

(collect-garbage)
(time-bounds-flbracketed bounds/lagrange2-impl2 (list-ref dqs 5) 100000)



