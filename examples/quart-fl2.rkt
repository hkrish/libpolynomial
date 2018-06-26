#lang racket

(require (for-syntax racket/syntax
                     syntax/parse)
         math
         math/flonum
         (only-in ffi/unsafe ffi-lib _fun _ptr _double _long _int _int32 ->)
         ffi/unsafe/define)

(define-ffi-definer define-polys (ffi-lib "libpolynomial5"))
(define-ffi-definer define-libm (ffi-lib "libm"))

(define-polys poly_quadroots (_fun _double _double _double
                                   (_double = -inf.0)
                                   (_double = +inf.0)
                                   (x1 : (_ptr o _double))
                                   (x2 : (_ptr o _double))
                                   -> (r : _int32)
                                   -> (take (list x1 x2) r)))

(define-polys poly_cubicroots (_fun _double _double _double _double
                                    (_double = -inf.0)
                                    (_double = +inf.0)
                                    (x1 : (_ptr o _double))
                                    (x2 : (_ptr o _double))
                                    (x3 : (_ptr o _double))
                                    -> (r : _int32)
                                    -> (take (list x1 x2 x3) r)))

(define-polys poly_quartroots (_fun _double _double _double _double _double
                                    (_double = -inf.0)
                                    (_double = +inf.0)
                                    (x1 : (_ptr o _double))
                                    (x2 : (_ptr o _double))
                                    (x3 : (_ptr o _double))
                                    (x4 : (_ptr o _double))
                                    -> (r : _int32)
                                    -> (take (list x1 x2 x3 x4) r)))

(define-libm frexp (_fun _double (e : (_ptr o _int))
                         -> (m : _double) -> e))

(define-libm scalbn (_fun _double _int -> _double))


;;; Example of solving a Quartic using double-double data types vs doubles

;; Polynomial is
;;   3.70944181855249*10^-22 - 2.17240456536774*10^-21 x
;;     + 4.72995343184905*10^-21 x^2 - 4.542145068966119*10^-21 x^3
;;     + 1.624401562*10^-21 x^4
;;
;; The roots are (from Mathematica):
;;  x == 0.5399708576229223
;;  x == 0.6990489875207541
;;  x == 0.7435933712274639
;;  x == 0.8135827337078411
;;

(define tests
  '((1.0 -2.63453369092816  11.47594580628995 -3.01193978276912 0.17801932206896)
    (1.0 -2.63453369092816 2.47594580628995 -0.969679092568218 0.133302415104222)
    (16.0 -16.62910939950919 8.30968497920564 -1.86739743049595 0.12844638447658)
    (1.0  6  -5  -10  -3)
    (1.0  -2.35220164231774431  2.00357193581024839  -0.732009944503505094 0.0963567270541595411)
    (4213952.415  -6257401.32508918  2907561.77975298  -488990.675178038  20927.9939160697)
    (57341.571  -161727.801652109  168247.656267706  -76482.4511413003  12801.3771995616)
    (1.624401562e-21  -4.54214506896612e-21  4.72995343184905e-21  -2.17240456536774e-21 3.70944181855249e-22)
    (1.324862655e-21  -2.25090952994387e-21  1.2704566968795e-21  -2.55885000154948e-22 9.1991026101414e-24)
    (3852.061343  -8484.23316822249  6538.97688624379  -2056.39000615966  212.050414780893)
    (3340.933089  -4990.31721011812  3726.99260868093  -1391.74258071197  194.890453343193)
    (3294.612452  -5742.33193538764  3753.21565167956  -1090.27543254508  118.768207298396)
    (1.63064755140304e-16  -2034.51169316325  1215.54655528735  -322.775442718126 32.1410930237635)
    (1606.268505  -3167.86593395418  2156.86431456392  -586.683298630047  53.9570421132115)
    (3063.145275  -4196.79839872473  2875.00513659746  -290.139062231677  5.32413681477471)
    (1836.085081  -3157.51102285065  2036.23649359839  -378.014005550609  27.7813065153978)
    (1347.375714  -2372.21108226666  1566.20719086287  -161.089864668723  5.63065402793846)
    (463.822363  -787.053152298914  500.826971130479  -141.64095202131  1.79148644473404)))

(define answers
  '((0.08809808973990368 0.1869183562478567)
    (0.332355001781926 0.4785910970082226)
    (0.117153125076721 0.3087847625899946)
    (-0.61803398874989484820 1.6180339887498948482 -6.5413812651491098445 -0.45861873485089015550)
    (0.33995335944143 0.55372069775034 0.58805041057944 0.87047717454653)
    (0.06385398737126084 0.267956069261246 0.3712311334375205 0.7818833436801197)
    (0.4897764366689155 0.6938173557314637 0.7051071274833584 0.9317275900380739)
    (0.5399708576228983 0.6990489875212664 0.743593371226479 0.8135827337080589)
    (0.04537099624177689 0.4012388080361285 0.5227712126269068 0.7295947383925314)
    (0.2018746976609952 0.55062936130567615948 0.55062947295979712970 0.8993841366044615)
    (0.3734224150244674)
    (0.4357365258470498 0.4357365258470498)
    (0.1991545128284868)
    (0.2023857558722718 0.3160990356192021 0.6699957312560749 0.783709011003024)
    (0.02374856039328096 0.0925818978834563)
    ()
    ()
    (0.01325662945878696 0.8351858078659913)))

(define test-no 9)
(define coeffs (list-ref tests test-no))
(define rts (list-ref answers test-no))

(define (quadroot2 a B c)
  (define b (lv/ B -2.))
  (define D (lv- (lv* b b) (lv* a c)))
  (define roots (cond
                  [(lv< (lvabs a) epsilon.0)
                   (if (lv< (lvabs B) epsilon.0) '(0.) (list (fl/ (fl- 0. c) B)))]
                  [(lv> D (- epsilon.0))
                   (define Q (if (lv< D 0.) 0. (lvsqrt D)))
                   (define R (lv+ b (lv* (sgn (car b)) Q)))
                   (cond
                     [(lv= R 0.)
                      (list (lv/ c a) (lv/ (lv- 0. c) a))]
                     [else
                      (list (lv/ R a) (lv/ c R))])]
                  [else
                   '()]))
  (set->list (list->set roots)))

(match-define (list a b c d) `(-3. 3. -1. (/ 10. 90.))) ; DEBUG
(define (cubroot2 a0 b0 c0 d0)
  (define a (if (list? a0) a0 (list a0 0.)))
  (define b (if (list? b0) b0 (list b0 0.)))
  (define c (if (list? c0) c0 (list c0 0.)))
  (define d (if (list? d0) d0 (list d0 0.)))
  (cond
    [(lv= 0. a) (quadroot2 b c d)]
    [(lv= 0. d) `(0. @(quadroot2 a b c))]
    [else
     (define x (lv/ (lv/ b a) -3.))
     (define ax (lv* a x))
     (define b1 (lv+ ax b))
     (define c2 (lv+ (lv* b1 x) c))
     (define qd (lv+ (lv* (lv+ ax b1) x) c2))
     (define q (lv+ (lv* c2 x) d))
     (define t (lv/ q a))
     (define r `(,(flexpt+ (car t) (cadr t) (/ 1. 3.)) 0.))
     (define s (if (lv< t 0.) -1. 1.))
     (set! t (lv/ (lv- 0. qd) a))
     (set! r (if (lv> t 0.) (lv* 1.324717957244746025 (lvmax r (lvsqrt t))) r))
     (define x0 (lv- x (lv* s r)))
     (if (not (lv= x0 x))
         (let ()
           (define (NR-rec x)
             (define ax (lv* a x))
             (define b1 (lv+ ax b))
             (define c2 (lv+ (lv* b1 x) c))
             (define qd (lv+ (lv* (lv+ ax b1) x) c2))
             (define q (lv+ (lv* c2 x) d))
             (define xn (if (lv= qd 0.) x (lv- x (lv/ q qd))))
             (if (> (car (lv* xn s)) (car (lv* x s)))
                 (NR-rec xn)
                 (values x b1 c2)))
           (define-values (x b1 c2) (NR-rec x0))
           (if (lv> (lv* (lv* (lvabs a) x) x) (lvabs (lv/ d x)))
               (let ([c2 (lv/ (lv- 0. d) x)]
                     [b1 (lv/ (lv- c2 c) x)])
                 `(x @(quadroot2 a b1 c2)))
               `(,x ,@(quadroot2 a b1 c2))))
         `(,x ,@(quadroot2 a b1 c2)))]))

(define (scale-coeffs cs)
  (define s (apply max (map abs cs)))
  (define p (+ 1 (- (frexp s))))
  (map (lambda (c) (scalbn c p)) cs))

(define-syntax-rule (vsl form) (let-values ([(a b) form]) (list a b)))
(define-syntax-rule (lvs form) (apply values form))

(define-syntax-rule (lv1 fn a) (match-let ([(list x1 x2) a])
                                 (vsl (fn x1 x2))))

(define-syntax-rule (lv2 fn a b) (match-let ([(list x1 x2) a]
                                             [(list y1 y2) b])
                                   (vsl (fn x1 x2 y1 y2))))

(define-syntax (define-fl2-fn stx)
  (syntax-case stx ()
    [(_ op a b res)
     (with-syntax ([fl2op (format-id stx "fl2~a" #'op)]
                   [name (format-id stx "lv~a" #'op)])
       #'(define (name a b)
           (match-let ([(list x1 x2) (if (list? a) a (list a 0.))]
                       [(list y1 y2) (if (list? b) b (list b 0.))])
             (res (fl2op x1 x2 y1 y2)))))]
    [(_ op a res)
     (with-syntax ([fl2op (format-id stx "fl2~a" #'op)]
                   [name (format-id stx "lv~a" #'op)])
       #'(define (name a)
           (match-let ([(list x1 x2) (if (list? a) a (list a 0.))])
             (res (fl2op x1 x2)))))]))

(define-fl2-fn = a b identity)
(define-fl2-fn > a b identity)
(define-fl2-fn < a b identity)
(define-fl2-fn >= a b identity)
(define-fl2-fn <= a b identity)
(define-fl2-fn ulp-error a b identity)

(define-fl2-fn + a b vsl)
(define-fl2-fn - a b vsl)
(define-fl2-fn * a b vsl)
(define-fl2-fn / a b vsl)
(define-fl2-fn sqr a vsl)
(define-fl2-fn sqrt a vsl)
(define-fl2-fn abs a vsl)
(define-fl2-fn prev a vsl)
(define-fl2-fn next a vsl)

(define (lvmax a b) (if (lv>= a b) a b))

(define (lvmax* . rs) (for/fold ([mx (car rs)])
                                ([r (in-list (cdr rs))])
                        (lvmax mx r)))

(define-syntax-rule (lv-re a) (match-let ([(list x1 x2) a]) (fl2->real x1 x2)))
(define-syntax-rule (lv-r a) (exact->inexact (match-let ([(list x1 x2) a]) (fl2->real x1 x2))))


(define (quartroot2 coeffs)
  (define scs (scale-coeffs coeffs))
  (match-define (list _ a b c d) (map (lambda (c) (vsl (fl//error c (car scs)))) scs))

  (define a2 (lvsqr a))
  (define a3 (lv* a2 a))
  (define a4 (lv* a3 a))
  (define a_4d (lv/ a 4.))

  (define e (lv- b (lv/ (lv* 3. a2) 8.)))
  (define f (lv+ c (lv- (lv/ a3 8.) (lv/ (lv* a b) 2.))))
  (define g (lv- (lv+ (lv- d (lv/ (lv* 3. a4) 256.)) (lv/ (lv* a2 b) 16.)) (lv/ (lv* a c) 4.)))

  (cond
    [(lv<= (lvabs g) epsilon.0)
     (define rs
       (cons (lv- 0. a_4d)
             (if (lv<= (lvabs f) epsilon.0)
                 (if (lv< e 0.)
                     (let ([qr1 (lvsqrt (lv- 0. e))])
                       (list (lv- qr1 a_4d) (lv- (lv- 0. qr1 a_4d))))
                     '())
                 (let ([ysc (cubroot2 1.0 0.0 e f)])
                   (map (lambda (r) (lv- r a_4d)) ysc)))))
     (sort rs < #:key car)]
    [(lv<= (lvabs f) epsilon.0)
     (void)]
    [else
     (define b1 (lv* 2. e))
     (define c1 (lv- (lv* e e) (lv* 4. g)))
     (define d1 (lv* (lv- 0. f) f))
     (define ysc (cubroot2 1. b1 c1 d1))
     (define h2 (apply lvmax* ysc))
     (define h (lvsqrt h2))
     (define j (lv/ (lv- (lv+ e h2) (lv/ f h)) 2.))
     (define qr1s (quadroot2 1. h j))
     (define r1s (map (lambda (r) (lv- r a_4d)) qr1s))
     (define qr2s (quadroot2 1. (lv- 0. h) (lv/ g j)))
     (define r2s (map (lambda (r) (lv- r a_4d)) qr2s))
     (define rs (map ((curry polish) coeffs 1) (append r1s r2s)))
     (sort rs < #:key car)]))


(define roots (apply poly_quartroots coeffs))
;; (displayln roots)

(define (abs-errors a b)
  (for/list ([r (in-list a)]
             [ro (in-list b)])
    (if (list? r)
        (lv- r ro)
        (- r ro))))

(define eval-c
  (curry (lambda (coeffs x)
           (for/fold ([v '(0. 0.)])
                     ([c (in-list coeffs)])
             (lv+ (lv* v x) c)))))

(define eval-e
  (curry (lambda (coeffs x)
           (define ex (if (list? x) (lv-re x) (inexact->exact x)))
           (exact->inexact (for/fold ([v 0])
                                     ([c (in-list (map inexact->exact coeffs))])
                             (+ (* v ex) c))))))

(define (print-coeffs cs) (displayln (string-join (map number->string cs) ", ")))

(define (polyD cs)
  (define d (sub1 (length cs)))
  (for/list ([i (in-range d 0 -1)]
             [c (in-list (drop-right cs 1))])
    (* i c)))

(define (polish coeffs n r)
  (if (= n 0)
      r
      (let ([qd (eval-c (polyD coeffs) r)])
        (if (lv= 0. qd)
            r
            (polish coeffs (sub1 n) (lv- r (lv/ (eval-c coeffs r) qd)))))))

(define rs (quartroot2 coeffs))
(define crs (apply poly_quartroots coeffs))
;; (abs-errors roots rts)

rts
crs
rs
