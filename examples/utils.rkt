#lang racket/base

(require racket/format
         racket/function
         racket/string)

(provide (all-defined-out))


(define (random-polys/quad N)
  (define cprg (current-pseudo-random-generator))
  (define scale 3200.0)
  (define count (min N 1000))
  (for/list ([i (in-range 0 count)])
    (define a (* scale (random cprg)))
    (define b (random cprg))
    (define A (* scale (random cprg)))
    (define B (* scale (random cprg)))
    (define ca (* A B))
    (define cb (- (+ (* A b) (* a B))))
    (define cc (* a b))
    (list ca cb cc)))


(define (random-polys/cubic N)
  (define cprg (current-pseudo-random-generator))
  (define scale 3200.0)
  (define count (min N 1000))
  (for/list ([i (in-range 0 count)])
    (define a (* (if (> 0.5 (random cprg)) scale (- scale)) (random cprg)))
    (define b (random cprg))
    (define c (random cprg))
    (define cb (- a b c))
    (define cc (+ (* a b) (* a c) (* b c)))
    (define cd (* (- a) b c))
    (list 1.0 cb cc cd)))


(define (random-polys/quart N)
  (define cprg (current-pseudo-random-generator))
  (define scale 32000.0)
  (define count (min N 1000))
  (for/list ([i (in-range 0 count)])
    (define a (random cprg))
    (define b (random cprg))
    (define c (random cprg))
    (define d (random cprg))
    (define A (* scale (random cprg)))
    (define c1 (* A (- a b c d)))
    (define c2 (* A (+ (* a b) (* a c) (* b c) (* a d) (* b d) (* c d))))
    (define c3 (* (- A) (+ (* a b c) (* a b d) (* a c d) (* b c d))))
    (define c4 (* (- A) a b c d))
    (list A c1 c2 c3 c4)))


(define (random-polys/quint N)
  (define cprg (current-pseudo-random-generator))
  (define scale 32000.0)
  (define count (min N 1000))
  (for/list ([i (in-range 0 count)])
    (define a (random cprg))
    (define b (random cprg))
    (define c (random cprg))
    (define d (random cprg))
    (define e (random cprg))
    (define A (* scale (random cprg)))
    (define ab (* a b))
    (define cd (* c d))
    (define fc (- A (* ab cd e)))
    (define ec (* A (+ (* ab cd) (* ab c e) (* ab d e) (* a cd e) (* b cd e))))
    (define dc (* A (- (* ab c) (* ab d) (* a cd) (* b cd) (* ab e)
                       (* a c e) (* b c e) (* a d e) (* b d e) (* cd e))))
    (define cc (* A (+ ab (* a c) (* b c) (* a d) (* b d) (* d c) (* a e)
                       (* b e) (* c e) (* d e))))
    (define bc (* A (- a b c d e)))
    (list A bc cc dc ec fc)))


;; (define testdata/cubic
;;   '(((0. 0. 0. 0.) ())
;;     ((-1. 0. 0. -1.) (-1.0))
;;     ((1. 3. 3. 1.) (-1.0))
;;     ((1. 0. -2. 5.) (-2.094551481542327))
;;     ((1. -3. 2. 0.) (2.0 1.0 0.0))
;;     ((1. -3. 2. 2.34e-89) (2.0 1.0 -1.17e-89))
;;     ((1. -7999999999. -8000000002. 16000000000.) (-2.0 1.0 8000000000.0))
;;     ((16000000000. -8000000002. -7999999999. 1.) (-0.5 1.25e-10 1.0))
;;     ((1. -99999.00001 -99999.00001 1.) (-1.0 1e-05 100000.0))
;;     ((0.01 -300. 2990000. -299.) (0.000100000001003344))
;;     ((-3. 3. -1. 0.1111111111) (0.3331786137701938))
;;     ((10000000000. -9999999998. -10000000000. 9999999998.) (-1.0 1.0 0.9999999998))
;;     ((9999999998. -10000000000. -9999999998. 10000000000.) (-1.0 1.0000000002 1.0))))


(define testdata/cubic
  '((0. 0. 0. 0.)
    (-1. 0. 0. -1.)
    (1. 3. 3. 1.)
    (1. 0. -2. 5.)
    (1. -3. 2. 0.)
    (1. -3. 2. 2.34e-89)
    (1. -7999999999. -8000000002. 16000000000.)
    (16000000000. -8000000002. -7999999999. 1.)
    (1. -99999.00001 -99999.00001 1.)
    (0.01 -300. 2990000. -299.)
    (-3. 3. -1. 0.1111111111)
    (10000000000. -9999999998. -10000000000. 9999999998.)
    (9999999998. -10000000000. -9999999998. 10000000000.)))


(define testdata/quart
  '((1.0 -2.63453369092816  11.47594580628995 -3.01193978276912 0.17801932206896)
    (1.0 -2.63453369092816 2.47594580628995 -0.969679092568218 0.133302415104222)
    (16.0 -16.62910939950919 8.30968497920564 -1.86739743049595 0.12844638447658)
    (1.0  6.  -5.  -10.  -3.)
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


(define testdata/quint
  '((1.0 -2.52910939950919052 2.30968497920563746 -0.896739743049594606 0.128446384476579654 -0.00275665806066778772)
    (1.0 -2.52910939950919052 2.30968497920563746 -0.937101799553697019 0.169404234782936457 -0.0110787163926311995)
    (1.0 -2.63453369092816168 2.47594580628994667 -0.969679092568217881 0.133302415104222411 -0.00159027602636815394)
    (1.0 -2.63453369092816168 2.47594580628994667 -1.01193978276912147 0.178019322068959658 -0.0112575475542043011)
    ;; The quintic from the paper "Isolator Polynomials -- Segerberg et. al."
    (1.0 0. -23. 6. 112. -96.)
    (1. -10. 15. 4. -16. 400.)
    (1. -3.5 4.875 -3.375 1.16015625 -0.158203125)
    (1. 3. -5. -15. 4. 12.)
    (1. -1.5 -15.29 23.895 -11.36 1.68)
    (1. -0.0037841796875 4.61935997009277e-6 -2.25554686039686e-9 4.40536496171262e-13 -2.77555756156289e-17)
    (1. -1.921875 1.02197265625 -0.095947265625 0.00285911560058594 -2.50339508056641e-5)
    (1.0 -2.52910939950919052 2.30968497920563746 -0.896739743049594606 0.128446384476579654 -0.00095665806066778772)
    (1.0 -2.5291093995091905 2.3096849792056375 -0.89673974304959461 0.12844638447657965 -0.0069566580606677877)
    (1. 1.470890600490809 0.6165349007762284 0.05472411292864265 -0.003358094522770555 -0.003797486716822673)
    (796.692571 -9.27278946051005 183.852120446035 -936.213438348927 1630927.19368701 -2068.9112857205)
    (1920000. -4800000. 4400802. -1801203. 307951.125 -13775.0625)
    (22228992000000.0 -55572480000000.0 57117707539200.0 -30104081308800.0 7961050818000.0 -815594524200.0)))

;;; FIX: DEBUG: Quint: (9) (16)

;; --------------------------------------------------
;; Printing differences in digits
;;

(define (digit-diff->strings a b)
  (define astr (if (number? a) (number->string a) ""))
  (define bstr (if (number? b) (number->string b) ""))
  (define-values (cab da db)
    (cond
      [(and a b)
       (define dsl (for/last ([i (in-naturals 0)]
                              [c1 (in-string astr)]
                              [c2 (in-string bstr)]
                              #:break (not (char=? c1 c2)))
                     i))
       (define ds (if dsl (add1 dsl) 0))
       (values (substring astr 0 ds) (substring astr ds) (substring bstr ds))]
      [(not a) (values "" "" bstr)]
      [(not b) (values "" astr "")]))
  (define red (~a (integer->char #x1b) "[1;31m"))
  (define clear (~a (integer->char #x1b) "[0m"))
  (define wid (+ (max (string-length astr) (string-length bstr))
                 (string-length red) (string-length clear)))
  (values (~a cab red da clear #:width wid)
          (~a cab red db clear #:width wid)))


(define (displayln-digit-diff a b #:prefix-a (prea #f) #:prefix-b (preb #f))
  (define-values (prefa prefb)
    (cond
      [(and prea preb) (values prea preb)]
      [prea (values prea (make-string (string-length prea) #\ ))]
      [preb (values (make-string (string-length preb) #\ ) preb)]
      [else (values "" "")]))
  (if (and (list? a) (list? b))
      (let* ([wa (length a)]
             [wb (length b)]
             [la (if (> wa wb) a (append a (build-list (- wb wa) (const #f))))]
             [lb (if (> wb wa) b (append b (build-list (- wa wb) (const #f))))])
        (define-values (sa sb) (for/fold ([stra '()]
                                          [strb '()])
                                         ([ad (in-list la)]
                                          [bd (in-list lb)])
                                 (define-values (ads bds) (digit-diff->strings ad bd))
                                 (values (cons ads stra) (cons bds strb))))
        (displayln (format "~a(~a)\n~a(~a)"
                           prefa (string-join (reverse sa) "  ")
                           prefb (string-join (reverse sb) "  "))))
      (let-values ([(da db) (digit-diff->strings a b)])
        (displayln (format "~a~a\n~a~a" prefa da prefb db)))))
