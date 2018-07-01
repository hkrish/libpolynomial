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

;; Returns 3 values :  digits agreeing, differing in A, and differing in B
(define (digit-diff a b)
  (define astr (number->string a))
  (define bstr (number->string b))
  (define ds (add1 (for/last ([i (in-naturals 0)]
                              [c1 (in-string astr)]
                              [c2 (in-string bstr)]
                              #:break (not (char=? c1 c2)))
                     i)))
  (values (substring astr 0 ds) (substring astr ds) (substring bstr ds)))


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
        (displayln (format "~a\n~a" da db)))))
