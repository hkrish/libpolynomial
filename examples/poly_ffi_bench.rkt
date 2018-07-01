#lang racket

(require (for-syntax racket/syntax
                     syntax/parse)
         math
         math/flonum
         "defs.rkt"
         "utils.rkt")

(define (test-perf/quad N)
  (define polys (random-polys/quad N))
  (define count (length polys))
  (define (test-rec polys N cnt)
    (if (<= N 0)
        cnt
        (let ([roots (for/list ([c (in-list polys)])
                       (apply poly_quadroots c))])
          (test-rec polys (- N count)
                    (+ cnt (caadr roots))
                    ))))
  (time (test-rec polys N 0)))

(define (test-perf/cubic N)
  (define polys (random-polys/cubic N))
  (define count (length polys))
  (define (test-rec polys N cnt)
    (if (<= N 0)
        cnt
        (let ([roots (for/list ([c (in-list polys)])
                       (apply poly_cubicroots c))])
          (test-rec polys (- N count) (+ cnt (caadr roots))))))
  (time (test-rec polys N 0)))

(define (test-perf/quart N)
  (define polys (random-polys/quart N))
  (define count (length polys))
  (define (test-rec polys N cnt)
    (if (<= N 0)
        cnt
        (let ([roots (for/list ([c (in-list polys)])
                       (apply poly_quartroots c))])
          (test-rec polys (- N count) (+ cnt (caadr roots))))))
  (time (test-rec polys N 0)))

(define (test-perf/quint N)
  (define polys (random-polys/quint N))
  (define count (length polys))
  (define (test-rec polys N cnt)
    (if (<= N 0)
        cnt
        (let ([roots (for/list ([c (in-list polys)])
                       (apply poly_quintroots c))])
          (test-rec polys (- N count) (+ cnt (caadr roots))))))
  (time (test-rec polys N 0)))


(test-perf/quad 1000000)
(test-perf/cubic 100000)
(test-perf/quart 100000)
(test-perf/quint 100000)
