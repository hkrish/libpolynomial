#lang racket

(require ffi/unsafe
         ffi/unsafe/define)

(provide poly_quadroots
         poly_cubicroots
         poly_quartroots
         poly_quintroots)


(define-ffi-definer define-polys (ffi-lib "libpolynomial5"))

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

(define-polys poly_quintroots (_fun _double _double _double _double _double _double
                                    (_double = -inf.0)
                                    (_double = +inf.0)
                                    (x1 : (_ptr o _double))
                                    (x2 : (_ptr o _double))
                                    (x3 : (_ptr o _double))
                                    (x4 : (_ptr o _double))
                                    (x5 : (_ptr o _double))
                                    -> (r : _int32)
                                    -> (take (list x1 x2 x3 x4 x5) r)))
