#lang racket

(require math/flonum
         "utils.rkt")

(define LDvD
  '((delta  0.630985788572369227    delta  0.6309857885723692)
    (xm     -0.505864157480360260   xm     -0.5058641574803605)
    (ym     -0.794938379186223911   ym     -0.7949383791862251)
    (a_d    0.796286130760589280    a_d    0.7962861307605892)
    (a_     -11.593550710039739447  a_     -11.593550710039747)
    (b_     8.580290246019681617    b_     8.580290246019688)
    (xden   0.706588485422656570    xden   0.7065884854226564)
    (x      -0.397288449331790305   x      -0.3972884493317904)
    (delta  0.261205040211899636    delta  0.2612050402118996)
    (xm     -0.136083409119890669   xm     -0.13608340911989084)
    (ym     -0.043577142903070158   ym     -0.04357714290307024)
    (a_d    0.136456146064200212    a_d    0.13645614606420017)
    (a_     -2.364864658057078129   a_     -2.3648646580570793)
    (b_     0.796838800952579638    b_     0.7968388009525802)
    (xden   0.592276677599120173    xden   0.5922766775991195)
    (x      -0.067392424577620591   x      -0.06739242457762064)
    (delta  0.096257027834814779    delta  0.09625702783481468)
    (xm     0.028864603257194188    xm     0.028864603257194053)
    (ym     0.000257573412718027    ym     0.00025757341271801654)
    (a_d    0.018530830815184614    a_d    0.018530830815184578)
    (a_     -0.731233865433552850   a_     -0.7312338654335538)
    (b_     0.101069208873953170    b_     0.10106920887395335)
    (xden   1.036220527938389862    xden   1.0362205279383883)
    (x      0.026361450519896194    x      0.02636145051989617)
    (delta  0.046876937548758393    delta  0.04687693754875841)
    (xm     -0.020515487028862198   xm     -0.020515487028862234)
    (ym     -0.005789617696380345   ym     -0.005789617696380351)
    (a_d    0.004394894547900389    a_d    0.004394894547900392)
    (a_     -1.051148807869987716   a_     -1.0511488078699878)
    (b_     0.173790618792495010    b_     0.1737906187924951)
    (xden   0.440483004785612187    xden   0.4404830047856122)
    (x      0.025738099568580338    x      0.02573809956858032)
    (delta  0.000311675475657928    delta  0.00031167547565792535)
    (xm     0.026049775044238266    xm     0.026049775044238242)
    (ym     0.000020498034279865    ym     2.0498034279863107e-05)
    (a_d    0.000000194283204253    a_d    1.9428320425318804e-07)
    (a_     -0.726360264231553449   a_     -0.7263602642340348)
    (b_     0.086252312309122675    b_     0.0862523123091228)
    (xden   1.003994718843192009    xden   1.0039947188432052)
    (x      0.025812596777325048    x      0.025812596777325048)
    (delta  0.000037248604372355    delta  3.7248604372364796e-05)
    (xm     0.025775348172952693    xm     0.025775348172952683)
    (ym     -0.000003226607894576   ym     -3.2266078945769605e-06)
    (a_d    0.000000002774917055    a_d    2.7749170553779078e-09)
    (a_     -0.728051264363099561   a_     -0.728051264499537)
    (b_     0.086651245176566891    b_     0.08665124517657367)
    (xden   0.999374071866454352    xden   0.999374071866337)
    (x      0.025812596547358183    x      0.025812596547358183)
    (delta  0.000000000114983432    delta  #f)
    (xm     0.025812596662341616    xm     #f)
    (ym     0.000000000009957221    ym     #f)
    (a_d    0.000000000000000000    a_d    #f)
    (a_     -0.736762661026417049   a_     #f)
    (b_     0.086597013260844078    b_     #f)
    (xden   1.000000001956545270    xden   #f)
    (x      0.025812596547358199    x      #f)))



(define (zip al bl #:pad (pad #f))
  (cond
    [(and (pair? al) (pair? bl))
     (cons (cons (car al) (car bl)) (zip (cdr al) (cdr bl) #:pad pad))]
    [(pair? al)
     (cons (cons (car al) pad) (zip (cdr al) bl #:pad pad))]
    [(pair? bl)
     (cons (cons pad (car bl)) (zip al (cdr bl) #:pad pad))]
    [else empty]))


(define round 0)

(for ([dat (in-list LDvD)])
  (define sym (first dat))
  (define nld (second dat))
  (define nd (fourth dat))
  (when (eq? sym 'delta)
    (displayln " ")
    (set! round (add1 round)))
  (displayln-digit-diff nld nd
                        #:prefix-a (~a round ")  "(symbol->string sym) ": " #:width 20)
                        #:prefix-b (if (and nd nld)
                                       (~a "   ulp:  " (flulp-error nd nld) #:width 20)
                                       "")))
