#lang racket
(provide (all-defined-out))
(require "DCT.rkt")

(define SIZE 512)
(define N-size 32)
(define S-block 7)
(define TL 8)
(define n (sqr TL))

(define read-order
  '(0 8 1 2 9 16 24 17 10 3 4 11 18 25 32 40 33 26 19 12 5 6 13 20 27 34 41 48 56 49 42 35
      28 21 14 7 15 22 29 36 43 50 57 58 51 44 37 30 23 31 38 45 52 59 60 53 46 39 47 54 61 62 55 63))

(define (matrix-get matrix i j) (vector-ref (vector-ref matrix i) j))
(define (matrix-set matrix i j val) (vector-set! (vector-ref matrix i) j val))
(define (flatten-matrix matrix) (vector->list (apply vector-append (vector->list matrix))))
(define (matrix->list matrix) (apply append matrix))

(define (get-matrix buffer)
  (for/vector ([i SIZE])
    (for/vector ([j (in-range (add1 (* i 4 SIZE)) (add1 (* (add1 i) 4 SIZE)) 4)])
      (bytes-ref buffer j))))

(define (matrix->bytes matrix)
  (list->bytes (apply append (map (λ(x) (list 255 x x x)) (flatten-matrix matrix)))))

(define (get-blocks matrix [nr (/ SIZE TL)] [size TL] [step TL])
  (matrix->list
   (for/list ([i (in-range 0 (* nr step) step)])
     (for/list ([j (in-range 0 (* nr step) step)])
       (for/vector ([a (in-range i (+ i size))])
         (for/vector ([b (in-range j (+ j size))])
           (matrix-get matrix a b)))))))

(define (blocks->image-matrix blocks)
  (define new-matrix (for/vector ([i SIZE]) (make-vector SIZE)))
  (for ([block blocks] [index (in-naturals)])
    (define row (quotient index (/ SIZE TL)))
    (define column (remainder index (/ SIZE TL)))
    (for ([i TL])
      (for ([j TL])
        (matrix-set new-matrix (+ (* TL row) i) (+ (* TL column) j) (exact-floor (matrix-get block i j))))))
  new-matrix)

(define (mean m)
  (define x (flatten-matrix m))
  (exact->inexact (/ (apply + x) (length x))))

(define (normalize x [bias 0])
  (set! x (exact-round (+ bias x)))
  (cond [(< x 0) 0] [(> x 255) 255] [else x]))

(define (debiasing block bias)
  (define delta (- bias (mean block)))
  (for/vector ([i 8])
    (for/vector ([j 8])
      (normalize (+ delta (matrix-get block i j))))))

(define (crop-block block [size N-size])
  (for/list ([i (take read-order size)])
    (matrix-get block (quotient i 8) (remainder i 8))))
    
(define (padd-block coefs)
  (define new-matrix (for/vector ([i 8]) (make-vector 8)))
  (for ([i coefs] [j (take read-order (length coefs))])
    (matrix-set new-matrix (quotient j 8) (remainder j 8) i))
  new-matrix)

(define (PSNR original decoded)
  (set! original (flatten-matrix original))
  (set! decoded (flatten-matrix decoded))
  (* 10 (log
         (/ (* SIZE SIZE (sqr (apply max original)))
            (apply + (map (λ(x y) (sqr (- x y))) original decoded)))
         10)))

;----------------------------------ENCODER------------------------------------------------

(define (get-ranges matrix [nr (/ SIZE TL)] [size TL] [step TL])
  (matrix->list
   (for/list ([i (in-range 0 (* nr step) step)])
     (for/list ([j (in-range 0 (* nr step) step)])
       (crop-block
        (DCT
         (for/vector ([a (in-range i (+ i size))])
           (for/vector ([b (in-range j (+ j size))])
             (matrix-get matrix a b)))))))))

(define (get-domains matrix [nr (/ SIZE (* 2 TL))] [size (* 2 TL)] [step (* 2 TL)])
  (matrix->list
   (for/list ([i (in-range 0 (* nr step) step)])
     (for/list ([j (in-range 0 (* nr step) step)])
       (crop-block
        (DCT
         (for/vector ([a (in-range TL)])
           (for/vector ([b (in-range TL)])
             (quotient (exact-round
                        (+ (matrix-get matrix (+ i (* 2 a)) (+ j (* 2 b)))
                           (matrix-get matrix (+ i (* 2 a)) (+ j (+ 1 (* 2 b))))
                           (matrix-get matrix (+ i (+ 1 (* 2 a))) (+ j (* 2 b)))
                           (matrix-get matrix (+ i (+ 1 (* 2 a))) (+ j (+ 1 (* 2 b))))))
                       4)))))))))

(define (search-range range domains)
  (let loop ([1st-error (expt 2 30)] [1st-index 0] [it 0] [1st-delta '()] [domains domains])
    (cond
      [(empty? domains) (list 1st-index 1st-delta)]
      [else
       (define ac-range (cdr range))
       (define domain (cdar domains))
       (define new-1st-error (apply + (map (λ(x y) (abs (- x y))) (take ac-range S-block) (take domain S-block))))
       (cond
         [(< new-1st-error 1st-error) (loop new-1st-error it (add1 it) (map - ac-range domain) (rest domains))]
         [else (loop 1st-error 1st-index (add1 it) 1st-delta (rest domains))])])))
 
(define (search-ranges ranges domains) (for/list ([i ranges]) (search-range i domains)))

(define (decode founds new-domains DCs)
  (set! new-domains (list->vector new-domains))
  (for/list ([i founds] [DC DCs])
    (define domain (cdr (vector-ref new-domains (first i))))
    (define dc-DCT (cons DC (map + domain (second i))))
    (IDCT (padd-block dc-DCT))))