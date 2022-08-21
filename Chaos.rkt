#lang racket
(provide (all-defined-out))
(require "DCT.rkt")

(define SIZE 512)
(define N-size 4)
(define TL 8)
(define n (sqr TL))

(define (matrix-get matrix i j)
  (vector-ref (vector-ref matrix i) j))

(define (matrix-set matrix i j val)
  (vector-set! (vector-ref matrix i) j val))

(define (flatten-matrix matrix)
  (vector->list (apply vector-append (vector->list matrix))))

(define (matrix->list matrix)
  (apply append matrix))

(define (get-blocks matrix [nr (/ SIZE TL)] [size TL] [step TL])
  (matrix->list
   (for/list ([i (in-range 0 (* nr step) step)])
     (for/list ([j (in-range 0 (* nr step) step)])
       (for/vector ([a (in-range i (+ i size))])
         (for/vector ([b (in-range j (+ j size))])
           (matrix-get matrix a b)))))))

(define (get-ranges matrix [nr (/ SIZE TL)] [size TL] [step TL])
  (matrix->list
   (for/list ([i (in-range 0 (* nr step) step)])
     (for/list ([j (in-range 0 (* nr step) step)])
       (list
        (crop-block
         (DCT
          (for/vector ([a (in-range i (+ i size))])
            (for/vector ([b (in-range j (+ j size))])
              (matrix-get matrix a b)))))
        (random 0 1024))))))

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
  (define ldomain (flatten-matrix (list-ref domains (second range))))
  (define lrange (flatten-matrix (first range)))
  (define delta (map - lrange ldomain))
  (list (second range) delta))

(define (search-ranges ranges domains [p-gauge #f])
  (for/list ([i ranges])
    (when p-gauge (send p-gauge set-value (add1 (send p-gauge get-value))))
    (search-range i domains)))

(define (get-decoding-domains matrix [nr (/ SIZE (* 2 TL))] [size (* 2 TL)] [step (* 2 TL)])
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

(define (decode founds new-domains)
  (set! new-domains (list->vector new-domains))
  (for/list ([i founds])
    (define domain (flatten-matrix (vector-ref new-domains (first i))))
    (define dc-DCT (map + domain (second i)))
    (set! dc-DCT
          (for/vector ([i N-size])
            (for/vector ([j N-size])
              (list-ref dc-DCT (+ j (* i N-size))))))
    (IDCT (padd-block dc-DCT))))

(define (blocks->image-matrix blocks)
  (define new-matrix (for/vector ([i SIZE]) (make-vector SIZE)))
  (for ([block blocks] [index (in-naturals)])
    (define row (quotient index (/ SIZE TL)))
    (define column (remainder index (/ SIZE TL)))
    (for ([i TL])
      (for ([j TL])
        (matrix-set new-matrix (+ (* TL row) i) (+ (* TL column) j) (exact-floor (matrix-get block i j))))))
  new-matrix)

(define (crop-block block [size N-size])
  (for/vector ([i size])
    (for/vector ([j size])
      (matrix-get block i j))))

(define (padd-block block [size N-size])
  (define new-matrix (for/vector ([i 8]) (make-vector 8)))
  (for ([i size])
    (for ([j size])
      (matrix-set new-matrix i j (matrix-get block i j))))
  new-matrix)

(define (list->matrix lst)
  (for/vector ([i 8])
    (for/vector ([j 8])
      (list-ref lst (+ j (* i 8))))))