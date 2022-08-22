#lang racket
(provide (all-defined-out))
(require "DCT.rkt")

(define SIZE 512)
(define N-size 16)
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

;----------------------------------ENCODER------------------------------------------------

(define (get-ranges matrix [nr (/ SIZE TL)] [size TL] [step TL])
  (map list
       (matrix->list
        (for/list ([i (in-range 0 (* nr step) step)])
          (for/list ([j (in-range 0 (* nr step) step)])
            (crop-block
             (DCT
              (for/vector ([a (in-range i (+ i size))])
                (for/vector ([b (in-range j (+ j size))])
                  (matrix-get matrix a b))))))))
       (build-list 4096 (λ _ (random 3969)))))

(define (get-domains matrix [nr (sub1 (/ SIZE TL))] [size (* 2 TL)] [step TL])
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
  (define ldomain (vector-ref (list->vector domains) (second range)))
  (define lrange (first range))
  (define delta (map - lrange ldomain))
  (list delta (second range)))

(define (search-ranges ranges domains)
  (for/list ([i ranges])
    (search-range i domains)))

;----------------------------------DECODER------------------------------------------------

(define (get-decoding-domains matrix [nr (sub1 (/ SIZE TL))] [size (* 2 TL)] [step TL])
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
    (define domain (vector-ref new-domains (second i)))
    (define dc-DCT (map + domain (first i)))
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

;----------------------------------ZIG-ZAG------------------------------------------------

(define/match (compare i j)
  [((list x y) (list a b)) (or (< x a) (and (= x a) (< y b)))])
 
(define/match (key i)
  [((list x y)) (list (+ x y) (if (even? (+ x y)) (- y) y))])
 
(define (zigzag-ht n)
  (define indexorder
    (sort (for*/list ([x n] [y n]) (list x y))
          compare #:key key))
  (for/hash ([(n i) (in-indexed indexorder)]) (values n i)))
 
(define (zigzag n)
  (define ht (zigzag-ht n))
  (for/list ([x n]) 
    (for/list ([y n])
      (hash-ref ht (list x y)))))

(define read-order
  (for/list ([i 64])
    (index-of (apply append (zigzag 8)) i)))

(define (crop-block block [size N-size])
  (for/list ([i (take read-order size)])
    (matrix-get block (quotient i 8) (remainder i 8))))
    
(define (padd-block coefs)
  (define new-matrix (for/vector ([i 8]) (make-vector 8)))
  (for ([i coefs] [j (take read-order (length coefs))])
    (matrix-set new-matrix (quotient j 8) (remainder j 8) i))
  new-matrix)