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
  (let loop ([error (expt 2 30)] [index 0] [it 0] [delta '()] [domains domains])
    (cond
      [(empty? domains) (list index delta)]
      [else
       (define new-error (apply + (map (λ(x y) (abs (- x y))) (take (cdr range) 5) (take (cdar domains) 5))))
       (if (< new-error error)
           (loop new-error it (add1 it) (map - range (car domains)) (rest domains))
           (loop error index (add1 it) delta (rest domains)))])))

(define (search-ranges ranges domains)
  (for/list ([i ranges])
    (search-range i domains)))

;----------------------------------DECODER------------------------------------------------

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
    (define domain (vector-ref new-domains (first i)))
    (define dc-DCT (map + domain (second i)))
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

;----------------------------------ZIG-ZAG SCAN------------------------------------------------

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