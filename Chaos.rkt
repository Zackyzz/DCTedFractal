#lang racket
(provide (all-defined-out))

(define SIZE 512)
(define N-size 2)
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
       (define sum 0)
       (define sum^2 0)
       (define block
         (for/vector ([a (in-range i (+ i size))])
           (for/vector ([b (in-range j (+ j size))])
             (define bi (matrix-get matrix a b))
             (set! sum (+ sum bi))
             (set! sum^2 (+ sum^2 (sqr bi)))
             bi)))
       (list block sum sum^2)))))

(define (get-domains matrix [nr (sub1 (/ SIZE TL))] [size (* 2 TL)] [step TL])
  (matrix->list
   (for/list ([i (in-range 0 (* nr step) step)])
     (for/list ([j (in-range 0 (* nr step) step)])
       (define sum 0)
       (define sum^2 0)
       (define block
         (for/vector ([a (in-range TL)])
           (for/vector ([b (in-range TL)])
             (define bi
               (quotient (exact-round
                          (+ (matrix-get matrix (+ i (* 2 a)) (+ j (* 2 b)))
                             (matrix-get matrix (+ i (* 2 a)) (+ j (+ 1 (* 2 b))))
                             (matrix-get matrix (+ i (+ 1 (* 2 a))) (+ j (* 2 b)))
                             (matrix-get matrix (+ i (+ 1 (* 2 a))) (+ j (+ 1 (* 2 b))))))
                         4))
             (set! sum (+ sum bi))
             (set! sum^2 (+ sum^2 (sqr bi)))
             bi)))
       (list block sum sum^2)))))

(define (search-range lrange domains)
  (define range (first lrange))
  (define sum-r (second lrange))
  (define sum-r^2 (third lrange))
  (let loop ([error (expt 2 30)] [index 0] [S 0] [O 0] [domains domains] [it 0])
    [cond
      [(empty? domains) (list index (exact->inexact S) (exact->inexact O))]
      [else
       (define domain (caar domains))
       (define sum-d (cadar domains))
       (define sum-d^2 (caddar domains))
       (define sum-rd (apply + (map * (flatten-matrix range) (flatten-matrix domain))))
       (define denom-s (- (* n sum-d^2) (sqr sum-d)))
       (define s (if (= 0 denom-s) 0
                     (/ (- (* n sum-rd) (* sum-r sum-d)) denom-s)))
       (define o (/ (- sum-r (* s sum-d)) n))
       (define Error (+ sum-r^2
                        (* s (+ (* s sum-d^2) (- (* 2 sum-rd)) (* 2 o sum-d)))
                        (* o (- (* o n) (* 2 sum-r)))))
       (if (< Error error)
           (loop Error it s o (rest domains) (add1 it))
           (loop error index S O (rest domains) (add1 it)))]]))

(define (search-ranges ranges domains [p-gauge #f])
  (for/list ([i ranges])
    (when p-gauge (send p-gauge set-value (add1 (send p-gauge get-value))))
    (search-range i domains)))

(define (get-decoding-domains matrix [nr (sub1 (/ SIZE TL))] [size (* 2 TL)] [step TL])
  (matrix->list
   (for/list ([i (in-range 0 (* nr step) step)])
     (for/list ([j (in-range 0 (* nr step) step)])
       (for/vector ([a (in-range TL)])
         (for/vector ([b (in-range TL)])
           (quotient (exact-round
                      (+ (matrix-get matrix (+ i (* 2 a)) (+ j (* 2 b)))
                         (matrix-get matrix (+ i (* 2 a)) (+ j (+ 1 (* 2 b))))
                         (matrix-get matrix (+ i (+ 1 (* 2 a))) (+ j (* 2 b)))
                         (matrix-get matrix (+ i (+ 1 (* 2 a))) (+ j (+ 1 (* 2 b))))))
                     4)))))))

(define (decode founds new-domains)
  (set! new-domains (list->vector new-domains))
  (for/vector ([i founds])
    (define domain (vector-ref new-domains (first i)))
    (define s (second i))
    (define o (third i))
    (for/vector ([i TL])
      (for/vector ([j TL])
        (+ o (* s (matrix-get domain i j)))))))

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