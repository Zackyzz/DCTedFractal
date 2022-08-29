#lang racket
(provide DCT IDCT)

(define (shift- vec)
  (vector-map (λ(x) (- x 128)) vec))

(define (shift+ vec)
  (vector-map (λ(x) (+ (exact-round x) 128)) vec))

(define (matrix-get matrix i j)
  (vector-ref (vector-ref matrix i) j))

(define T
  (for/vector ([i 8])
    (for/vector ([j 8])
      (if (= i 0)
          (/ 1.0 (sqrt 8))
          (* 0.5 (cos (/ (* (+ 1 (* 2 j)) i pi) 16)))))))

(define (lm->mm m) (list->vector (map list->vector m)))
(define (mm->lm m) (vector->list (vector-map vector->list m)))

(define (mult m1 m2)
  (set! m1 (mm->lm m1))
  (set! m2 (mm->lm m2))
  (lm->mm
   (for/list ([r m1])
     (for/list ([c (apply map list m2)])
       (apply + (map * r c))))))

(define (transpose m) (lm->mm (apply map list (mm->lm m))))

(define (quality vec)
  (vector-map (λ(x) (exact-round (* (- 2 (* 2 0.5)) x))) vec))

(define q-table
  (vector-map quality
              (vector (vector 16 11 10 16 24 40 51 61)
                      (vector 12 12 14 19 26 58 60 55)
                      (vector 14 13 16 24 40 57 69 56)
                      (vector 14 17 22 29 51 87 80 62)
                      (vector 18 22 37 56 68 109 103 77)
                      (vector 24 35 55 64 81 104 113 92)
                      (vector 49 64 78 87 103 121 120 101)
                      (vector 72 92 95 98 112 100 103 99))))

(define (quantize block)
  (for/vector ([i 8])
    (for/vector ([j 8])
      (exact-round (/ (matrix-get block i j) (matrix-get q-table i j))))))

(define (dequantize block)
  (for/vector ([i 8])
    (for/vector ([j 8])
      (* (matrix-get block i j) (matrix-get q-table i j)))))

(define (DCT block) (quantize (mult (mult T (vector-map shift- block)) (transpose T))))
(define (IDCT block) (vector-map shift+ (mult (mult (transpose T) (dequantize block)) T)))