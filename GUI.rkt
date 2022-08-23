#lang racket/gui
(require "Chaos.rkt")

(define (get-matrix buffer)
  (for/vector ([i SIZE])
    (for/vector ([j (in-range (add1 (* i 4 SIZE)) (add1 (* (add1 i) 4 SIZE)) 4)])
      (bytes-ref buffer j))))

;------------------------------------GUI----------------------------------------

(define frame
  (new frame%
       [label "Fractal Coding"]
       [x 250] [y 150]
       [width 1150] [height 670]))

(send frame show #t)

(define main-panel
  (new horizontal-panel%
       [parent frame]))

;------------------------------------ENCODE PANEL----------------------------------------

(define encode-panel
  (new vertical-panel%
       [parent main-panel]))

(define encode-bitmap (make-bitmap SIZE SIZE))
(define encode-dc (send encode-bitmap make-dc))
(send encode-dc set-background (make-color 0 0 0))
(send encode-dc clear)

(define encode-canvas
  (new canvas%
       [parent encode-panel]
       [min-width SIZE]
       [paint-callback
        (λ (canvas dc)
          (send dc draw-bitmap encode-bitmap 20 20))]))

(define ranges #f)
(define domains #f)
(define encode-buffer (make-bytes (* SIZE SIZE 4)))
(define original-matrix #f)

(define load-button-encode
  (new button%
       [parent encode-panel]
       [label "Load file"]
       [callback
        (λ (button event)
          (define path (get-file #f #f "../FractalDCT/utils" #f #f null))
          (time
           (when path
             (set! encode-bitmap (read-bitmap path))
             (send encode-canvas on-paint)
             (send encode-bitmap get-argb-pixels 0 0 SIZE SIZE encode-buffer)
             (set! original-matrix (get-matrix encode-buffer))
             (set! ranges (get-ranges original-matrix))
             (set! domains (get-domains original-matrix)))))]))

(define founds #f)
(define process-button
  (new button%
       [parent encode-panel]
       [label "Process"]
       [callback
        (λ (button event)
          (when (and ranges domains)
            (time (set! founds (search-ranges ranges domains)))))]))
                   
;------------------------------------DECODE PANEL----------------------------------------

(define decode-panel
  (new vertical-panel%
       [parent main-panel]))

(define decode-bitmap (make-bitmap SIZE SIZE))
(define decode-dc (send decode-bitmap make-dc))
(send decode-dc set-background (make-color 0 0 0))
(send decode-dc clear)

(define decode-canvas
  (new canvas%
       [parent decode-panel]
       [min-width SIZE]
       [paint-callback
        (λ (canvas dc)
          (send dc draw-bitmap decode-bitmap 20 20))]))

(define (normalize x [bias 0])
  (set! x (exact-round (+ bias x)))
  (cond [(< x 0) 0] [(> x 255) 255] [else x]))

(define (matrix->bytes matrix)
  (list->bytes (apply append (map (λ(x) (list 255 x x x)) (flatten-matrix matrix)))))

(define (mean m)
  (define x (flatten-matrix m))
  (exact->inexact (/ (apply + x) (length x))))

(define final-matrix (for/vector ([i SIZE]) (make-vector SIZE)))
(define decode-button
  (new button%
       [parent decode-panel]
       [label "Decode"]
       [callback
        (λ (button event)
          (time
           (when founds
             (for ([i 15])
               (define blocks (decode founds (get-decoding-domains final-matrix)))
               (set! final-matrix (blocks->image-matrix blocks))
               (define bias (- (mean original-matrix) (mean final-matrix)))
               (define refinal-matrix
                 (for/vector ([i SIZE])
                   (for/vector ([j SIZE])
                     (normalize (matrix-get final-matrix i j) bias))))
               (send psnr-field set-value (number->string (PSNR original-matrix refinal-matrix)))
               (send mae-field set-value (number->string (MAE original-matrix refinal-matrix)))
               (send decode-bitmap set-argb-pixels 0 0 SIZE SIZE (matrix->bytes refinal-matrix))
               (send decode-canvas on-paint)))))]))

(define (PSNR original decoded)
  (set! original (flatten-matrix original))
  (set! decoded (flatten-matrix decoded))
  (* 10 (log
         (/ (* SIZE SIZE (sqr (apply max original)))
            (apply + (map (λ(x y) (sqr (- x y))) original decoded)))
         10)))

(define (MAE original decoded)
  (set! original (flatten-matrix original))
  (set! decoded (flatten-matrix decoded))
  (/ (apply + (map (λ(x y) (abs (- x y))) original decoded)) (sqr 512.0)))

(define psnr-field
  (new text-field%
       [parent decode-panel]
       [label "PSNR:"]
       [horiz-margin 150]
       [init-value ""]))

(define mae-field
  (new text-field%
       [parent decode-panel]
       [label "MAE: "]
       [horiz-margin 150]
       [init-value ""]))