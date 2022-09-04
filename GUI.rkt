#lang racket/gui
(require "Chaos.rkt")

(define frame
  (new frame%
       [label "Fractal Coding"]
       [x 250] [y 150]
       [width 1150] [height 650]))

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
            (time (set! founds (search-ranges ranges domains)))
            (for ([i 50])
              (printf "~a ~a ~a\n"
                      i
                      (+ (length (filter (λ(x) (= i x)) (apply append (map second founds))))
                      (length (filter (λ(x) (= (- i) x)) (apply append (map second founds)))))
                      (+ (length (filter (λ(x) (= i x)) (apply append ranges)))
                      (length (filter (λ(x) (= (- i) x)) (apply append ranges))))))))]))
                   
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
  
(define final-matrix (for/vector ([i SIZE]) (make-vector SIZE)))
(define decode-button
  (new button%
       [parent decode-panel]
       [label "Decode"]
       [callback
        (λ (button event)
          (time
           (when founds
             (for ([i 25])
               (define blocks (decode founds (get-domains final-matrix)))
               (set! final-matrix (blocks->image-matrix blocks))
               (set! final-matrix
                     (for/vector ([i 512])
                       (for/vector ([j 512])
                         (normalize (matrix-get final-matrix i j)))))
               (send psnr-field set-value (number->string (PSNR original-matrix final-matrix)))
               (send decode-bitmap set-argb-pixels 0 0 SIZE SIZE (matrix->bytes final-matrix))
               (send decode-canvas on-paint)))))]))

(define psnr-field
  (new text-field%
       [parent decode-panel]
       [label "PSNR:"]
       [horiz-margin 150]
       [init-value ""]))