;;;; package.lisp

(defpackage #:cl-convolve
  (:use #:cl #:incudine #:incudine.util  #:ana #:cl-plot)

  (:export #:unzip-fft
           #:convolve
           #:cross-correlate))


