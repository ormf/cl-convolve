;;;; package.lisp

(defpackage #:cl-convolve
  (:use #:cl #:incudine #:incudine.util #:incudine.analysis #:incudine.vug
        #:ana)

  (:export #:unzip-fft
           #:convolve
           #:cross-correlate
           #:sweep-log))


