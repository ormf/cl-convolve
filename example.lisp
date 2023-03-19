;;; 
;;; example.lisp
;;;
;;; **********************************************************************
;;; Copyright (c) 2020 Orm Finnendahl <orm.finnendahl@selma.hfmdk-frankfurt.de>
;;;
;;; Revision history: See git repository.
;;;
;;; This program is free software; you can redistribute it and/or
;;; modify it under the terms of the Gnu Public License, version 2 or
;;; later. See https://www.gnu.org/licenses/gpl-2.0.html for the text
;;; of this agreement.
;;; 
;;; This program is distributed in the hope that it will be useful,
;;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
;;; GNU General Public License for more details.
;;;
;;; **********************************************************************

(ql:quickload "cl-convolve")

(in-package :cl-convolve)

;;; (use-package :incudine.scratch)
;;;(use-package :incudine.util)


(defparameter *snd1*
  (incudine:make-buffer (expt 2 6)
                        :fill-function (gen:partials '(1))))


(defparameter *snd2*
  (incudine:make-buffer (expt 2 6)
                        :fill-function (gen:partials '(0 0 1))))

(plot *snd1*)
(plot *snd2*)

(plot (convolve
       (incudine:make-buffer (expt 2 6)
                             :fill-function (gen:partials '(1)))
       (incudine:make-buffer (expt 2 6)
                             :fill-function (gen:partials '(0 0 1)))))

(defparameter *a* (make-buffer 16))
(defparameter *b* (make-buffer 16))

(clear-buf *a*)
(clear-buf *b*)
(setf (buffer-value *a* 15) 1.0d0)
(setf (buffer-value *a* 13) 0.6d0)
(setf (buffer-value *b* 11) 1.0d0)

(plot *a*)
(plot *b*)

(plot (convolve *a* *b*))

(plot (cross-correlate *a*))

(defparameter *sweep1l*
  (buffer-load (pathname "~/work/programmieren/lisp/cl-convolve/snd/Position1_01.wav")))

(defparameter *sweep1-bw*
  (buffer-load (pathname "~/work/programmieren/lisp/cl-convolve/snd/Position1_01-bw.wav")))

(defparameter *impulse-response* (convolve *sweep1* *sweep1-bw*))

(defparameter *sweep1l*
  (buffer-load (pathname "~/work/programmieren/lisp/cl-convolve/snd/Position1_00.wav")))

(defparameter *sweep1r*
  (buffer-load (pathname "~/work/programmieren/lisp/cl-convolve/snd/Position1_01.wav")))

(plot *impulse-response*)

(defparameter *impulse-response1l* (cross-correlate *sweep1l*))
(defparameter *impulse-response1r* (cross-correlate *sweep1r*))
(buffer-save *impulse-response1l* "/tmp/ir-pos1l.wav")
(buffer-save *impulse-response1r* "/tmp/ir-pos1r.wav")

(dotimes (i 6)
  (let ((bufl (buffer-load
               (pathname (format nil "/home/orm/work/programmieren/lisp/cl-convolve/snd/Position~d_00.wav" (+ i 1)))))
        (bufr (buffer-load
               (pathname (format nil "/home/orm/work/programmieren/lisp/cl-convolve/snd/Position~d_01.wav" (+ i 1))))))
    (format t "~&processing Position~d...~%" (+ 1 i))
    (buffer-save (cross-correlate bufl)
                 (pathname (format nil "/home/orm/work/programmieren/lisp/cl-convolve/snd/ir-new-pos~dl.wav" (+ i 1))))
    (buffer-save (cross-correlate bufr)
                 (pathname (format nil "/home/orm/work/programmieren/lisp/cl-convolve/snd/ir-new-pos~dr.wav" (+ i 1))))))


(plot *impulse-response*)
(plot *impulse-response2*)
(buffer-save *impulse-response* "/tmp/impulse-response.wav")
(buffer-save *impulse-response2* "/tmp/impulse-response2.wav")

(sweep (convolve (incudine:make-buffer (expt 2 6)
                                      :fill-function (gen:partials '(1)))
                (incudine:make-buffer (expt 2 6)
                                      :fill-function (gen:partials '(0 0 1)))))

(defparameter result-size nil)
(defparameter fft-size nil)
(defparameter fft1 nil)
(defparameter fft2 nil)
(defparameter ifft nil)
(defparameter buf1 nil)
(defparameter buf2 nil)

(setf i-tmp (+ (* x v) (* y x)))
(setf r-tmp (- (* x x) (* y v)))

u = x
v = -y

(setf i-tmp (+ (* x -y) (* y x)))
(setf r-tmp (- (* x x) (* y -y)))

x^2+y^2

(progn
  (setf buf1 (make-buffer 4))
  (setf (buffer-value buf1 1) 1.0d0)
  (setf buf2 (make-buffer 4))
  (setf (buffer-value buf2 1) 1.0d0)
  (setf result-size (+ (buffer-frames buf1)
                       (buffer-frames buf2)
                       ))
  (setf fft-size (* 2 (next-power-of-two result-size)))
  (setf fft1 (make-fft fft-size :window-function #'rectangular-window))
  (setf fft2 (make-fft fft-size :window-function #'rectangular-window))
  (setf ifft (make-ifft fft-size :window-function #'rectangular-window)))


(progn
  (setf buf1 *a*)
  (setf buf2 *b*)
  (setf result-size (+ (buffer-frames buf1)
                       (buffer-frames buf2)))
  (setf fft-size (* 2 (next-power-of-two result-size)))
  (setf fft1 (make-fft fft-size :window-function #'rectangular-window))
  (setf fft2 (make-fft fft-size :window-function #'rectangular-window))
  (setf ifft (make-ifft fft-size :window-function #'rectangular-window)))

(compute-fft (buffer->fft-input buf1 fft1) t)
(compute-fft (buffer->fft-input buf2 fft2) t)

(plot buf1)
(plot buf2)
(plot (unzip-fft (compute-fft (buffer->fft-input buf1 fft1) t) :input))


(plot (unzip-fft (compute-fft (buffer->fft-input buf1 fft1) t) :output))
(plot (unzip-fft (compute-fft (buffer->fft-input buf2 fft2) t) :output))

(fft-mult fft1 fft2)

(plot (unzip-fft fft1 :input))
(plot (unzip-fft fft1 :output))
(plot (unzip-fft fft2 :output))

(clear-fft fft2 :output)

(loop for i below (analysis-output-buffer-size fft2)
      collect (smp-ref (analysis-output-buffer fft2) i))

(1.0d0 0.0d0 -1.0d0 -0.0d0 1.0d0 -0.0d0 -1.0d0 0.0d0 1.0d0 0.0d0)

(1.0d0 0.0d0 0.0d0 -1.0d0 -1.0d0 0.0d0 0.0d0 1.0d0 1.0d0 0.0d0)

(loop for i below (analysis-output-buffer-size fft1)
      collect (smp-ref (analysis-output-buffer fft1) i))

(1.0d0 0.0d0 0.0d0 -1.0d0 -1.0d0 0.0d0 0.0d0 1.0d0 1.0d0 0.0d0)

(defparameter *fft* (make-fft 8 :window-function #'rectangular-window))

;; The method WINDOW-FUNCTION is SETFable, therefore the alternative is
;; (setf (window-function *fft*) #'rectangular-window)

(progn
  (setf *fft* (make-fft 8 :window-function #'rectangular-window))
  (dotimes (i (fft-size *fft*))
    (setf (fft-input *fft*) (if (= i 0) 1d0 0d0)))
  (plot (unzip-fft (compute-fft *fft* t) :input))
  ;; 5 complex bins.
  (loop for i below (analysis-output-buffer-size *fft*)
        collect (smp-ref (analysis-output-buffer *fft*) i)))

(plot (unzip-fft *fft* :input))

(plot (unzip-fft *fft* :output))

(plot (convolve  *a* *b*))

(defun convolve (buf1 buf2)
  (let* ((result-size
           (+ (buffer-frames buf1)
              (buffer-frames buf2)))
         (fft-size (* 2 (next-power-of-two result-size)))
         (result (incudine:make-buffer result-size))
         (fft1 (make-fft fft-size :window-function #'rectangular-window))
         (fft2 (make-fft fft-size :window-function #'rectangular-window))
         (ifft (make-ifft fft-size :window-function #'rectangular-window)))
    (fft-output->buffer
     (compute-ifft
      (fft-output->ifft-input
       (fft-mult
        (compute-fft (buffer->fft-input buf1 fft1) t)
        (compute-fft (buffer->fft-input buf2 fft2) t))
       ifft)
      nil t)
     result)))


(defparameter *a* (make-buffer 16))
(setf (buffer-value *a* 0) 1.0d0)

(defparameter *fft* (make-fft 32))
(buffer->fft-input *a* *fft*)
(unzip-fft *fft* :input)
(calc-fft (buffer->fft-input *a* *fft*))

(plot (unzip-fft *fft* :input))
(plot (unzip-fft *fft* :output))

(defparameter *a* (make-buffer 64))
(defparameter *fft* (make-fft 128))
(defparameter *ifft* (make-fft 128))
(setf (buffer-value *a* 4) 1.0d0)
(calc-fft (buffer->fft-input *a* *fft*))

(plot (unzip-fft *fft* :input))
(plot (unzip-fft *fft* :output))

(clear-fft *fft* :input)
(calc-ifft *fft*)

(plot (unzip-fft *fft* :input))
(plot (unzip-fft *fft* :output))

(rt-start)
(defparameter *test* (compute-fft *fft* t))

(compute-fft *fft* t)

(plot (unzip-fft *test* :output))

(plot (convolve *a* *b*))

(ana::make-fft)

(clear-fft *fft1* :input)
(clear-fft *fft2* :input)

(plot-fft *my-fft* :input)
(buffer->fft-input *test-snd* *my-fft*)
(calc-fft *my-fft*)
(plot-fft *my-fft* :output)

(calc-ifft *my-fft*)
(plot-fft *my-fft* :input)

;;; (compute-fft *my-fft* t)
(plot-fft *my-fft* :output)
(cl-plot:plot *test-snd*)

(cl-plot:plot
 (loop for x below (ana::fft-size *my-fft*)
       collect (incudine:smp-ref (analysis-input-buffer *my-fft*) x)))


(plot *snd1*)
(plot *snd2*)

(clear-fft *fft1* :input)
(clear-fft *fft2* :input)
(buffer->fft-input *snd1* *fft1*)
(buffer->fft-input *snd2* *fft2*)
(calc-fft *fft1*)
(calc-fft *fft2*)

(plot-fft *fft1* :input)
(plot-fft *fft1* :output)
(plot-fft *fft2* :input)
(plot-fft *fft2* :output)

(fft-mult *fft1* *fft2* *fft-res*)
(plot-fft *fft-res* :output)
(calc-ifft *fft-res*)
(plot-fft *fft-res* :input)

(fft-input->buffer *fft-res* *snd-res*)

(plot *snd-res*)




(plot (convolve *snd1* *snd2*))

(setf *buf* (make-buffer 4))
(setf *fft* (make-fft 8))
(ql:quickload "incudine")
(in-package :scratch)

(defvar *buf* (make-buffer 4))
(setf (buffer-value *buf* 0) 1.0d0)

(in-package :scratch)

(defvar *fft* (make-fft 8))

(defun setup-fft (fft)
  "setup fft by clearing its input and output buffers and setting the
real part of the first bin to the unit impulse. Return the contents of
the fft-input buffer as list."
  (dotimes (i (ana:fft-size fft))
    (setf (smp-ref (analysis-output-buffer fft) i) 0.0d0)
    (setf (smp-ref (analysis-input-buffer fft) i) 0.0d0))
  (setf (smp-ref (analysis-input-buffer fft) 0) 1.0d0)
  (loop for i below (ana:fft-size *fft*)
      collect (smp-ref (analysis-input-buffer *fft*) i)))

(setup-fft *fft*)

;;; -> (1.0d0 0.0d0 0.0d0 0.0d0 0.0d0 0.0d0 0.0d0 0.0d0)

;;; this works as expected:

(progn
  (setup-fft *fft*)
  (ana::fft-execute
   (ana::fft-plan *fft*) (ana::fft-input-buffer *fft*) (ana::fft-output-buffer *fft*))
  (loop for i below (ana:fft-size *fft*)
        collect (smp-ref (analysis-output-buffer *fft*) i)))

;;; -> (1.0d0 0.0d0 1.0d0 0.0d0 1.0d0 0.0d0 1.0d0 0.0d0)

;;;; this doesn't seem to work:

(progn
  (setup-fft *fft*)
  (ana:compute-fft *fft* t)
  (loop for i below (ana:fft-size *fft*)
        collect (smp-ref (analysis-output-buffer *fft*) i)))

;;; -> (1.6927740957517413d185 0.0d0 -1.6927740957517413d185 -9.156735986014867d165
;;;     1.6927740957517413d185 0.0d0 -1.6927740957517413d185 9.156735986014867d165

(if *fft* (free *fft*))
(defparameter *fft* (make-fft 8 :window-function #'rectangular-window))

;; The method WINDOW-FUNCTION is SETFable, therefore the alternative is
;; (setf (window-function *fft*) #'rectangular-window)

(progn
  (dotimes (i (fft-size *fft*))
    (setf (smp-ref (analysis-input-buffer *fft*) i)
          1d0))
   (loop for i below (analysis-input-buffer-size *fft*)
         do (format t "~a " (smp-ref (analysis-input-buffer *fft*) i)))
   (compute-fft *fft* t)
  ;; 5 complex bins.
  (list


   (loop for i below (analysis-output-buffer-size *fft*)
         collect (smp-ref (analysis-output-buffer *fft*) i))))



;; => (1.0d0 0.0d0 1.0d0 0.0d0 1.0d0 0.0d0 1.0d0 0.0d0 1.0d0 0.0d0)
  



(defparameter *ifft* (make-ifft 8 :window-function #'rectangular-window))

(progn
  (dotimes (i (analysis-input-buffer-size *ifft*))
    (setf (smp-ref (analysis-input-buffer *ifft*) i)
          (if (evenp i) 1d0 0d0)))
  (compute-ifft *ifft* nil t)
  ;; The scale factor in this case is 1 because COMPUTE-IFFT is called
  ;; without ABUFFER (the IFFT object doesn't know the FFT object).
  (loop for i below 8 collect (/ (ifft-output *ifft*) (ifft-size *ifft*))))

;; => (1.0d0 0.0d0 0.0d0 0.0d0 0.0d0 0.0d0 0.0d0 0.0d0)




(progn
  (cl-convolve::fft-output->ifft-input *fft* *ifft*)
  (compute-ifft *ifft* nil t)
  (loop for i below 8 collect (/ (ifft-output *ifft*) (ifft-size *ifft*))))

(analysis-output-buffer-size *fft*)


(1- (incudine.analysis::analysis-nbins *fft*))

(loop for x below 10 by 2 collect x)


(defparameter *fft* (make-fft 8 :window-function #'incudine.analysis::blackman-harris-3-1))

(incudine.analysis:compute-fft *fft* t)



(defun fft->r-i-arrays (fft r-array i-array &optional (dir :input))
  "unzip complex values into left and right half."
  (let* ((vector-size (case dir
                 (:output (analysis-output-buffer-size fft))
                 (otherwise (analysis-input-buffer-size fft))))
         (buf (case dir
                (:output (analysis-output-buffer fft))
                (otherwise (analysis-input-buffer fft)))))
    (dotimes (n (ash vector-size -1))
      (setf (aref r-array n)
            (incudine:smp-ref buf (ash n 1)))
      (setf (aref i-array n)
            (incudine:smp-ref buf (1+ (ash n 1)))))
    (values r-array i-array)))

(fft->r-i-arrays *fft* (make-array 4) (make-array 4) :input)



(fft-num-bins *fft*)

(analysis-output-buffer-size *fft*)

(analysis-input-buffer-size *ifft*)

(cl-plot:plot
 (cl-convolve:unzip-fft *fft* :output))

(cl-plot:plot
 (cl-convolve:unzip-fft *ifft* :input))


(make-fft 4096 :flags +fft-plan-fast+ :window-function (incudine.gen::blackman-harris-3-1))

(incudine.gen:sine-window)
