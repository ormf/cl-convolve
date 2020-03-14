;;;
;;; cl-convolve.lisp
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

(in-package #:cl-convolve)

(defun copy-ir (left-ir right-ir result)
  "copy the right half of left-ir and right-ir into 2-channel result
buffer of half the left-ir size."
  (let ((ir-size (ash (buffer-size left-ir) -1)))
    (dotimes (idx ir-size)
      (setf (buffer-value result (ash idx 1))
            (buffer-value left-ir (+ idx ir-size)))
      (setf (buffer-value result (1+ (ash idx 1)))
            (buffer-value right-ir (+ idx ir-size))))))

(defun clear-fft (fft dir &optional (offset 0))
  "set all sampes of fft :input or :output to zero."
  (let ((size (case dir
               (:input (analysis-input-buffer-size fft))
               (:output (analysis-output-buffer-size fft))))
        (buf (case dir
               (:input (analysis-input-buffer fft))
               (:output (analysis-output-buffer fft)))))
    (if buf
        (loop
          for x from offset below size
          do (setf (incudine:smp-ref buf x) incudine.util:+sample-zero+)))))

(defun clear-buf (buffer)
  "set all samples of buffer to zero."
  (dotimes (i (buffer-size buffer))
    (setf (buffer-value buffer i) incudine.util:+sample-zero+))
  buffer)

(defun buffer->fft-input (buffer fft)
  "copy the contents of buffer to the real part of the fft input buffer. 
Set all other values to 0d0. Return the fft."
  (dotimes (i (analysis-input-buffer-size fft))
    (setf (fft-input fft)
          (if (and (zerop (logand i 1))
                   (< (ash i -1) (buffer-size buffer)))
              (buffer-value buffer (ash i -1))
              0.0d0)))
  fft)

(defun fft-input->buffer (fft buffer)
  "copy the real part of the fft input buffer to the buffer. Return the buffer."
  (dotimes (i (min (buffer-size buffer) (analysis-input-buffer-size fft)))
    (setf (buffer-value buffer i) (smp-ref (analysis-input-buffer fft) (* 2 i))))
  buffer)

(defun fft-output->buffer (fft buffer)
  "copy the real part of the fft input buffer to the buffer. Return the buffer."
  (dotimes (i (min (buffer-size buffer) (analysis-output-buffer-size fft)))
    (setf (buffer-value buffer i) (smp-ref (analysis-output-buffer fft) (* 2 i))))
  buffer)

;;; (incudine:smp-ref (analysis-output-buffer fft) i)

(defun fft-output->ifft-input (fft ifft)
  (let ((size (analysis-input-buffer-size ifft)))
    (if (= (analysis-output-buffer-size fft)
           size)
        (incudine::foreign-copy-samples
         (analysis-input-buffer ifft)
         (analysis-output-buffer fft)
         size)
        (error "buffer sizes don't match ~a ~a" fft ifft))
    ifft))

(defun fft-output->ifft-input-norm (fft ifft)
  (let ((size (analysis-input-buffer-size ifft))
        (factor (float (/ (ana:fft-size fft)) 1.0)))
    (if (= (analysis-output-buffer-size fft)
           size)
        (dotimes (i size)
          (setf (smp-ref (analysis-input-buffer ifft) i)
                (* factor (smp-ref (analysis-output-buffer fft) i))))
        (error "buffer sizes don't match ~a ~a" fft ifft))
    ifft))

(defun fft-mult (fft1 fft2)
  "in place complex multiplication of tha analysis output buffers of
fft1 and fft2. Store result into the analysis output buffer of fft2
and return it."
  (let ((num-fft-bins (analysis-output-buffer-size fft1)))
    (if (= (analysis-output-buffer-size fft2)
           num-fft-bins)
        (progn
          (dotimes (i (ash num-fft-bins -1))
            (let* ((cidx (ash i 1))
                   (x (smp-ref (analysis-output-buffer fft1) cidx))
                   (y (smp-ref (analysis-output-buffer fft1) (1+ cidx)))
                   (u (smp-ref (analysis-output-buffer fft2) cidx))
                   (v (smp-ref (analysis-output-buffer fft2) (1+ cidx)))
                   r-tmp i-tmp)
              (setf i-tmp (+ (* x v) (* y u)))
              (setf r-tmp (- (* x u) (* y v)))
;;;              (break "~&~d: ~a ~a, ~a ~a, ~a ~a" i x y u v r-tmp i-tmp)
              (setf (smp-ref (analysis-output-buffer fft2) cidx) r-tmp)
              (setf (smp-ref (analysis-output-buffer fft2) (1+ cidx)) i-tmp)))
          fft2)
        (error "fft sizes don't match: ~a ~a." fft1 fft2))))

(defun fft-cc-mult (fft)
  "in place complex multiplication of fft output and its complex
conjugate. To calculate the cross-correlation take the ifft and
average the values of the time domain signal. If using the fft of a
sine sweep into some system, this functions calculates the fft of the
impulse response of the system."
  (let ((num-fft-bins (analysis-output-buffer-size fft)))
    (progn
      (dotimes (i (ash num-fft-bins -1))
        (let* ((cidx (ash i 1))
               (real (smp-ref (analysis-output-buffer fft) cidx))
               (img (smp-ref (analysis-output-buffer fft) (1+ cidx))))
          (setf (smp-ref (analysis-output-buffer fft) cidx) (+ (* real real) (* img img)))
          (setf (smp-ref (analysis-output-buffer fft) (1+ cidx)) 0.0d0)))
      fft)))

(defun unzip-fft (fft &optional (dir :input))
  "unzip complex values into left and right half."
  (let* ((vector-size (case dir
                 (:output (analysis-output-buffer-size fft))
                 (otherwise (analysis-input-buffer-size fft))))
         (vector (make-array vector-size))
         (offset (ash vector-size -1))
         (buf (case dir
                (:output (analysis-output-buffer fft))
                (otherwise (analysis-input-buffer fft))))
         )
    (dotimes (i (ash vector-size -1))
      (setf (aref vector i)
            (incudine:smp-ref buf (ash i 1)))
      (setf (aref vector (+ offset i))
            (incudine:smp-ref buf (1+ (ash i 1)))))
    vector))

(defun convolve (buf1 buf2)
  (let* ((result-size
           (+ (buffer-frames buf1)
              (buffer-frames buf2)))
         (fft-size (* 2 (next-power-of-two result-size)))
         (result (incudine:make-buffer result-size))
         (fft1 (make-fft fft-size :flags +fft-plan-fast+))
         (fft2 (make-fft fft-size :flags +fft-plan-fast+))
         (ifft (make-ifft fft-size :flags +fft-plan-fast+)))
    (fft-output->buffer
     (compute-ifft
      (fft-output->ifft-input-norm
       (fft-mult
        (compute-fft (buffer->fft-input buf1 fft1) t)
        (compute-fft (buffer->fft-input buf2 fft2) t))
       ifft)
      nil t)
     result)))


(defun convolve2 (buf1 buf2 fft1 fft2 ifft)
  (let* ((result-size
           (+ (buffer-frames buf1)
              (buffer-frames buf2)))
         (result (incudine:make-buffer result-size)))
    (fft-output->buffer
     (compute-ifft
      (fft-output->ifft-input-norm
       (fft-mult
        (compute-fft (buffer->fft-input buf1 fft1) t)
        (compute-fft (buffer->fft-input buf2 fft2) t))
       ifft)
      nil t)
     result)))

(defun calc-ir (buf orig-rev-fft)
  (let* ((result-size (* 2 (buffer-frames buf)))
         (fft-size (* 2 (next-power-of-two result-size)))
         (result (incudine:make-buffer result-size))
         (fft (make-fft fft-size :flags +fft-plan-fast+))
         (ifft (make-ifft fft-size :flags +fft-plan-fast+)))
    (fft-output->buffer
     (compute-ifft
      (fft-output->ifft-input-norm
       (fft-mult
        (compute-fft (buffer->fft-input buf fft) t)
        orig-rev-fft)
       ifft)
      nil t)
     result)))

(defun cross-correlate (buf)
  (if buf
      (let* ((result-size (ash (buffer-frames buf) 1))
             (fft-size (ash (next-power-of-two result-size) 1))
             (result (incudine:make-buffer result-size))
             (fft (make-fft fft-size :flags +fft-plan-fast+))
             (ifft (make-ifft fft-size :flags +fft-plan-fast+)))
        (fft-output->buffer
         (compute-ifft
          (fft-output->ifft-input-norm
           (fft-cc-mult
            (compute-fft (buffer->fft-input buf fft) t))
           ifft)
          nil t)
         result))))
