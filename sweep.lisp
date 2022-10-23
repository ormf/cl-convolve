;;; 
;;; sweep.lisp
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

(in-package :cl-convolve)

(defun mtof (m)
  (* 440.0 (expt 2 (/ (- m 69) 12))))

(defun ftom (f)
  (+ 69 (* 12 (log (/ f 440) 2))))

(defparameter *sweep-env* (make-envelope '(0 1 1 0) '(0.05 .9 .05)))

(dsp! sweep-log (start end dur amp)
  (reduce-warnings
    (out (* (envelope *sweep-env* 1 dur :done-action #'free)
            (sine (mtof (line (ftom start) (ftom end) dur :done-action #'free))
                  amp 0)))))

(dsp! sweep-lin (start end dur amp)
  (out (* (envelope *sweep-env* 1 dur :done-action #'free)
          (sine (line start end dur)
                amp 0))))


;;; (sweep-log 50 20000 5 0.2)

