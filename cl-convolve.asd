;;;; cl-convolve.asd

(asdf:defsystem #:cl-convolve
  :description "Describe cl-convolve here"
  :author "Orm Finnendahl <orm.finnendahl@selma.hfmdk-frankfurt.de>"
  :license  "gpl 2.0 or higher"
  :version "0.0.1"
  :serial t
  :depends-on (:incudine :incudine-plot)
  :components ((:file "package")
               (:file "cl-convolve")))
