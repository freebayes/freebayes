;; To use this file to build HEAD of freebayes:
;; guix build -f guix.scm
;;
;; TODO:
;; Read previous version instead of hardcoding it
;; Run the test suite

(use-modules
  ((guix licenses) #:prefix license:)
  (guix gexp)
  (guix packages)
  (guix git-download)
  (guix build-system cmake)
  (gnu packages base)
  (gnu packages compression)
  (gnu packages bioinformatics)
  (gnu packages build-tools)
  (gnu packages curl)
  (gnu packages ninja)
  (gnu packages perl)
  (gnu packages pkg-config)
  (srfi srfi-1)
  (ice-9 popen)
  (ice-9 rdelim))

(define %source-dir (dirname (current-filename)))

(define %git-commit
    (read-string (open-pipe "git show HEAD | head -1 | cut -d ' ' -f 2" OPEN_READ)))

(define-public freebayes-git
  (package
    (inherit freebayes)
    (name "freebayes-git")
    (version (git-version "1.3.2" "HEAD" %git-commit))
    (source (local-file %source-dir #:recursive? #t))
    (build-system cmake-build-system)
    (arguments
     `(#:phases
       (modify-phases %standard-phases
         (add-after 'unpack 'prepare-build
           (lambda _
             ;; Set SHELL for SeqLib's configure script.
             (setenv "CONFIG_SHELL" (which "sh"))
             ;; Stash our build version in the executable.
             (substitute* "src/version_release.txt"
               (("v1.0.0") ,version))
             #t))
         (add-before 'install 'patch-install-location
           (lambda* (#:key outputs #:allow-other-keys)
             (let ((out (assoc-ref outputs "out")))
               (substitute* "Makefile"
                 (("/usr/local") out))
               (mkdir-p (string-append out "/bin"))
               #t))))
       #:tests? #f
       #:make-flags (list "CC=gcc")))
    (propagated-inputs
     `(("perl" ,perl)         ; for testing
       ("grep" ,grep)         ; for testing
       ("samtools" ,samtools) ; for testing
       ("which" ,which)       ; for version
       ))
    (native-inputs
     `(
       ("meson" ,meson)
       ("ninja" ,ninja)
       ("pkg-config" ,pkg-config)
       ))
    (inputs
     `(
       ("curl" ,curl)
       ("zlib" ,zlib)
       ))))

freebayes-git
