;; To use this file to build HEAD of freebayes:
;;
;;   guix build -f guix.scm
;;
;; To get a development container
;;
;;   guix environment -C -l guix.scm
;;
;; For the tests you need /usr/bin/env. In a container create it with
;;
;;   mkdir -p /usr/bin
;;   ln -s $GUIX_ENVIRONMENT/bin/env /usr/bin/env

(use-modules
  ((guix licenses) #:prefix license:)
  (guix gexp)
  (guix packages)
  (guix git-download)
  (guix build-system meson)
  (gnu packages algebra)
  (gnu packages base)
  (gnu packages compression)
  (gnu packages bioinformatics)
  (gnu packages build-tools)
  (gnu packages curl)
  (gnu packages llvm)
  (gnu packages ninja)
  (gnu packages parallel)
  (gnu packages perl)
  (gnu packages perl6)
  (gnu packages pkg-config)
  (srfi srfi-1)
  (ice-9 popen)
  (ice-9 rdelim))

(define %source-dir (dirname (current-filename)))

(define %git-commit
    (read-string (open-pipe "git show HEAD | head -1 | cut -d ' ' -f 2" OPEN_READ)))

(define-public freebayes-git
  (package
    (name "freebayes-git")
    (version (git-version "1.3.4" "HEAD" %git-commit))
    (source (local-file %source-dir #:recursive? #t))
    (build-system meson-build-system)
    (propagated-inputs
     `(("perl" ,perl)         ; for testing
       ("grep" ,grep)         ; for testing
       ("samtools" ,samtools) ; for testing
       ("vcflib" ,vcflib)     ; for freebayes-parallel
       ("which" ,which)       ; for version
       ("htslib" ,htslib)     ; does work, but lacks codecs
       ("tabixpp" ,tabixpp)
        ))
    (native-inputs
     `(
       ("meson" ,meson)
       ("ninja" ,ninja)
       ("pkg-config" ,pkg-config)
       ))
    (inputs
     `(
       ;; ("clang" ,clang)      ; add this to test clang builds
       ;; ("lld" ,lld)          ; add this to test clang builds
       ("bc" ,bc)               ; for tests
       ("coreutils" ,coreutils) ; for echo and env in tests
       ("curl" ,curl)
       ("perl6-tap-harness" ,perl6-tap-harness) ; for tests
       ("parallel" ,parallel) ; for freebayes-parallel
       ("zlib" ,zlib)
       ("xz" ,xz)          ; liblzma part of htslib
       ("bzip2" ,bzip2)    ; libz2 part of htslib
       ))
    (arguments
     `(#:phases (modify-phases %standard-phases
       ;; add timeout extension for slower processors
       (replace 'check
                (lambda _
                 (invoke "meson" "test" "--timeout-multiplier" "5"))))))
     (synopsis "freebayes haplotype-based genetic variant caller")
     (description
      "freebayes is a Bayesian genetic variant detector designed to find small
polymorphisms, specifically SNPs (single-nucleotide polymorphisms), indels
(insertions and deletions), MNPs (multi-nucleotide polymorphisms), and complex
events (composite insertion and substitution events) smaller than the length of
a short-read sequencing alignment.")
     (home-page "https://github.com/freebayes/freebayes")
     (license license:expat)))

freebayes-git
