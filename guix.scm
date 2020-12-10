;; To use this file to build HEAD of freebayes:
;;
;;   guix build -f guix.scm
;;
;; To get a development container
;;
;;   guix environment -C -l guix.scm

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
    (version (git-version "1.3.3" "HEAD" %git-commit))
    (source (local-file %source-dir #:recursive? #t))
    (build-system meson-build-system)
    (arguments
     `(#:phases
       (modify-phases %standard-phases
         (add-after 'unpack 'prepare-build
           (lambda _
             ;; Stash our build version in the executable.
             (substitute* "src/version_release.txt"
               (("v1.0.0") ,version))
             #t))
         (add-before 'install 'patch-install-location
           (lambda* (#:key outputs #:allow-other-keys)
             (let ((out (assoc-ref outputs "out")))
               (mkdir-p (string-append out "/bin"))
               #t))))
         ; (replace 'check
         ;   (lambda* (#:key outputs #:allow-other-keys)
         ;     (with-directory-excursion "../source/test"
         ;       (invoke "prove" "-e" "bash" "t/00_region_and_target_handling.t")
         ;      (invoke "prove" "-e" "bash" "t/01_call_variants.t")
         ;      (invoke "prove" "-e" "bash" "t/02_multi_bam.t")
         ;      (invoke "prove" "-e" "bash" "t/03_reference_bases.t"))
         ;    #t)))
       #:tests? #f))
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
       ;; ("clang" ,clang)      ; add this to test clang builds
       ;; ("lld" ,lld)          ; add this to test clang builds
       ("bc" ,bc)               ; for tests
       ("coreutils" ,coreutils) ; for echo in test
       ("curl" ,curl)
       ("htslib" ,htslib)
       ("perl6-tap-harness" ,perl6-tap-harness) ; for tests
       ;; ("vcflib" ,vcflib)       ; no longer for tests
       ("zlib" ,zlib)
       ;; ("lz4" ,lz4) ; not used for CRAM
       ("xz" ,xz) ; for liblzma
       ("bzip2" ,bzip2) ; for libz2
       ))
     (synopsis "freebayes haplotype-based genetic variant caller")
     (description
      "freebayes is a Bayesian genetic variant detector designed to find small
polymorphisms, specifically SNPs (single-nucleotide polymorphisms), indels
(insertions and deletions), MNPs (multi-nucleotide polymorphisms), and complex
events (composite insertion and substitution events) smaller than the length of
a short-read sequencing alignment.")
     (home-page "https://github.com/ekg/freebayes")
     (license license:expat)))

freebayes-git
