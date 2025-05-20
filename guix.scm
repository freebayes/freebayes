;; To use this file to build HEAD of freebayes:
;;
;;   guix build -f guix.scm
;;
;; To get a development container
;;
;;   guix shell -C -D -F -f guix.scm
;;
;; Build with
;;
;;   meson setup build/ --buildtype debug --wipe
;;   cd build
;;   ninja
;;   meson test -t 2
;;
;; For a static build use
;;
;;   guix build -L . freebayes-static-git [--tune=native]
;;   meson setup build/ --buildtype debug --wipe -Dprefer_system_deps=false -Dstatic=true
;;
;; and for debugging
;;
;;   guix shell -C -D -F -L . freebayes-debug

(define-module (guix)
  #:use-module (ice-9 popen)
  #:use-module (ice-9 rdelim)
  #:use-module (srfi srfi-1)
  #:use-module ((guix licenses) #:prefix license:)
  #:use-module (guix build-system cmake) ; for vcflib
  #:use-module (guix build-system meson)
  #:use-module (guix download)
  #:use-module (guix gexp)
  #:use-module (guix git-download)
  #:use-module (guix packages)
  #:use-module (guix utils)
  #:use-module (gnu packages algebra)
  #:use-module (gnu packages assembly)
  #:use-module (gnu packages base)
  #:use-module (gnu packages bioinformatics)
  #:use-module (gnu packages build-tools)
  #:use-module (gnu packages check)
  #:use-module (gnu packages compression)
  #:use-module (gnu packages curl)
  #:use-module (gnu packages gcc)
  #:use-module (gnu packages gdb)
  #:use-module (gnu packages haskell-xyz) ; pandoc for help files
  #:use-module (gnu packages llvm)
  #:use-module (gnu packages ninja)
  #:use-module (gnu packages parallel)
  #:use-module (gnu packages perl)
  #:use-module (gnu packages perl6)
  #:use-module (gnu packages pkg-config)
  #:use-module (gnu packages python)
  #:use-module (gnu packages python-xyz) ; for pybind11
  #:use-module (gnu packages ruby)
  #:use-module (gnu packages time)
  #:use-module (gnu packages tls)
  #:use-module (gnu packages zig))

(define %source-dir (dirname (current-filename)))

(define %git-commit
    (read-string (open-pipe "git show HEAD | head -1 | cut -d ' ' -f 2" OPEN_READ)))

(define-public vcflib-github ;; should update upstream
  (let ((commit
         "9e8c0192d677bddd68eb2743ff9ec194995e92b7"
                ))
    (package
     (name "vcflib-github")
     (version (string-append "1.0.14-" (string-take commit 7)))
     (source
      (origin
       (method git-fetch)
       (uri (git-reference
             (url "https://github.com/vcflib/vcflib/")
             (commit commit)
             (recursive? #t)))
       (file-name (string-append name "-" version "-checkout"))
       (sha256
        (base32
         "0jvy98q4zy7md14c2h40xpd9hng8a2070g8b3lw0rigjnfkwziyl"
        ))))
    (build-system cmake-build-system)
    (arguments
     `(#:tests? #f
       #:configure-flags
       ,#~(list
           "-DCMAKE_BUILD_TYPE=Debug"
           "-DZIG=OFF")))
    (inputs
     (list
       curl
       fastahack  ;; dev version not in Debian
       htslib ;; disable to test local build - with local package below
       pandoc ; for man pages
       perl
       python
       python-pytest
       pybind11
       ruby ; for man pages
       smithwaterman
       tabixpp
       time ; for tests
       wfa2-lib ; alternative:  cmake  -DCMAKE_BUILD_TYPE=Debug -DWFA_GITMODULE=ON -DZIG=ON ..
       xz
       ;; zig-0.14
       ))
    (native-inputs
     `(("pkg-config" ,pkg-config)))
    (home-page "https://github.com/vcflib/vcflib/")
    (synopsis "Library for parsing and manipulating VCF files")
    (description "Vcflib")
    (license license:expat))))

;; Guix does not come with a static version of libdeflate
(define-public libdeflate-static
  (package
    (inherit libdeflate)
    (name "libdeflate-static")
    (version "1.19")
    (arguments
     (list #:configure-flags
           #~(list "-DLIBDEFLATE_BUILD_STATIC_LIB=YES"
                   "-DLIBDEFLATE_BUILD_TESTS=YES")))))

;; A minimal static version of htslib that does not depend on curl and openssl. This
;; reduces the number of higher order dependencies in static linking.
(define-public htslib-static
  (package
    (inherit htslib)
    (name "htslib-static")
    (version "1.19")
    (source (origin
            (method url-fetch)
            (uri (string-append
                  "https://github.com/samtools/htslib/releases/download/"
                  version "/htslib-" version ".tar.bz2"))
            (sha256
             (base32
              "0dh79lwpspwwfbkmllrrhbk8nkvlfc5b5ib4d0xg5ld79w6c8lc7"))))
    (arguments
     (substitute-keyword-arguments (package-arguments htslib)
       ((#:configure-flags flags ''())
        ''())))
    (inputs
     (list bzip2 xz))))


(define-public vcflib-static-github
  "Optimized for latest AMD architecture build and static deployment. These binaries can be copied to HPC."
  (package
    (inherit vcflib-github)
    (name "vcflib-static-github")
    (arguments
     `(#:tests? #f
       #:configure-flags
       ,#~(list
           "-DBUILD_STATIC=ON"
           "-DZIG=OFF"
           "-DCMAKE_BUILD_TYPE=Generic" ;; to optimize use guix -- tune=march-type (e.g. --tune=native)
           "-DCMAKE_INSTALL_RPATH=")))   ; force cmake static build and do not rewrite RPATH
    (inputs
     (modify-inputs (package-inputs vcflib-github)
                    (prepend
                     `(,bzip2 "static")
                     `(,zlib "static")
                     `(,xz "static")
                     libdeflate-static
                     htslib-static)))))


(define-public freebayes-git
  (package
    (name "freebayes-git")
    (version (git-version "1.3.10" "HEAD" %git-commit))
    (source (local-file %source-dir #:recursive? #t))
    (build-system meson-build-system)
    (propagated-inputs
     `(
       ;; for the libs also see contrib/README.md
       ("bzip2-static" ,bzip2 "static")    ; libz2 part of htslib for static builds
       ("fastahack" ,fastahack) ; shared lib used by vcflib; bundle for Debian
       ("grep" ,grep)          ; for testing
       ("htslib" ,htslib)      ; does work, but lacks codecs
       ("intervaltree" ,intervaltree)
       ("perl" ,perl)          ; for testing
       ("python" ,python)      ; for testing
       ("samtools" ,samtools)  ; for testing
       ("simde" ,simde)
       ("smithwaterman" ,smithwaterman) ; vcflib shared lib dependency ; bundle for Debian
       ("tabixpp" ,tabixpp)    ; for htslib
       ("vcflib-github" ,vcflib-github)  ; for includes and testing freebayes-parallel
       ("wfa2-lib" ,wfa2-lib)  ; vcflib shared lib dependency
       ("which" ,which)))        ; for version
    (native-inputs
     `(
       ("meson" ,meson)
       ("ninja" ,ninja)
       ("gdb" ,gdb)
       ("pkg-config" ,pkg-config)
       ))
    (inputs
     `(
       ;; ("clang" ,clang)      ; add this to test clang builds
       ;; ("lld" ,lld)          ; add this to test clang builds
       ;; ("gcc" ,gcc-13)
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
     (list #:phases
           #~(modify-phases %standard-phases
                  (add-after 'unpack 'includes
                    (lambda _
                      (substitute* "meson.build"
                                   (("vcflib_inc = files\\(\\)")
                                    (string-append "vcflib_inc = include_directories('" #$vcflib-github "/include/vcflib')")))))
                  ;; add timeout extension for slower processor
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

(define-public freebayes-static-git
  "Optimized for latest AMD architecture build and static deployment. These binaries can be copied to HPC."
  (package
    (inherit freebayes-git)
    (name "freebayes-static-git")
    (arguments
     (list
      #:tests? #f
      #:configure-flags
      #~(list
         "-Dprefer_system_deps=false"
         "-Dstatic=true")  ; force static build and do not rewrite RPATH
      #:phases
        #~(modify-phases %standard-phases
                  (add-after 'unpack 'includes
                    (lambda _
                      (substitute* "meson.build"
                                   (("vcflib_inc = files\\(\\)")
                                    (string-append "vcflib_inc = include_directories('" #$vcflib-static-github "/include/vcflib')"))))))))
    (inputs
     (modify-inputs (package-inputs freebayes-git)
                    (prepend
                     `(,bzip2 "static")
                     `(,zlib "static")
                     `(,xz "static")
                     ;; ("xz-static" ,xz "static")     ; for static builds
                     ;; ("zlib-static" ,zlib "static")))
                     libdeflate-static
                     vcflib-static-github
                     htslib-static)))))

(define-public freebayes-debug
  (package
    (inherit freebayes-static-git)
    (name "freebayes-debug")
    (build-system meson-build-system)
    ;; (inputs
    ;;  (modify-inputs (package-inputs freebayes-git)
    ;;                 (prepend
    ;;                  meson)))
    (arguments
     `(#:tests? #f
       #:phases (modify-phases %standard-phases
                               (delete 'configure)
                               (delete 'build)
                               (delete 'package)
                               (delete 'check)
                               (delete 'install))))))

freebayes-debug
