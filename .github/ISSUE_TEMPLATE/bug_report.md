---
name: Bug report 🐞
about: Create a report to help us improve
title: ''
labels: bug
assignees: ''
---
**Only bug reports!**

The C++ version of Freebayes is in *maintenance* mode. Use the github issue
tracker to report bugs *only*. For comments, questions and features,
please use the google group mailing list as stated on the
[README](https://github.com/freebayes/freebayes/blob/master/README.md)!

**Describe the bug**

A clear and concise description of what the bug is.

**To Reproduce**

Include all steps to reproduce the behavior and paste any complete
errors from the terminal.

**Expected behavior**

A clear and concise description of what you expected to happen.

**Screenshots**

If applicable, add screenshots to help explain your problem.

**Additional context**

Add any other context about the problem here.

Include a set of BAM/BED files to reproduce the issue

+ bonus points if you try to minimize the test case yourself, as issues are often localized:
  - try to use sambamba or samtools slice to first extract the reference where the error occurs
  - if that succeeds (the error is still reproducible), continue to crop the file in binary-search fashion

**Finally**

Please check the README and docs carefully. Everyone working on freebayes is doing that for free. Please respect our time (too).
