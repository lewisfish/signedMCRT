---
project: signedMCRT
src_dir: ./src
         ./app
output_dir: ./docs
media_dir: ./images
project_github: https://github.com/lewisfish/signedMCRT
project_download: https://github.com/lewisfish/signedMCRT/releases/latest
summary: ![signedMCRT](|media|/sMCRT_logo.png)<br> Monte Carlo radiation transfer (MCRT) using Signed Distance functions (SDF)
author: Lewis McMillan
github: https://github.com/lewisfish
predocmark_alt: >
predocmark: <
docmark_alt:
docmark: !
display: public
         private
         protected
source: true
graph: false
search: false
sort: alpha
fpp_extensions: fpp
preprocess: true
extra_mods: iso_fortran_env:https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html
md_extensions: markdown.extensions.toc
               markdown_checklist.extension
               ford.md_striped_table
page_dir: ./old_docs
---

--------------------
[TOC]

Brief description
-----------------

A Monte Carlo radiation transfer code with signed distance functions representing the geometry, written in modern Fortran.


License
-------

The signedMCRT source code and related files and documentation are
distributed under a permissive free software license (MIT).

