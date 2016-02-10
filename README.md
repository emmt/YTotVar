# YTotVar - Edge preserving regularization for Yorick

This plug-in implements relaxed total variation (TV) and edge-preserving
regularization for [Yorick](http://yorick.github.com/) in 2-D, 3-D or 4-D.


## Installation

To compile and install from the source directory:
````~shell
  shell> yorick -batch make.i
  shell> make clean
  shell> make
  shell> make install
````

It is also possible to compile and install from another directory:
````~shell
  shell> mkdir -p ${BUILD_DIR}
  shell> cd ${BUILD_DIR}
  shell> ${SRC_DIR}/configure
  shell> make clean
  shell> make
  shell> make install
````
where `${BUILD_DIR}` is the building directory and `$SRC_DIR}` is the source
directory.  The `configure` script has some options which can be displayed
with:
````~shell
  shell> ${SRC_DIR}/configure --help
````



