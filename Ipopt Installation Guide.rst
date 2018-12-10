Example Installation on Ubuntu 18.04 With Custom Compiled IPOPT
--------------------------------------------------

Install system wide dependencies:

  $ sudo apt install pkg-config python-dev wget
  $ sudo apt build-dep coinor-libipopt1v5

Install pip so all Python packages can be installed via pip:

$ sudo apt install python-pip

Then use pip to install the following packages:

$ pip install --user numpy cython six future

Compile Ipopt

The Ipopt compilation instructions are derived from https://www.coin-or.org/Ipopt/documentation/node14.html. If you get errors, start there for help.

Download Ipopt source code. Choose the version that you would like to have from <https://www.coin-or.org/download/source/Ipopt/>. For example:

$ cd ~
$ wget https://www.coin-or.org/download/source/Ipopt/Ipopt-3.12.11.tgz

Extract the Ipopt source code:

$ tar -xvf Ipopt-3.12.11.tgz

Create a temporary environment variable pointing to the Ipopt directory:

export IPOPTDIR=~/Ipopt-3.12.11

To use linear solvers other than the default mumps, e.g. ma27, ma57, ma86 solvers, the HSL package are needed. HSL can be downloaded from its official website <http://www.hsl.rl.ac.uk/ipopt/>.

Extract HSL source code after you get it. Rename the extracted folder to coinhsl and copy it in the HSL folder: Ipopt-3.12.11/ThirdParty/HSL

Build Ipopt:

$ mkdir $IPOPTDIR/build
$ cd $IPOPTDIR/build
$ ../configure
$ make
$ make test

Add make install if you want a system wide install.

Set environment variables:

$ export IPOPT_PATH="~/Ipopt-3.12.11/build"
$ export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$IPOPT_PATH/lib/pkgconfig
$ export PATH=$PATH:$IPOPT_PATH/bin

Get help from this web-page if you get errors in setting environments:

https://stackoverflow.com/questions/13428910/how-to-set-the-environmental-variable-ld-library-path-in-linux

Now compile cyipopt. Download the cyipopt source code from PyPi, for example:

$ cd ~
$ wget https://files.pythonhosted.org/packages/05/57/a7c5a86a8f899c5c109f30b8cdb278b64c43bd2ea04172cbfed721a98fac/ipopt-0.1.9.tar.gz
$ tar -xvf ipopt-0.1.8.tar.gz
$ cd ipopt

Compile cyipopt:

$ python setup.py build

If there is no error, then you have compiled cyipopt successfully

Check that everything linked correctly with ldd

$ ldd build/lib.linux-x86_64-2.7/cyipopt.so
linux-vdso.so.1 (0x00007ffe895e1000)
libipopt.so.1 => /home/<username>/Ipopt-3.12.11/build/lib/libipopt.so.1 (0x00007f74efc2a000)
libc.so.6 => /lib/x86_64-linux-gnu/libc.so.6 (0x00007f74ef839000)
libcoinmumps.so.1 => /home/<username>/Ipopt-3.12.11/build/lib/libcoinmumps.so.1 (0x00007f74ef4ae000)
libcoinhsl.so.1 => /home/<username>/Ipopt-3.12.11/build/lib/libcoinhsl.so.1 (0x00007f74ef169000)
liblapack.so.3 => /usr/lib/x86_64-linux-gnu/liblapack.so.3 (0x00007f74ee8cb000)
libblas.so.3 => /usr/lib/x86_64-linux-gnu/libblas.so.3 (0x00007f74ee65e000)
libdl.so.2 => /lib/x86_64-linux-gnu/libdl.so.2 (0x00007f74ee45a000)
libstdc++.so.6 => /usr/lib/x86_64-linux-gnu/libstdc++.so.6 (0x00007f74ee0d1000)
libm.so.6 => /lib/x86_64-linux-gnu/libm.so.6 (0x00007f74edd33000)
/lib64/ld-linux-x86-64.so.2 (0x00007f74f02c0000)
libgcc_s.so.1 => /lib/x86_64-linux-gnu/libgcc_s.so.1 (0x00007f74edb1b000)
libcoinmetis.so.1 => /home/<username>/Ipopt-3.12.11/build/lib/libcoinmetis.so.1 (0x00007f74ed8ca000)
libgfortran.so.4 => /usr/lib/x86_64-linux-gnu/libgfortran.so.4 (0x00007f74ed4eb000)

Install cyipopt (prepend sudo if you want a system wide install):

$ python setup.py install

To use cyipopt you will need to set the LD_LIBRARY_PATH to point to your Ipopt install if you did not install it to a standard location. For example:

$ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/Ipopt-3.12.11/build/lib

You can add this to your shell's configuration file if you want it set every time you open your shell, for example the following line can it can be added to your ~/.bashrc

$ echo 'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/Ipopt-3.12.11/build/lib' >> ~/.bashrc

Now you should be able to run a cyipopt example:

$ cd test
$ python -c "import ipopt"
$ python examplehs071.py

If it could be run successfully, the optimization will start with the following descriptions:

******************************************************************************
This program contains Ipopt, a library for large-scale nonlinear optimization.
 Ipopt is released as open source code under the Eclipse Public License (EPL).
         For more information visit http://projects.coin-or.org/Ipopt
******************************************************************************

This is Ipopt version 3.12.11, running with linear solver ma27.
...

