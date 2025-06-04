# Instructions to be able to run this code on LINUX

$ sudo apt update
$ sudo apt install build-essential
$ sudo apt install coinor-libipopt-dev
$ pip3 install casadi

$ sudo apt install manpages-dev software-properties-common
$ sudo add-apt-repository ppa:ubuntu-toolchain-r/test
$ sudo apt update && sudo apt install gcc-13 g++-13
$ sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-13 110
$ sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-13 110

$ sudo apt install gfortran liblapack-dev pkg-config --install-recommends
$ sudo apt install swig
$ cd
$ git clone https://github.com/casadi/casadi.git -b master casadi
$ cd casadi
$ mkdir build
$ cd build
$ cmake -DWITH_PYTHON=ON -DWITH_IPOPT=ON -DWITH_OPENMP=ON -DWITH_THREAD=ON ..
$ make
$ sudo make install

$ cd RocketLab/CppExercise/
$ mkdir build
$ cd build
$ cmake ..
$ make