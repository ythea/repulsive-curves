# Ubuntu 20.04: clang-10
# Ubuntu 22.04: clang-11/12/13 ~
CLANG_VER=18
#sudo apt update
#sudo apt install libc++-$CLANG_VER-dev
#sudo apt install libc++abi-$CLANG_VER-dev
#sudo apt install libomp-$CLANG_VER-dev
#git pull --recursive origin master
#git submodule update --init --recursive
mkdir -p build
pushd build
export CC=clang-18
export CXX=clang++-18
BUILD_TYPE=Debug #Release #RelWithDebInfo #Debug
/home/ducthan/CMAKE/3.30.3/bin/cmake .. \
		 -DCMAKE_INSTALL_PREFIX=../release \
		 -DCMAKE_BUILD_TYPE=$BUILD_TYPE

cmake --build . -j16
cmake --install .
popd
