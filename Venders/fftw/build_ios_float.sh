#!/bin/sh

# build_ios_float.sh
# build an armv7 / arm64 / i386 / x86_64 lib of fftw3
# make sure to check that all the paths used in this script exist on your system
#
# adopted from:
# http://robertcarlsen.net/2009/07/15/cross-compiling-for-iphone-dev-884
# changed by Nickun
# original:
# http://stackoverflow.com/questions/3588904/how-to-link-third-party-libraries-like-fftw3-and-sndfile-to-an-iphone-project-in
#
# changed by ibireme for xcode6 and arm64

export IOS_SDK_VERSION=8.1
export OSX_SDK_VERSION=10.10

# make sure we start out clean
make distclean

# this is the folder where the results of our efforts will end up:
export RESULT_DIR=ios-library

# Select toolchains folder
export XCODE_TOOLCHAINS=/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain
# Select the desired iPhone SDK
export DEVROOT_IOS=/Applications/Xcode.app/Contents/Developer/Platforms/iPhoneOS.platform/Developer
export SDKROOT_IOS=$DEVROOT_IOS/SDKs/iPhoneOS$IOS_SDK_VERSION.sdk
# Select the OSX SDK
export DEVROOT_OSX=/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer
export SDKROOT_OSX=$DEVROOT_OSX/SDKs/MacOSX$OSX_SDK_VERSION.sdk

# ------------------------ armv7---------------------------
# Set up relevant environment variables
export CPPFLAGS="-I$SDKROOT_IOS/usr/include/ -mfpu=neon"
export CFLAGS="$CPPFLAGS -arch armv7 -mfpu=neon -no-cpp-precomp -miphoneos-version-min=5.0 -isysroot $SDKROOT_IOS"
export LD=$XCODE_TOOLCHAINS/usr/bin/ld
export CXX="$XCODE_TOOLCHAINS/usr/bin/clang -O3 -x c++ -arch armv7 -std=gnu++11 -stdlib=libc++ -mfpu=neon"
export CC="$XCODE_TOOLCHAINS/usr/bin/clang -O3 -x c -arch armv7 -std=gnu99 -mfpu=neon"
export CXXFLAGS="$CFLAGS"

# TODO: add custom flags as necessary for package
#  remove '--enable-float' for double precision
#  and take a 'libfftw3.a' file instead
./configure --host=arm-apple-darwin --enable-float --enable-neon

make -j2

# Copy the ARM library to a temporary location
mkdir $RESULT_DIR
cp .libs/libfftw3f.a $RESULT_DIR/libfftw3f_armv7.a

# Copy the header file too, just for convenience
cp api/fftw3.h $RESULT_DIR/fftw3.h

# ------------------------ arm64---------------------------
# Do it all again for i386
make distclean

# Restore default environment variables
unset CPPFLAGS CFLAGS CPP LD CXX LDFLAGS CXXFLAGS

# Set up relevant environment variables
export CPPFLAGS="-I$SDKROOT_IOS/usr/include/ -mfpu=neon"
export CFLAGS="$CPPFLAGS -arch arm64 -mfpu=neon -no-cpp-precomp -miphoneos-version-min=5.0 -isysroot $SDKROOT_IOS"
export LD=$XCODE_TOOLCHAINS/usr/bin/ld
export CXX="$XCODE_TOOLCHAINS/usr/bin/clang -O3 -x c++ -arch arm64 -std=gnu++11 -stdlib=libc++ -mfpu=neon"
export CC="$XCODE_TOOLCHAINS/usr/bin/clang -O3 -x c -arch arm64 -std=gnu99 -mfpu=neon"
export CXXFLAGS="$CFLAGS"

# TODO: add custom flags as necessary for package
#  remove '--enable-float' for double precision
#  and take a 'libfftw3.a' file instead
./configure --host=arm-apple-darwin --enable-float --enable-neon

make -j2

# Copy the ARM library to a temporary location
mkdir $RESULT_DIR
cp .libs/libfftw3f.a $RESULT_DIR/libfftw3f_arm64.a


# ------------------------ i386 ---------------------------
# Do it all again for i386
make distclean

# Restore default environment variables
unset CPPFLAGS CFLAGS CPP LD CXX LDFLAGS CXXFLAGS

export CPPFLAGS="-I$SDKROOT_OSX/usr/include/"
export CFLAGS="$CPPFLAGS -arch i386 -no-cpp-precomp -mmacosx-version-min=10.7 -isysroot $SDKROOT_OSX"
export LD="$XCODE_TOOLCHAINS/usr/bin/ld"
export CXX="$XCODE_TOOLCHAINS/usr/bin/clang -O3 -x c++ -arch i386 -std=gnu++11 -stdlib=libc++"
export CC="$XCODE_TOOLCHAINS/usr/bin/clang -O3 -x c -arch i386 -std=gnu99"
export CXXFLAGS="$CFLAGS"

# TODO: error checking
./configure --enable-float
make -j2

# Copy the FAT native library to the temporary location
cp .libs/libfftw3f.a $RESULT_DIR/libfftw3f_i386.a

# ------------------------ x86_64 ---------------------------
# Do it all again for x86_64
make distclean

# Restore default environment variables
unset CPPFLAGS CFLAGS CPP LD CXX LDFLAGS CXXFLAGS

export CPPFLAGS="-I$SDKROOT_OSX/usr/include/"
export CFLAGS="$CPPFLAGS -arch x86_64 -no-cpp-precomp -mmacosx-version-min=10.7 -isysroot $SDKROOT_OSX"
export LD="$XCODE_TOOLCHAINS/usr/bin/ld"
export CXX="$XCODE_TOOLCHAINS/usr/bin/clang -O3 -x c++ -arch x86_64 -std=gnu++11 -stdlib=libc++"
export CC="$XCODE_TOOLCHAINS/usr/bin/clang -O3 -x c -arch x86_64 -std=gnu99"
export CXXFLAGS="$CFLAGS"

# TODO: error checking
./configure --enable-float
make -j2

# Copy the FAT native library to the temporary location
cp .libs/libfftw3f.a $RESULT_DIR/libfftw3f_x86_64.a



# Create fat lib by combining the two versions
lipo -arch armv7 $RESULT_DIR/libfftw3f_armv7.a -arch arm64 $RESULT_DIR/libfftw3f_arm64.a -arch i386 $RESULT_DIR/libfftw3f_i386.a -arch x86_64 $RESULT_DIR/libfftw3f_x86_64.a -create -output $RESULT_DIR/libfftw3_float.a

# Remove intermediate binaries
#rm $RESULT_DIR/libfftw3_armv7.a
#rm $RESULT_DIR/libfftw3_i386.a
#rm $RESULT_DIR/libfftw3_x86_64.a

# Unset used environment variables
unset CPPFLAGS CFLAGS CPP LD LDFLAGS CXX CXXFLAGS