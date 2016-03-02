outdir=$PREFIX/share/$PKG_NAME-$PKG_VERSION-$PKG_BUILDNUM

./configure --prefix=$outdir

make

make install

ln -s $outdir/libexec/countFastaBases $PREFIX/bin

ln -s $outdir/libexec/strelka2 $PREFIX/bin
