#!/bin/bash

outdir=$PREFIX/share/$PKG_NAME-$PKG_VERSION-$PKG_BUILDNUM

mkdir -p $outdir

mkdir -p $PREFIX/bin

cp -R ./* $outdir/

cp $RECIPE_DIR/snpeff.py $outdir/snpeff.py

chmod +x $outdir/snpeff.py

ln -s $outdir/snpeff.py $PREFIX/bin/bw-snpeff

cd $outdir

bw-snpeff download GRCh37.75
