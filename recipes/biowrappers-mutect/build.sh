outdir=$PREFIX/share/$PKG_NAME-$PKG_VERSION-$PKG_BUILDNUM

sh $RECIPE_DIR/build_java.sh $outdir

# Copy files
mkdir -p $outdir

cp mutect-$PKG_VERSION.jar $outdir/mutect.jar

cp $RECIPE_DIR/mutect.py $outdir/mutect.py

chmod +x $outdir/mutect.py

mkdir -p $PREFIX/bin

ln -s $outdir/mutect.py $PREFIX/bin/bw-mutect
