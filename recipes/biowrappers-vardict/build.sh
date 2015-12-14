outdir=$PREFIX/share/$PKG_NAME-$PKG_VERSION-$PKG_BUILDNUM

export TERM=dumb

./gradlew clean installApp

mkdir -p ${outdir}

mv build/install/VarDict/* ${outdir}

mv VarDict ${outdir}/scripts

cp $RECIPE_DIR/vardict.py ${outdir}

chmod +x ${outdir}/vardict.py

mkdir -p $PREFIX/bin

ln -s ${outdir}/vardict.py $PREFIX/bin/bw-vardict