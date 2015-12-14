out_dir=$PREFIX/share/$PKG_NAME-$PKG_VERSION-$PKG_BUILDNUM

$PYTHON setup.py install --single-version-externally-managed --record=record.txt

feature_file=$RECIPE_DIR/src/normal_tumour_features.h5

mkdir -p $out_dir

museq train_model --in_file ${feature_file} --out_file ${out_dir}/normal_tumour_model.pickle

cp $RECIPE_DIR/museq.py $out_dir/museq.py

chmod +x $out_dir/museq.py

mkdir -p $PREFIX/bin

ln -s $out_dir/museq.py $PREFIX/bin/bw-museq