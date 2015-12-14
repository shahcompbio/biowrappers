out_dir=$PREFIX/share/$PKG_NAME-$PKG_VERSION-$PKG_BUILDNUM

$PYTHON setup.py install

feature_file=normal_tumour_features_$PKG_VERSION.h5

rsync -auv monco:/home/aroth/software/museq/${feature_file} ${feature_file}

mkdir -p $out_dir

museq train_model --in_file ${feature_file} --out_file ${out_dir}/normal_tumour_model.pickle