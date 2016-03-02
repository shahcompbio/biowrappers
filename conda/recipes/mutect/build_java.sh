#/bin/bash -eu

java_dir=$1/java

URL=http://download.oracle.com/otn-pub/java/jdk/7u79-b15/jdk-7u79-linux-x64.tar.gz

NAME=jdk1.7.0_79

TAR_FILE=jdk-7u79-linux-x64.tar.gz

wget --no-check-certificate --no-cookies --header "Cookie: oraclelicense=accept-securebackup-cookie" $URL -O $TAR_FILE

tar -zxvf $TAR_FILE

rm $TAR_FILE

mkdir -p $java_dir

mv $NAME/jre/* $java_dir

rm -rf $NAME

rm -rf $java_dir/*.zip 
