GUPPY_VER='3.4.5'

echo "INSTALLING GUPPY VERSION ${GUPPY_VER}"
wget https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy_${GUPPY_VER}_linux64.tar.gz 
tar -zvxf ont-guppy_${GUPPY_VER}_linux64.tar.gz
mv ont-guppy mop_preprocess/bin/
cd mop_preprocess/bin
ln -s ont-guppy/bin/guppy_* .
ln -s ont-guppy/lib/* .
cd ../../
if [ ! -e "mop_preprocess/bin/ont-guppy/lib/libz.so" ] ; then
        unlink mop_preprocess/bin/ont-guppy/lib/libz.so
        cd mop_preprocess/bin/ont-guppy/lib/
        ln -s libz.so.1 libz.so
        cd ../../../../
fi
rm ont-guppy_${GUPPY_VER}_linux64.tar.gz

echo 'copying custom models'
for i in custom_models/*.gz; do zcat $i > mop_preprocess/bin/ont-guppy/data/`basename $i .gz`; done


