export PATH=$PATH:$(pwd)

gzip -d HMM-lib.tar.gz
tar xf HMM-lib.tar
rm HMM-lib.tar

gzip -d test_dataset.tar.gz
tar xf test_dataset.tar
rm test_dataset.tar
