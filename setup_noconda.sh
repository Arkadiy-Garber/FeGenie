export PATH=$PATH:$(pwd)

gzip -d test_dataset.tar.gz
tar xf test_dataset.tar
rm test_dataset.tar
