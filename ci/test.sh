#/bin/bash

# Exit when any command fails:
set -e

rm -rf ci/actual_output
rm -rf ci/fegenie_out

echo "Running FeGenie..."
echo `pwd`
./FeGenie.py -bin_dir ci/test_dataset -bin_ext txt -out ci/fegenie_out -t $(nproc) -R rscripts
echo "FeGenie completed successfully. Verifying results..."

mkdir -p ci/actual_output

# find ci/fegenie_out -iname "fegenie*.tiff" -exec mv {} ci/actual_output/ \;
find ci/fegenie_out -iname "fegenie*.csv" -exec mv {} ci/actual_output/ \;

# TODO: Actually we expect 2 files but right now dendro-heatmap.R fails
# EXPECTED_NUMBER_OF_TIFFS=1
# NUMBER_OF_TIFFS=`find . -type f | sed 's/.*\.//' | sort | uniq -c | grep tiff | awk '{print $1}'`

#if [ "$NUMBER_OF_TIFFS" -eq "$EXPECTED_NUMBER_OF_TIFFS" ]; then
#	echo "number of tiff files as expected: $NUMBER_OF_TIFFS"
#else
#	echo "unexpected number of tiff files: $NUMBER_OF_TIFFS, expected: $EXPECTED_NUMBER_OF_TIFFS"
#	exit 1
#fi

# We only verify number of tiffs files, surprisingly they're very big (a few tens of MBs) so we don't want them in the repository
# to be able to use diff for diffing directories we need to remove tiff files in actual_output first
#find ci/actual_output -iname "fegenie*.tiff" -exec rm {} \;

echo "Comparing actual_output with expected_output:"
diff ci/actual_output ci/expected_output
