#!/bin/bash

docker run -it -v $(pwd):/data fegenie /bin/bash /data/ci/container-test.sh
