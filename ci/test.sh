#!/bin/bash

docker run -it -v $(pwd):/data note/fegenie-deps /bin/bash /data/ci/container-test.sh
