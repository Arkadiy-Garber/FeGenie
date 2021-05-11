#!/bin/bash

docker run -it -v $(pwd):/data --env iron_hmms=/data/hmms/iron --env rscripts=/data/rscripts note/fegenie-deps /bin/bash /data/ci/container-test.sh
