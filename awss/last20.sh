#!/bin/bash
echo "...working"
echo "caco-02"
aws s3 ls s3://cmgp-coastcam/cameras/caco-02/products --recursive --profile "coastcam" | grep "timex" | sort | tail -n 20 > last20.txt

cat last20.txt
exit 0
