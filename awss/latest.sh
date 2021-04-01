#!/bin/bash
echo "...working"
echo "caco-02"
echo "latest"
aws s3 ls s3://cmgp-coastcam/cameras/caco-02/latest --recursive --profile "coastcam" | sort | tail -n 10 > latest.txt
echo "products"
aws s3 ls s3://cmgp-coastcam/cameras/caco-02/products --recursive --profile "coastcam" | sort | tail -n 10 >> latest.txt
echo ""
echo "caco-01"
echo "latest"
aws s3 ls s3://cmgp-coastcam/cameras/caco-01/latest --recursive --profile "coastcam" | sort | tail -n 10 >> latest.txt
echo "products (may take a while):"
aws s3 ls s3://cmgp-coastcam/cameras/caco-01/products --recursive --profile "coastcam" | sort | tail -n 12 >> latest.txt

cat latest.txt
exit 0
