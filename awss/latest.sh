#!/bin/bash
echo "caco-01"
echo "latest:"
aws s3 ls s3://cmgp-coastcam/cameras/caco-01/latest --recursive --profile "coastcam" | sort | tail -n 10
echo "products:"
aws s3 ls s3://cmgp-coastcam/cameras/caco-01/products --recursive --profile "coastcam" | sort | tail -n 12
echo ""
echo "caco-02"
echo "latest:"
aws s3 ls s3://cmgp-coastcam/cameras/caco-02/latest --recursive --profile "coastcam" | sort | tail -n 10
echo "products:"
aws s3 ls s3://cmgp-coastcam/cameras/caco-02/products --recursive --profile "coastcam" | sort | tail -n 10
