#!/bin/sh 

source activate rh_classifier

#python Rh_classifier-train_test.py

python Rh_classifier-assess_results.py

python Rh_classifier-generate_plots.py