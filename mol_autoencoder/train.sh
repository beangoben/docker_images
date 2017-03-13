#!/bin/bash
source activate python2
cd autoencoder
python train_autoencoder.py \
        ../data/250k_rndm_zinc_drugs_clean.smi \
        ../data/zinc_char_list.json \
        -l5000