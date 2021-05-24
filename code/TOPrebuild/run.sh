#! /bin/bash

#make clean
#make

##background 
./DircRec data/24.00_45.00_2.00_pion.root result/24.00_45.00_2.00_pion.root
./DircRec data/24.00_45.00_2.00_kaon.root result/24.00_45.00_2.00_kaon.root
./DircRec data/27.00_45.00_2.00_pion.root result/27.00_45.00_2.00_pion.root
./DircRec data/27.00_45.00_2.00_kaon.root result/27.00_45.00_2.00_kaon.root
./DircRec data/30.00_45.00_2.00_pion.root result/30.00_45.00_2.00_pion.root
./DircRec data/30.00_45.00_2.00_kaon.root result/30.00_45.00_2.00_kaon.root
./DircRec data/33.00_45.00_2.00_pion.root result/33.00_45.00_2.00_pion.root
./DircRec data/33.00_45.00_2.00_kaon.root result/33.00_45.00_2.00_kaon.root
./DircRec data/24.00_45.00_2.00_bkg_pion.root result/24.00_45.00_2.00_bkg_pion.root
./DircRec data/24.00_45.00_2.00_bkg_kaon.root result/24.00_45.00_2.00_bkg_kaon.root
./DircSep result/pi_K_Background.root 24.00_45.00_2.00 27.00_45.00_2.00 30.00_45.00_2.00 33.00_45.00_2.00 24.00_45.00_2.00_bkg
