#! /usr/bin/env bash
# Author: Izaak Neutelings (October 2025)
# Description: Produce the plots for the DPS note
INFILE=output/tauAlgoEffs_2024_DY.root
python3 plotTauROCCurves.py $INFILE -m roc -t _DY -v1
python3 plotTauPt.py $INFILE -v1
python3 plotTauIDVennDiagrams.py $INFILE -v1