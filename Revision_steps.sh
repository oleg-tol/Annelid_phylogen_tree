#!/usr/bin/env bash
# Revision analyses (reviewer comments): add two C. teleta r-opsins, remove redundant
# H. robusta entries, re-infer ML + Bayesian trees, AU test for O. fusiformis placement.
# Run from the folder containing opsin_newA4.fa (final published 7TM alignment).
set -euo pipefail
mkdir -p revision_capitella && cd revision_capitella

# 1) Add C. teleta r-opsin2/3 to the existing alignment (columns and existing rows unchanged)
#    capitella_new_coded2.fa: >Annelid_ropsins_23/1-206 (AXB38864.1), >Annelid_ropsins_24/1-391 (AXB38865.1)
mafft --add capitella_new_coded2.fa --keeplength --mapout --localpair --maxiterate 1000 --thread 4 \
      ../opsin_newA4.fa > opsin_newA4_rev.fa

# 2) Remove redundant H. robusta entries (same genes as cDNA + genome model, >=98.6% identity):
#    Annelid_ropsins_10 = XP_009022296.1 (= opsin B AID66634.1), Annelid_ropsins_12 = AID66633.1 (= XP_009024616.1)
awk '/^>/{p=!($0 ~ /^>Annelid_ropsins_(10|12)\//)} p' opsin_newA4_rev.fa > opsin_newA4_rev2.fa   # 52 seqs x 319

# 3) Maximum likelihood (same settings as the original analysis)
iqtree2 -s opsin_newA4_rev2.fa -m LG+F+R8 -bb 2000 -alrt 1000 -nt 4 -pre opsin_newA4_rev2

# 4) Bayesian inference: NEXUS (names without /coordinates), MrBayes 3.2.7a
awk '/^>/{if(n)seq[n]=s; n++; h=substr($0,2); sub(/\/.*/,"",h); name[n]=h; s=""; next}{s=s $0}
     END{seq[n]=s; printf "#NEXUS\nBEGIN DATA;\n  DIMENSIONS NTAX=%d NCHAR=%d;\n  FORMAT DATATYPE=PROTEIN MISSING=? GAP=-;\nMATRIX\n", n, length(seq[1]);
         for(i=1;i<=n;i++) printf "%-30s %s\n", name[i], seq[i]; print ";\nEND;"}' \
     opsin_newA4_rev2.fa > oopsin_newA4_rev2.nex
cat > mb_opsin_rev2.nex << 'EOF'
#NEXUS
begin mrbayes;
   set autoclose=yes nowarn=yes;
   execute oopsin_newA4_rev2.nex;
   prset aamodelpr=fixed(lg);
   lset rates=gamma ngammacat=6;
   prset statefreqpr=fixed(empirical);
   mcmcp nruns=2 nchains=4 samplefreq=500 printfreq=500 diagnfreq=5000;
   mcmcp burninfrac=0.25 relburnin=yes filename=opsin_newA4_rev2_mb;
   mcmc ngen=15000000;
   sumt relburnin=yes burninfrac=0.25;
   sump relburnin=yes burninfrac=0.25;
quit;
end;
EOF
# (run used: 5e6 generations, then continued from checkpoint with append=yes to 1.5e7)
mb mb_opsin_rev2.nex > mb_rev2_output.log 2>&1

# 5) AU test: ML tree vs best tree with O. fusiformis CAH1779377 (Annelid_ropsins_16)
#    constrained as sister to the r-opsin1 + r-opsin3 clade
names=$(grep -oE "[A-Za-z]+_[A-Za-z]+_[0-9]+(/1-[0-9]+)?" opsin_newA4_rev2.treefile)
pick(){ for p in "$@"; do echo "$names" | grep -E "^${p}(/|$)"; done | paste -sd, ; }
CL=$(pick Annelid_ropsins_1 Annelid_ropsins_8 Mfu_opsins_3 Mfu_opsins_4 Annelid_ropsins_2 Annelid_ropsins_3 \
          Annelid_ropsins_6 Annelid_ropsins_7 Annelid_ropsins_9 Annelid_ropsins_11 Annelid_ropsins_13 \
          Annelid_ropsins_14 Annelid_ropsins_15)
OWN=$(pick Annelid_ropsins_16)
MOL=$(pick BrachiopodaMolluscs_ropsins_16 BrachiopodaMolluscs_ropsins_17 BrachiopodaMolluscs_ropsins_18)
echo "(($CL),$OWN,$MOL);" > constr_preDup.nwk
iqtree2 -s opsin_newA4_rev2.fa -m LG+F+R8 -g constr_preDup.nwk -pre constr_preDup -nt 2
cat opsin_newA4_rev2.treefile constr_preDup.treefile > au_trees.nwk
iqtree2 -s opsin_newA4_rev2.fa -m LG+F+R8 -z au_trees.nwk -zb 10000 -au -n 0 -pre au_test -nt 2 -redo
grep -A8 "USER TREES" au_test.iqtree

# 6) Residues distinguishing r-opsin1 from r-opsin3 (Sup. Fig. 3)
python3 ../Diagnostic_residues.py opsin_newA4_rev2.fa opsins_listA_Tree_Labels_rev.csv \
        SupFig3_diagnostic ../listA_filtered_add.fa
