#!/bin/bash

allnames=""
steering=$(dirname $0)/material_budget_steering.txt
outdir=$1

# directory to store all the results
mkdir -p $outdir

for (( phi=-180; phi<=180; phi+=2 )); do
    cp $steering material_budget_tmp.txt
    sed -i "s/^phi [0-9]*/phi ${phi}/" material_budget_tmp.txt
    sed -i "s/^rootfile [^ ]*/rootfile material_budget_phi${phi}.root/" material_budget_tmp.txt
    $(which materialBudget) xml/OpenDataDetector.xml material_budget_tmp.txt
    mv material_budget_phi${phi}.root ${outdir}/.
    allnames+=" ${outdir}/material_budget_phi${phi}.root"
done

rm material_budget_tmp.txt

python3 $(dirname $steering)/material_budget_dd4hep_plots.py -i ${allnames} -o ${outdir}/material_budget_results.root

rm ${allnames}
