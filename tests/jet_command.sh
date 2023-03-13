java -Xmx1000m \
    -cp .:./jet/extLibs/vecmath.jar \
    jet.JET \
    -c custom.conf \
    -i data/adrb2/adrb2.pdb \
    -o `pwd` \
    -p J \
    -r input \
    -f data/adrb2/adrb2_A.fasta \
    -d chain \
    -n 1 -l data/adrb2_jet.log > data/adrb2_jet_test.out