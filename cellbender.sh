cellbender remove-background \
--input /d2/studies/scanPy/VM_NDB_Stress/D1_4_NDB_ctrl/outs/raw_feature_bc_matrix.h5 \
--output /d2/studies/scanPy/VM_NDB_Stress/CellBender/VM_D1_NDB_Ctrl_CellBender_Expect13k_Total50k/VM_NDB_Ctrl_CellBender_Matrix.h5 \
--cuda --expected-cells 13000 \
--total-droplets-included 50000 \
--fpr .01 \
--epochs 150
