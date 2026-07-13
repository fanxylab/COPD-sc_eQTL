import pandas as pd
import numpy as np
import torch
import tensorqtl
from tensorqtl import genotypeio, pgen, cis, trans, post
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"torch: {torch.__version__} (CUDA {torch.version.cuda}), device: {device}")
print(f"pandas: {pd.__version__}")

# Define 34 cell type for sc-eQTL mapping
cell_types = ['gCap','NKT_cell','Memory_CD4_T_cell','Naive_CD4_T_cell','Non-classical_monocytes','AT2a','CD8T_cell','NK_cell','cDC2','AT1',
              'Classical_monocytes','Aerocyte','Adventitial_fibroblast','Arterial','Alveolar_macrophage','Culb_1','Venous','Interstitial_macrophages',
              'Lymphatic','AT2b','Culb_2','Mast_cell','Goblet','Ciliated','Transitional_AT2','XCL1+_T_cell','B_cell','Neutrophils','Treg_T_cell', 
              'Plasma_cell','Alveolar_fibroblast','Proliferating_T_cells','Fibroblast','Basal'
             ]

plink_prefix_path = '/datg/xuxiaopeng/sc_eQTL/01_tensorqtl/genotype/COPD_maf0.1_R0.8_filtered'

# PLINK reader for genotypes
pgr = pgen.PgenReader(plink_prefix_path)
genotype_df = pgr.load_genotypes()
variant_df = pgr.variant_df
snp_df = pgr.pvar_df

for cell_type in cell_types:
    # define paths to data
    expression_bed = f'/datg/xuxiaopeng/sc_eQTL/01_tensorqtl/phenotype/{cell_type}.bed'
    covariates_file = f'/datg/xuxiaopeng/sc_eQTL/01_tensorqtl/covariates/{cell_type}.covariates.txt'
    
    prefix = f'{cell_type}_nominal'
    
    # load phenotypes and covariates
    phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
    covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T
    interaction_df = covariates_df[["Disease"]].copy()
    
    ######################################################### step1: 数据预处理 ###################################################
    base_cols = ['nCells', 'Sex', 'Age']
    pc_cols = [f'PC{i}' for i in range(1, 36)]
    gpc_cols = [f'gPC{i}' for i in range(1, 4)]
    
    df = covariates_df[base_cols + pc_cols + gpc_cols].copy()
    
    sex_mapping = {'Female': 0, 'Male': 1, np.nan: np.nan}
    disease_mapping = {'Control': 0, 'Case': 1, np.nan: np.nan}
    
    df['Sex'] = df['Sex'].map(sex_mapping).astype('Int64')
    interaction_df['Disease'] = interaction_df['Disease'].map(disease_mapping).astype('Int64')
    
    df.fillna(-1, inplace=True)
    interaction_df.fillna(-1, inplace=True)
    
    df = df.apply(pd.to_numeric, errors='coerce').astype('float64')
    interaction_df = interaction_df.apply(pd.to_numeric, errors='coerce').astype('float64')
    
    # genotype and phenotype match
    common_cols = sorted(set(phenotype_df.columns) & set(genotype_df.columns))
    phenotype_df_2 = phenotype_df.loc[:, common_cols]
    genotype_df_2 = genotype_df.loc[:, common_cols]
    covariates_df_2 = df.loc[common_cols, :]
    interaction_df_2 = interaction_df.loc[common_cols, :]
    
    ######################################################### step2: eQTL mapping ###################################################
    # map all cis-associations (results for each chromosome are written to file)
    cis.map_nominal(genotype_df_2, variant_df,
                    phenotype_df_2, phenotype_pos_df,
                    prefix, 
                    covariates_df = covariates_df_2,
                    interaction_df = interaction_df_2, maf_threshold_interaction=0.05,
                    run_eigenmt = True, write_top = True, write_stats = True,
                    output_dir='/datg/xuxiaopeng/sc_eQTL/01_tensorqtl/interaction_result'
                   )
    
  
