import pandas as pd
import numpy as np
import torch
import tensorqtl
from tensorqtl import genotypeio, pgen, cis, trans, post
from tqdm import tqdm
import time
import logging

# Initialize logging
logging.basicConfig(
  filename='eQTL_pipeline.log',
  level=logging.INFO,
  format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger()

# Device configuration
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"torch: {torch.__version__} (CUDA {torch.version.cuda}), device: {device}")
print(f"pandas: {pd.__version__}")

# Define cell types
cell_types = ['gCap','NKT_cell','Memory_CD4_T_cell','Naive_CD4_T_cell','Non-classical_monocytes','AT2a','CD8T_cell','NK_cell','cDC2','AT1',
              'Classical_monocytes','Aerocyte','Adventitial_fibroblast','Arterial','Alveolar_macrophage','Culb_1','Venous','Interstitial_macrophages',
              'Lymphatic','AT2b','Culb_2','Mast_cell','Goblet','Ciliated','Transitional_AT2','XCL1+_T_cell','B_cell','Neutrophils','Treg_T_cell', 
              'Plasma_cell','Alveolar_fibroblast','Proliferating_T_cells','Fibroblast','Basal'
]

# Load genotype data
plink_prefix_path = '/datg/xuxiaopeng/sc_eQTL/01_tensorqtl/genotype/COPD_maf0.1_R0.8'
pgr = pgen.PgenReader(plink_prefix_path)
genotype_df = pgr.load_genotypes()
variant_df = pgr.variant_df
snp_df = pgr.pvar_df

# Global timer
global_start = time.time()
logger.info(f"Pipeline started at {time.strftime('%Y-%m-%d %H:%M:%S')}")

# Print startup banner
print(f"\n{'='*60}")
print(f"🏁 Starting eQTL analysis pipeline | Total cell types: {len(cell_types)}")
print(f"⏰ Start time: {time.strftime('%Y-%m-%d %H:%M:%S')}")
print(f"{'='*60}\n")

# Main processing loop
for cell_idx, cell_type in enumerate(tqdm(cell_types, desc="Overall progress", unit="cell"), 1):
  cell_start = time.time()
logger.info(f"Processing cell type {cell_idx}/{len(cell_types)}: {cell_type}")

try:
  # ================== Data Loading Phase ==================
print(f"\n🌀 [{cell_idx:02d}/{len(cell_types)}] Processing {cell_type}")
print(f"🕒 Current time: {time.strftime('%H:%M:%S')}")
print("▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬▬")

# Load phenotype and covariates
expression_bed = f'/datg/xuxiaopeng/sc_eQTL/01_tensorqtl/phenotype/{cell_type}.bed'
covariates_file = f'/datg/xuxiaopeng/sc_eQTL/01_tensorqtl/covariates/{cell_type}.covariates.txt'
prefix = f'{cell_type}_nominal'

phenotype_df, phenotype_pos_df = tensorqtl.read_phenotype_bed(expression_bed)
covariates_df = pd.read_csv(covariates_file, sep='\t', index_col=0).T

# ================== Data Preprocessing ==================
base_cols = ['nCells', 'Sex', 'Age', 'Disease']
pc_cols = [f'PC{i}' for i in range(1, 36)]
gpc_cols = [f'gPC{i}' for i in range(1, 4)]

df = covariates_df[base_cols + pc_cols + gpc_cols].copy()
sex_mapping = {'Female': 1, 'Male': 2, np.nan: np.nan}
disease_mapping = {'Control': 2, 'Case': 1, np.nan: np.nan}

df['Sex'] = df['Sex'].map(sex_mapping).astype('Int64')
df['Disease'] = df['Disease'].map(disease_mapping).astype('Int64')
df.fillna(-1, inplace=True)
df = df.apply(pd.to_numeric, errors='coerce').astype('float64')

# Data alignment
common_cols = sorted(set(phenotype_df.columns) & set(genotype_df.columns))
phenotype_df_2 = phenotype_df.loc[:, common_cols]
genotype_df_2 = genotype_df.loc[:, common_cols]
covariates_df_2 = df.loc[common_cols, :]

print(f"|-- Common samples: {len(common_cols)}")

# ================== eQTL Mapping ==================
cis.map_nominal(
  genotype_df_2, variant_df,
  phenotype_df_2, phenotype_pos_df,
  prefix, 
  covariates_df=covariates_df_2,
  output_dir='/datg/xuxiaopeng/sc_eQTL/01_tensorqtl/nominal_result'
)

# ================== Result Merging ==================
combined_df = pd.DataFrame()
print(f"|-- Merging results from 22 chromosomes...")

# Chromosome processing progress
chr_bar = tqdm(range(1,23), desc="Chromosomes processing", leave=False)
for chr in chr_bar:
  chr_bar.set_postfix({'Current': f'chr{chr}'})
pairs_df = pd.read_parquet(
  f'/datg/xuxiaopeng/sc_eQTL/01_tensorqtl/nominal_result/{prefix}.cis_qtl_pairs.{chr}.parquet'
)
output_df = pairs_df.merge(
  snp_df.set_index('id')[['chrom', 'pos', 'ref', 'alt']], 
  left_on='variant_id', 
  right_index=True, 
  how='left'
)
output_df['maf'] = np.minimum(output_df['af'], 1 - output_df['af'])
output_df['n'] = len(common_cols)
combined_df = pd.concat([combined_df, output_df], axis=0)
chr_bar.close()

# ================== Result Saving ==================
output_path = f"/datg/xuxiaopeng/sc_eQTL/01_tensorqtl/nominal_result/{cell_type}.nominal.all.txt"
combined_df.to_csv(output_path, sep='\t', index=False)
logger.info(f"Saved {combined_df.shape[0]} records to {output_path}")

# ================== Completion Report ==================
cell_elapsed = time.time() - cell_start
total_elapsed = time.time() - global_start
print(f"\n✅ Completed [{cell_idx:02d}/{len(cell_types)}] {cell_type}")
print(f"⏱ Current duration: {cell_elapsed//60:.0f}m {cell_elapsed%60:.2f}s")
print(f"⌛️ Total elapsed: {total_elapsed//3600:.0f}h {total_elapsed%3600//60:.0f}m")
print(f"📊 Records generated: {combined_df.shape[0]:,} associations")
print(f"▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔▔\n")

except Exception as e:
  error_msg = f"❌ [{cell_type}] Processing failed: {str(e)}"
print(f"\n{error_msg}")
logger.error(error_msg, exc_info=True)
continue

# ================== Final Report ==================
total_time = time.time() - global_start
print(f"\n{'='*60}")
print(f"🎉 All tasks completed! Total duration: {total_time//3600:.0f}h {total_time%3600//60:.0f}m")
print(f"📂 Output directory: /datg/xuxiaopeng/sc_eQTL/01_tensorqtl/nominal_result")
print(f"📈 Statistics: {len(cell_types)} cell types × 22 chromosomes")
print(f"{'='*60}")
