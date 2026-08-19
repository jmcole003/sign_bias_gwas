## Get lab measurements

import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

dataset = %env WORKSPACE_CDR
bucket = os.getenv("WORKSPACE_BUCKET")


# get lab measurements
labms = pd.read_gbq(f'''SELECT DISTINCT
  m.person_id,
  LOWER(m.measurement_source_value) AS measurement_source_value,
  m.value_as_number,
  COALESCE(u.concept_name, CAST(m.unit_source_value AS STRING)) AS unit_concept_name
FROM `{dataset}.measurement` AS m
LEFT JOIN `{dataset}.concept` AS u
  ON u.concept_id = m.unit_concept_id
WHERE LOWER(m.measurement_source_value) IN ('706-2','770-8','5905-5','6690-2','789-8','785-6')''')


# loinc names
lnames = {'706-2': 'basophil_percentage','770-8': 'neutrophil_percentage','5905-5': 
'monocyte_percentage', '6690-2': 'white_blood_cell_count', '789-8': 'red_blood_cell_count', 
'785-6': 'mean_corpuscular_hemoglobin'}

labms['lab'] = labms['measurement_source_value'].map(lnames)

# raw dists
raw = (labms
       .groupby('lab', dropna=False)['value_as_number']
       .agg(n='size',
            n_missing=lambda x: x.isna().sum(),
            min='min',
            median='median',
            max='max')
       .sort_values('max', ascending=False))


####################################################################################
## percentages

plabs = [
    'basophil_percentage',
    'neutrophil_percentage',
    'monocyte_percentage'
]

punits = [
    'percent',
    'percentage unit',
    'percent of white blood cells'
]


pct = labms[
    labms['lab'].isin(plabs) &
    labms['unit_concept_name'].isin(punits)
].copy()


# remove missing and values outside 0-100
pct = pct[
    pct['value_as_number'].between(0,100, inclusive='both')
].copy()


# median measurements
pctmed = (pct
          .groupby(['person_id','lab'])['value_as_number']
          .median()
          .reset_index())

pctw = (pctmed
        .pivot(index='person_id',
               columns='lab',
               values='value_as_number')
        .reset_index())


pctw = pctw.rename(columns={
    'basophil_percentage': 'basophil_percentage_median',
    'monocyte_percentage': 'monocyte_percentage_median',
    'neutrophil_percentage': 'neutrophil_percentage_median'
})


####################################################################################
## WBC

wbc = labms[labms['lab'] == 'white_blood_cell_count'].copy()

wbc['unit_concept_name'] = (
    wbc['unit_concept_name']
    .astype(str)
    .str.lower()
    .str.strip()
)

wbc_1 = {
    'thousand per microliter',
    'thousand per cubic millimeter',
    'x10(3)/mcl',
    'billion per liter'
}
wbc_1000 = {
    'cells per microliter',
    'cells/ul',
    'per microliter',
    'per cubic millimeter',
    '/mm3'
}


wbc['val'] = np.nan

wbc.loc[
    wbc['unit_concept_name'].isin(wbc_1),
    'val'
] = wbc['value_as_number']

wbc.loc[
    wbc['unit_concept_name'].isin(wbc_1000),
    'val'
] = wbc['value_as_number'] * 1e-3


# bounds
wbc = wbc[wbc['val'].between(0.1,500.0, inclusive='both')].copy()


wbcmed = (wbc
          .groupby('person_id')['val']
          .median()
          .reset_index(name='white_blood_cell_count_median'))


####################################################################################
## RBC

rbc = labms[labms['lab'] == 'red_blood_cell_count'].copy()

rbc['unit_concept_name'] = (
    rbc['unit_concept_name']
    .astype(str)
    .str.lower()
    .str.strip()
)


rbc_million = {
    'million per microliter',
    'million per cubic millimeter'
}
rbc_cells = {
    'cells per microliter',
    'cells/ul',
    'per microliter',
    'per cubic millimeter',
    '/mm3'
}

rbc_thousand = {
    'thousand per microliter',
    'thousand per cubic millimeter',
    'billion per liter'
}


rbc['val'] = np.nan

rbc.loc[
    rbc['unit_concept_name'].isin(rbc_million),
    'val'
] = rbc['value_as_number']

rbc.loc[
    rbc['unit_concept_name'].isin(rbc_cells),
    'val'
] = rbc['value_as_number'] * 1e-6

rbc.loc[
    rbc['unit_concept_name'].isin(rbc_thousand),
    'val'
] = rbc['value_as_number'] * 1e-3

rbc.loc[
    rbc['unit_concept_name'] == 'million per liter',
    'val'
] = rbc['value_as_number'] * 1e-6


# bounds
rbc = rbc[rbc['val'].between(0.5,12.0, inclusive='both')].copy()


rbcmed = (rbc
          .groupby('person_id')['val']
          .median()
          .reset_index(name='red_blood_cell_count_median'))


####################################################################################
## MCH

mchunits = {
    'picogram',
    'picogram per cell',
    'pg/cell'
}


mch = labms[
    (labms['lab'] == 'mean_corpuscular_hemoglobin') &
    (labms['unit_concept_name'].isin(mchunits))
].copy()


# only keep plausible values
mch = mch[
    mch['value_as_number'].between(5,100, inclusive='both')
].copy()


mchmed = (mch
          .groupby('person_id')['value_as_number']
          .median()
          .reset_index(name='mean_corpuscular_hemoglobin_median'))


####################################################################################
## combine

labsout = pctw.merge(
    mchmed,
    on='person_id',
    how='outer'
)

labsout = labsout.merge(
    rbcmed,
    on='person_id',
    how='outer'
)

labsout = labsout.merge(
    wbcmed,
    on='person_id',
    how='outer'
)


# column order
labsout = labsout[[
    'person_id',
    'basophil_percentage_median',
    'mean_corpuscular_hemoglobin_median',
    'monocyte_percentage_median',
    'neutrophil_percentage_median',
    'red_blood_cell_count_median',
    'white_blood_cell_count_median'
]]


# write out
labsout.to_csv("lab_measures_all_v2.tsv", sep="\t", index=False,na_rep="NA")

# copy to bucket
!gsutil -m -u $GOOGLE_PROJECT cp lab_measures_all_v2.tsv '{bucket}/aou_gwas/pheno/'