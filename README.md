# BSC-MDD
Key words: Major Depressive Disorder; Diagnosis and Classification; Gender Differences.

Collection of data sources and scripts for "Shifts in the Brain Sex Continuum in Major Depressive Disorder: Evidence for a Persistent Neurobiological Marker"

## Scripts

Scripts were developed using R and Matlab.

### Step 1: Multi-site Harmonization
Harmonization of individual rsFCs in the PKU and XYH cohort with reference to the UKBiobank dataset using ComBat function in “sva” package in R. (1 Multi-site harmonization.R)

### Step 2: BSC Calculation
Computing the Brain Sex Continuum (BSC) score and split BSC scores within and between 6 networks. (2_BSC_calculation.m)

### Step 3: Classification Accuracy and Group Comparison
Calculation of BSC classification accuracy by sex, and comparison of BSC scores between depressed patients and healthy controls. (3 Classification accuracy and group comparison.R)

### Step 4: Multivariate Classifier
Classifier used to distinguish depressed patients from healthy controls. (4 Multivariate_classifier.R)

### Step 5: Associations between the BSC score and the Symptom Severity
Evaluation of the relationship between BSC scores and symptom severity, as measured by HAMD at baseline and at the 8-week follow-up. (5 Associations between the BSC score and the Symptom Severity.R)

## Data sources

The XYH data is part of the ZIB-MDD data, which are available by application to the ZIB data management committee after evaluation according to an established procedure (https://zib.fudan.edu.cn/). The PKU data are available upon request to the project PI Prof. Tianmei Si.

## Contact
```
Qiang Luo, PhD, Visiting Fellow at the Clare Hall, Cambridge
Associate Principal Investigator
Institute of Science and Technology for Brain-Inspired Intelligence (ISTBI)
Fudan University
Email: mrqiangluo@gmail.com;      qluo@fudan.edu.cn
Office Phone: 86-21-65648454
https://sites.google.com/site/qluochina/
http://homepage.fudan.edu.cn/qiangluo/
```

## Team Members:
```
Yinhan Chen, BSc, yinhanchen23@m.fudan.edu.cn
Linlin Zhu, Ph.D., zhulinlin@bjmu.edu.cn
Yi Zhang, Ph.D., yz669@cam.ac.uk
Jitao Li, Ph.D., ljt_102124@163.com
Christelle Langley, Ph.D., cl798@medschl.cam.ac.uk
Barbara J. Sahakian, Ph.D., bjs1001@cam.ac.uk
Jianfeng Feng, Ph.D., jffeng@fudan.edu.cn
Xiang Wang, Ph.D., wang0916xia@gmail.com
Tianmei Si, Ph.D., si.tian-mei@163.com
```

The current version was relesed on 22 Feb 2026.

Citation of this package: Chen et al. Shifts in the Brain Sex Continuum in Major Depressive Disorder: Evidence for a Persistent Neurobiological Marker.
