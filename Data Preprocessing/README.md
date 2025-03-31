GEMs folder: Metabolic models generated in the pipeline described in Silva-lace et al 2024 in COBRA Matlab format:
	- iEC3006 GEM: Reference model
	- Patient-specific models: The models are named following the format modelPatient_XX_YYY.mat where:
		- XX: Patient number ID
		- YYY: Boundaries relative to Healthy individual group Maximum (YYY = Max),  Mean (YYY = Meam) or  Manimum (YYY = Min)

matlab_data_preprocessing.ipynb: Code is designed to preprocess and analyze datasets stored in MATLAB .mat files, with the goal of cleaning and transforming the data for further analysis. It performs several essential steps, including data validation, handling missing values, removing outliers, numerical stabilization, and dimensionality reduction using Principal Component Analysis (PCA). The preprocessed data is then saved in CSV format for easy access and further use.

References:
1. Silva-Lance F, Montejano-Montelongo I, Bautista E, Nielsen LK, Johansson PI, Marin de Mas I. Integrating Genome-Scale Metabolic Models with Patient Plasma Metabolome to Study Endothelial Metabolism In Situ. Int J Mol Sci. 2024 May 15;25(10):5406. doi: 10.3390/ijms25105406. PMID: 38791446; PMCID: PMC11121795.