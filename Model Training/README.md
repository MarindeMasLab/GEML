# **Patient Classification Pipeline Using XGBoost and SHAP-Based Feature Selection**

## 📌 Overview  
This repository contains a multi-step pipeline for patient classification based on metabolomic data. The process includes data preprocessing, feature selection using SHAP values, and predictive modeling using XGBoost.

## 🛠 Workflow  

### **1️⃣ Data Preprocessing (`1_patient_classification_pipeline.py`)**  
- Processes raw patient data and outputs `.csv` files, each corresponding to an individual patient.  
- These preprocessed files serve as input for the modeling stage.  

### **2️⃣ Modeling and SHAP Analysis (`2_modeling_100_rounds_shap_analysis.ipynb`)**  
- Loads the preprocessed patient data.  
- Assigns target labels based on `Patient_Trauma_Groups.csv`.  
- Trains an XGBoost model on the complete dataset.  
- Uses SHAP (SHapley Additive exPlanations) values to determine feature importance.  
- Stores the most significant features in `Ranking_global_sharp.csv`.  

### **3️⃣ Feature-Selected Modeling (`3_top_features_modeling.py`)**  
- Models the dataset using only the top 50, 30, 20, and 10 most important features identified by SHAP.  
- Runs multiple iterations of training and evaluation to assess performance across different feature subsets.  
- Generates classification reports, confusion matrices, and feature importance visualizations.  

## 📂 Key Files  
| File | Description |
|------|------------|
| `1_patient_classification_pipeline.py` | Preprocesses raw patient data into structured `.csv` files. |
| `2_modeling_100_rounds_shap_analysis.ipynb` | Loads preprocessed data, assigns target labels, performs XGBoost modeling, and extracts feature importance using SHAP values. |
| `3_top_features_modeling.py` | Conducts modeling using only the most important features identified in SHAP rankings. |

## 🔧 Key Functions  
- **`extract_patient_index(file_path)`** → Extracts patient indices from file names.  
- **`load_preprocessed_data(preprocessed_path, num_patients)`** → Loads preprocessed patient data from `.csv` files.  
- **`combine_patient_data(patients, test_indices, target)`** → Merges individual patient datasets into a structured DataFrame with target labels.  
- **`conf_matrix(clf, X_train, X_test, y_train, y_test)`** → Generates confusion matrices for model evaluation.  

## 📊 Why Use SHAP-Based Feature Selection?  
SHAP values offer a robust method for interpreting machine learning models, ensuring that only the most informative features are used for classification. This enhances both model performance and interpretability.  
