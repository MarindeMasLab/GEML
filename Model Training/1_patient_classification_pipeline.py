# Patient Metabolic Data Analysis and Classification Pipeline

"""
Title: Integrating Machine Learning with Metabolic Models for Precision Trauma Care:
Personalized ENDOTYPE Stratification and Metabolic Target Identification

Authors:
- Igor Marin de Mas (Copenhagen University Hospital, Rigshospitalet)
- Lincoln Moura (Universidade Federal do Ceará)
- Fernando Luiz Marcelo Antunes (Universidade Federal do Ceará)
- Josep Maria Guerrero (Aalborg University)
- Pär Ingemar Johansson (Copenhagen University Hospital, Rigshospitalet)

Description:
This script implements a complete pipeline for the analysis and classification of metabolic data from patients.
The main steps include:

1. Loading preprocessed data: Patient data stored in CSV files are loaded and organized.
2. Data combination: Data from different patients are combined into a single DataFrame, including target classification groups.
3. Preprocessing: Checks for missing values and adjusts the data to ensure integrity.
4. Model training: Uses XGBoost to perform multiclass classification of metabolic groups.
5. Model evaluation: Performs cross-validation using Stratified K-Fold and presents metrics such as accuracy and confusion matrix.

This script is ideal for exploring metabolic patterns among patients and evaluating the performance of classification models.
"""

# Standard libraries
import os
import glob

# Third-party libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
from sklearn.metrics import classification_report
from xgboost import XGBClassifier
import scikitplot as skplt
import warnings

# Suppress warnings
warnings.filterwarnings("ignore")

# -----------------------------
# Utility Functions
# -----------------------------

def extract_patient_index(file_path):
    """
    Extract the patient index from the file name using regex.

    Args:
        file_path (str): Path to the file.

    Returns:
        int: Patient index extracted from the file name.

    Raises:
        ValueError: If the index cannot be extracted from the file name.
    """
    import re
    match = re.search(r'(\d+)', os.path.basename(file_path))
    if match:
        return int(match.group(1)) - 1
    else:
        raise ValueError(f"Cannot extract index from file name: {file_path}")


def load_preprocessed_data(preprocessed_path, num_patients):
    """
    Load preprocessed patient data from CSV files.

    Args:
        preprocessed_path (str): Path to the directory with preprocessed data.
        num_patients (int): Number of patients.

    Returns:
        dict: Dictionary of patient DataFrames indexed by patient number.
        list: List of indices for successfully loaded patients.
    """
    patients = ["patient_" + str(x) for x in range(num_patients)]
    test_indices = []

    for file_path in glob.glob(preprocessed_path + '*.csv*'):
        try:
            index = extract_patient_index(file_path)
            print(f"Index = {index}")
            test_indices.append(index)
            patients[index] = pd.read_csv(file_path, index_col=0)
        except ValueError as e:
            print(f"[ERROR] {e}")
            continue
        except Exception as e:
            print(f"[ERROR] Could not load file: {file_path}. Error: {e}")
            continue

    return patients, test_indices


def combine_patient_data(patients, test_indices, target):
    """
    Combine individual patient data into a single DataFrame with target labels.

    Args:
        patients (dict): Dictionary of patient DataFrames.
        test_indices (list): List of indices for patients.
        target (list): Target labels for patients.

    Returns:
        pd.DataFrame: Combined DataFrame with all patient data and target labels.
    """
    dataframes = []

    for i in test_indices:
        temp = patients[i].T
        temp["target"] = target[i]
        dataframes.append(temp)

    df = pd.concat(dataframes, axis=0).reset_index(drop=True)

    # Check for missing values
    if df.isna().sum().sum() > 0:
        print("[WARNING] Missing values detected. Filling with column means.")
        df.fillna(df.mean(), inplace=True)

    return df


def plot_explainability_histogram(explainability):
    """
    Plot a histogram of explained variance from PCA data.

    Args:
        explainability (pd.DataFrame): DataFrame containing explained variance.
    """
    ax = explainability.hist(figsize=(10, 5))
    plt.title("")
    plt.xlabel("Explained Variance Using PCA with 600 Components (%)", fontsize=12)
    plt.ylabel("Number of patients", fontsize=12)
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.gca().set_xticklabels([f'{x:.2%}' for x in plt.gca().get_xticks()])
    plt.show()


def conf_matrix(clf, X_train, X_test, y_train, y_test):
    """
    Plot confusion matrices for training and testing predictions.

    Args:
        clf: Trained classifier.
        X_train: Training data.
        X_test: Testing data.
        y_train: Training labels.
        y_test: Testing labels.
    """
    Y_train_pred = clf.predict(X_train)
    Y_test_pred = clf.predict(X_test)

    fig = plt.figure(figsize=(15, 6))
    ax1 = fig.add_subplot(121)
    skplt.metrics.plot_confusion_matrix(Y_train_pred, y_train, normalize=False, title="Confusion Matrix", cmap="Oranges", ax=ax1)
    ax2 = fig.add_subplot(122)
    skplt.metrics.plot_confusion_matrix(Y_test_pred, y_test, normalize=False, title="Confusion Matrix", cmap="Purples", ax=ax2)
    plt.show()


def train_and_evaluate_model(X, y, n_splits=10, n_rounds=5):
    """
    Train and evaluate the model using cross-validation and stratified k-fold.

    Args:
        X (pd.DataFrame): Feature data.
        y (pd.Series): Target labels.
        n_splits (int): Number of splits for k-fold.
        n_rounds (int): Number of rounds for cross-validation.

    Returns:
        pd.DataFrame: DataFrame containing evaluation metrics for each round.
    """
    results = []

    for t in range(n_rounds):
        print(f"[INFO] Round {t}")
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.25, stratify=y)
        X_train.reset_index(drop=True, inplace=True)
        y_train.reset_index(drop=True, inplace=True)

        kfold = StratifiedKFold(n_splits=n_splits, shuffle=True)
        model = XGBClassifier(objective='multi:softmax')
        scores = cross_val_score(model, X_train, y_train, cv=kfold, scoring="accuracy")
        results.append({
            "round": t,
            "mean_score": np.mean(scores),
            "std_score": np.std(scores)
        })

        # Final test with all data
        model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        print(classification_report(y_test, y_pred))
        conf_matrix(model, X_train, X_test, y_train, y_test)

    return pd.DataFrame(results)


# -----------------------------
# Main Script
# -----------------------------

def main():
    # Load target data
    target_df = pd.read_csv("Patient_Trauma_Groups.csv", delimiter=";")
    target = target_df["Metabo-group"].values

    # Load explainability data
    explainability = pd.read_csv("explainability.csv", delimiter=",", index_col=0).T
    plot_explainability_histogram(explainability)

    # Load preprocessed patient data
    preprocessed_path = "preprocess_PCA/"
    num_patients = 95
    patients, test_indices = load_preprocessed_data(preprocessed_path, num_patients)

    # Combine patient data into a single DataFrame
    df = combine_patient_data(patients, test_indices, target)

    # Prepare features and target
    X, y = df.drop('target', axis=1), df['target'] - 1

    # Train and evaluate model
    results_df = train_and_evaluate_model(X, y)
    print(results_df)


if __name__ == "__main__":
    main()
