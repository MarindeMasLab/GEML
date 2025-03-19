########################
## 0. Load libraries ##
########################
import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import pymc as pm
import arviz as az
import matplotlib.pyplot as plt

####################################
## 1. Load data and pre-treatment ##
####################################

# Updated DataFrame with Age
data = pd.DataFrame({
    'Cluster': [3, 3, 3, 4, 3, 4, 3, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 1, 3, 4, 3, 3, 3, 4, 3, 3, 3, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3],
    'Sex': [0, 1, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0],
    'Age': [63, 51, 60, 65, 43, 51, 39, 65, 59, 47, 63, 22, 48, 45, 57, 71, 54, 52, 42, 59, 32, 34, 38, 65, 31, 30, 50, 65, 41, 48, 37, 45, 51, 60, 42, 61, 35, 63, 48, 55, 41],
    'Improvement': [0, 37.8921568627451, 0, 40.1764705882353, 0, 0, 72.5392156862745, 0, 0, 13.4117647058824, 76.0588235294118, 68.8235294117647, 66.8725490196078, 0, 0.245098039215685, 0, 65.7647058823529, 0, 7.70588235294117, 4.35294117647059, 0, 56.9411764705882, 0, 0, 70.2156862745098, 2.11764705882352, 43.2941176470588, 2.98039215686274, 18.4705882352941, 15.6764705882353, 6.82352941176471, 0, 0.117647058823522, 0.823529411764712, 0, 2.82352941176471, 0.117647058823522, 2.94117647058824, 0, 24.5882352941177, 0],
    'Patient_ID': [4, 5, 6, 10, 12, 13, 15, 17, 22, 24, 25, 26, 27, 28, 32, 33, 35, 44, 46, 47, 54, 55, 57, 60, 65, 66, 67, 68, 69, 72, 74, 75, 78, 81, 83, 84, 87, 88, 89, 91, 92]
})

# Data rounding
data['Improvement_rounded'] = data['Improvement'].round(2)

# Standardize Age
data['Age_scaled'] = StandardScaler().fit_transform(data[['Age']])

data['Age_centered'] = data['Age'] - data['Age'].mean()

kmeans = KMeans(n_clusters=2, random_state=0)
data['Age_group_cluster'] = kmeans.fit_predict(data[['Age']])

threshold = data['Age'].median()
data['Age_group'] = (data['Age'] > threshold).astype(int)


# Convert Cluster to categorical
data['Cluster'] = data['Cluster'].astype('category')
data['Cluster'] = data['Cluster'].replace({1: 2, 4: 3}) # Combine small clusters (Cluster 1 with 2, and Cluster 3 with 4)


# Preprocessing: Convert categorical variables to numeric
data['Cluster'] = data['Cluster'].astype('category').cat.codes
data['Sex'] = data['Sex'].astype('category').cat.codes

# Handle missing data by dropping rows with NaN values (if any)
data = data.dropna(subset=['Improvement', 'Cluster', 'Sex', 'Age_scaled', 'Age_group'])

########################
## 2. Bayesian Models ##
########################

# Define the Bayesian model
with pm.Model() as bayesian_model:
    # Priors
    cluster_effect = pm.Normal('Cluster_Effect', mu=0, sigma=10, shape=2)  # Assuming 2 clusters    
    sex_effect = pm.Normal('Sex_Effect', mu=0, sigma=10)
    age_effect = pm.Normal('Age_Effect', mu=0, sigma=10)    
    # Interaction terms
    age_cluster_interaction = pm.Normal('Age_Cluster_Interaction', mu=0, sigma=10, shape=2)  # Age x Cluster interaction
    age_sex_interaction = pm.Normal('Age_Sex_Interaction', mu=0, sigma=10)  # Age x Sex interaction
    cluster_sex_interaction = pm.Normal('Cluster_Sex_Interaction', mu=0, sigma=10, shape=2)  # Cluster x Sex interaction    
    sigma = pm.HalfNormal('sigma', sigma=10)    
    # Model Equation including all interactions    
    # Option 1: Include age as a continuous variable (use Age_scaled or Age_centered)  
    improvement = pm.Normal(
        'Improvement',
        mu=cluster_effect[data['Cluster'].values] + 
           sex_effect * data['Sex'].values + 
           age_effect * data['Age_centered'].values + 
           age_cluster_interaction[data['Cluster'].values] * data['Age_centered'].values +  # Age x Cluster interaction
           age_sex_interaction * data['Sex'].values * data['Age_centered'].values +  # Age x Sex interaction
           cluster_sex_interaction[data['Cluster'].values] * data['Sex'].values,  # Cluster x Sex interaction
        sigma=sigma,
        observed=data['Improvement']
    )
    # Sampling
    trace = pm.sample(2000, return_inferencedata=True, tune=1000)

# Posterior predictive sampling
with bayesian_model:
    ppc = pm.sample_posterior_predictive(trace)

# Extract the posterior predictive samples for the "Improvement" variable
improvement_ppc = ppc.posterior_predictive["Improvement"].values
print(improvement_ppc.shape)  # This should show (chains, draws, number of observations)

# Optionally, convert the posterior predictive samples to a DataFrame
ppc_df = pd.DataFrame(improvement_ppc.reshape(-1, improvement_ppc.shape[-1]).T)
print(ppc_df.head())  # Show the first few rows of the DataFrame

# Convert trace to dictionary format compatible with ArviZ
posterior_dict = {var: trace.posterior[var].values for var in trace.posterior.keys()}

# Create InferenceData object using the extracted posterior dictionary
az_combined = az.from_dict(
    posterior=posterior_dict,
    posterior_predictive={"Improvement": improvement_ppc},
    observed_data={"Improvement": data["Improvement"].values}
)


##################################
## 3. Compute WAIC, LOOIC, RMSE ##
##################################

# Compute WAIC (Widely Applicable Information Criterion)
# WAIC can be computed from the posterior predictive check (ppc)
# The ppc already contains the predictions we need.

# Compute the log-likelihood for each posterior sample
sigma_samples = trace.posterior['sigma'].values # Extract posterior samples for sigma
log_likelihood_samples = [] # Initialize a list to store log-likelihood samples

# Loop over each chain
for i in range(ppc.posterior_predictive["Improvement"].shape[0]):  # Loop over chains
    improvement_sample = ppc.posterior_predictive["Improvement"][i]  # Shape: (2000, 41) for the i-th chain    
    # Loop over each draw (i.e., each set of predictions)
    for j in range(improvement_sample.shape[0]):  # Loop over draws (2000)
        # improvement_sample[j, :] corresponds to a set of predictions for all 41 observations        
        # Ensure the correct sigma value is used for the current draw (from sigma_samples)
        sigma = sigma_samples[i, j]  # This corresponds to the sigma for the i-th chain, j-th draw        
        # Compute the log-likelihood for each data point
        log_likelihood = -0.5 * np.sum(np.log(2 * np.pi * sigma ** 2) + 
                                       ((data['Improvement'].values - improvement_sample[j, :]) ** 2) / (sigma ** 2))
        log_likelihood_samples.append(log_likelihood)
        
log_likelihood_samples = np.array(log_likelihood_samples) # Convert the list of log-likelihoods to a numpy array

# WAIC computation
waic = -2 * np.sum(np.log(np.mean(np.exp(log_likelihood_samples), axis=0)))
print(f"WAIC: {waic}")

# Compute LOOIC (Leave-One-Out Information Criterion)
# LOOIC is computed using the log-likelihood for each data point across posterior samples.
looic = -2 * np.sum(np.log(np.mean(np.exp(log_likelihood_samples), axis=0)))
print(f"LOOIC: {looic}")

# Compute RMSE (Root Mean Square Error)
# To compute RMSE, we will use the observed data and the posterior predictive samples
# Calculate RMSE for each posterior sample
rmse_samples = np.sqrt(np.mean((ppc.posterior_predictive["Improvement"] - data['Improvement'].values) ** 2, axis=1))
rmse_mean = np.mean(rmse_samples)
print(f"RMSE: {rmse_mean}")


#################
## 4. Plotting ##
#################
## Posteriors
# Posterior Predictive Check (PPC)
az.plot_ppc(az_combined)
plt.show()

# Trace plot for parameters
'''
This generates trace plots for each parameter in your model.
It shows the posterior distributions (histograms/density plots) on one side and the MCMC sampling traces on the other.
You can interpret the posterior distributions to see the range of likely parameter values and their uncertainty.
'''
az.plot_trace(trace)
plt.show()

# Rank Plots
'''
Rank plots assess the convergence and mixing of the chains by showing rank-order distributions.
Useful for debugging MCMC convergence issues.
'''
az.plot_rank(trace)
plt.show()

# Posterior Summary Plots
'''
This shows only the posterior distributions (no traces), including summary statistics like the mean, median, and credible intervals for each parameter.
It’s ideal for focusing on the inference and credible intervals.
'''
az.plot_posterior(trace)
plt.show()

# Pair Plots
'''
Visualizes pairwise relationships between parameters to understand how they are correlated.
'''
az.plot_pair(trace, var_names=["Cluster_Effect", "Sex_Effect", "Age_Effect"], kind="kde", marginals=True)
plt.show()

# Density Overlay
'''
To compare the posterior densities of specific parameters:
Displays credible intervals as horizontal bars for parameters, making it easy to compare them.
'''
az.plot_forest(trace, var_names=["Age_Effect", "Sex_Effect", "Cluster_Effect"], combined=True)
plt.show()


#####################
## 5. Save results ##
#####################
# Summary results
print(az.summary(trace)) # Optional: Summary of the model parameters
summary = az.summary(trace) # Generate the summary (assuming `trace` exists)
summary_df = pd.DataFrame(summary) # Convert the summary to a DataFrame
summary_df.to_csv('bayesian_model_summary.csv', index=True) # Save the DataFrame to a CSV file
print("Summary saved to 'bayesian_model_summary.csv'.")

# Save the WAIC, LOOIC, and RMSE to a CSV file
results = {
    'WAIC': waic,
    'LOOIC': looic,
    'RMSE': rmse_mean
}
# Convert results to a DataFrame and save
results_df = pd.DataFrame([results])
results_df.to_csv('model_performance_metrics.csv', index=False)

print("Performance metrics saved to 'model_performance_metrics.csv'.")


