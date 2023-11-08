import pandas as pd
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
import numpy as np

def load_data(file_path):
    # Load data from the Excel file
    data = pd.read_excel(file_path)
    return data[["Identity(%)", "Match_Length "]].values

def plot_initial_data(X):
    # Scatter plot of the data
    plt.figure(figsize=(12, 6))
    plt.scatter(X[:, 0], X[:, 1], s=30, cmap='viridis')
    plt.xlabel('Identity(%)')
    plt.ylabel('Match_Length')
    plt.title('Scatter Plot of Identity(%) vs Match_Length')
    plt.grid(True)
    plt.show()

def plot_elbow_method(X):
    # Determine the optimal number of clusters using the "elbow" method
    n_components = np.arange(1, 21)
    models = [GaussianMixture(n, random_state=0).fit(X) for n in n_components]
    plt.figure(figsize=(12, 6))
    plt.plot(n_components, [m.bic(X) for m in models], label='BIC', marker='o')
    plt.plot(n_components, [m.aic(X) for m in models], label='AIC', marker='o')
    plt.legend(loc='best')
    plt.xlabel('n_components')
    plt.title('BIC and AIC for different number of clusters')
    plt.grid(True)
    plt.show()

def apply_gmm(X, n_clusters):
    # Fitting the GMM with the specified number of clusters
    gmm = GaussianMixture(n_components=n_clusters, random_state=0).fit(X)
    labels = gmm.predict(X)
    
    # Plotting the results
    plt.figure(figsize=(12, 6))
    plt.scatter(X[:, 0], X[:, 1], c=labels, s=30, cmap='viridis')
    plt.xlabel('Identity(%)')
    plt.ylabel('Match_Length')
    plt.title(f'Clusters from GMM (n_clusters={n_clusters})')
    plt.grid(True)
    plt.colorbar()
    plt.show()

def main():
    # Load data
    file_path = "path_to_your_excel_file.xlsx"
    X = load_data(file_path)
    
    # Plot initial data
    plot_initial_data(X)
    
    # Plot elbow method
    plot_elbow_method(X)
    
    # Input desired number of clusters
    n_clusters = int(input("Enter the desired number of clusters based on the elbow method: "))
    
    # Apply GMM and plot results
    apply_gmm(X, n_clusters)

if __name__ == "__main__":
    main()
