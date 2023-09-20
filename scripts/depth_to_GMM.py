#Under Construction

read_depth_data = pd.read_csv("/g/data/xl04/ka6418/chromosome_graph/rTilRug_HiC_pctg_PacBio_sorted.bam.binned.fixed.depth.csv")
def plot_gaussian_distributions(data, n_clusters, contig_id):
    """
    Plots Gaussian distributions for read depth data based on a specified number of clusters.
    
    Parameters:
    - data: DataFrame with read depth data.
    - n_clusters: Integer specifying the number of Gaussian components or clusters.
    - contig_id: ID of the contig for visualization.
    """
    # Filter data for the specified contig
    contig_depth_data = data[data['Chromosome'] == contig_id]
    read_depth_values = contig_depth_data['AverageDepth'].values.reshape(-1, 1)

    # Fitting a GMM with specified number of components
    gmm = GaussianMixture(n_components=n_clusters, random_state=0)
    gmm.fit(read_depth_values)
    labels = gmm.predict(read_depth_values)

    # Plotting
    plt.figure(figsize=(14, 7))
    
    colors = plt.cm.viridis(np.linspace(0, 1, n_clusters))
    
    for i in range(n_clusters):
        cluster_data = read_depth_values[labels == i]
        
        # Calculating summary statistics for the cluster
        cluster_mean = cluster_data.mean()
        cluster_std = cluster_data.std()

        # Define the domain for plotting
        x = np.linspace(read_depth_values.min(), read_depth_values.max(), 1000)

        # Calculate the Gaussian PDF for the cluster
        pdf_cluster = stats.norm.pdf(x, cluster_mean, cluster_std)

        # Plot histogram
        plt.hist(cluster_data, bins=50, density=True, alpha=0.5, color=colors[i], label=f'Cluster {i} Data')

        # Plot Gaussian PDF
        plt.plot(x, pdf_cluster, color=colors[i], linestyle='-', linewidth=2, label=f'Cluster {i} Gaussian')

    plt.xlabel('Read Depth')
    plt.xlim(0,200)
    plt.ylabel('Density')
    plt.title(f'Gaussian Distributions of Read Depth for Contig {contig_id}')
    plt.legend()
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout()
    plt.show()
    

# Example use of the function for the ninth contig with 2 clusters
plot_gaussian_distributions(read_depth_data, 1, ninth_contig_id)
