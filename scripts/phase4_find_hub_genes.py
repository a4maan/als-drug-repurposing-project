import pandas as pd
import networkx as nx

# Load STRING network
network = pd.read_csv(
    "results/cluster1_string_network.csv"
)

print("Edges:", len(network))

# Build graph
G = nx.Graph()

for _, row in network.iterrows():
    G.add_edge(
        row["protein1"],
        row["protein2"]
    )

print("Nodes:", G.number_of_nodes())

# Degree centrality
centrality = nx.degree_centrality(G)

hub_df = pd.DataFrame({
    "protein": list(centrality.keys()),
    "centrality": list(centrality.values())
})

hub_df = hub_df.sort_values(
    "centrality",
    ascending=False
)

hub_df.to_csv(
    "results/cluster1_hub_genes.csv",
    index=False
)

print("\nTop 30 hubs:")
print(hub_df.head(30))