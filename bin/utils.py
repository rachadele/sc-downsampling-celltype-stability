import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import numpy as np
import scanpy as sc
import random
import cellxgene_census
import cellxgene_census.experimental
import os
import anndata as ad
import scvi
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import *
from sklearn.preprocessing import label_binarize
from collections import defaultdict
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import subprocess
from scipy.stats import median_abs_deviation
from statsmodels.formula.api import ols

def setup(organism="homo_sapiens", version="2024-07-01"):
    organism=organism.replace(" ", "_") 
    #census = cellxgene_census.open_soma(census_version=version)
    outdir = f"scvi-{organism}-{version}"  # Concatenate strings using f-string
    
    # Create the output directory if it doesn't exist
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Check if the model file exists
    model_file_path = os.path.join(outdir, "model.pt")
    #if not os.path.exists(model_file_path):
        # Get scVI model info
    scvi_info = cellxgene_census.experimental.get_embedding_metadata_by_name(
            embedding_name="scvi",
            organism=organism,
            census_version=version,
        )

        # Extract the model link
    model_link = scvi_info["model_link"]
    date = model_link.split("/")[5]
    url = os.path.join("https://cellxgene-contrib-public.s3.us-west-2.amazonaws.com/models/scvi/", date, organism, "model.pt")

    # Download the model using wget if it doesn't already exist
    subprocess.run(["wget", "--no-check-certificate", "-q", "-O", model_file_path, url])
# else:
     #   print(f"File already exists at {model_file_path}, skipping download.")

    return(outdir)

def subsample_and_save(dataset_path, n_cells=1000):
    dataset = ad.read_h5ad(dataset_path)
    subsampled_dataset = dataset[np.random.choice(dataset.n_obs, size=n_cells, replace=False), :] if dataset.n_obs > n_cells else dataset
    dir_name, base_name = os.path.split(dataset_path)
    file_name, ext = os.path.splitext(base_name)
    output_path = os.path.join(dir_name, f"{file_name}_subsampled_{n_cells}{ext}")
    subsampled_dataset.write_h5ad(output_path)
    print(f"Saved to {output_path}")
    
# Subsample x cells from each cell type if there are n>x cells present
#ensures equal representation of cell types in reference
def subsample_cells(data, filtered_ids, subsample=500, relabel_path=None, 
                    ref_keys=["subclass","class","family"], seed=42):
    random.seed(seed)         # For `random`
    np.random.seed(seed)      # For `numpy`
    scvi.settings.seed = seed # For `scvi`
    
    # Filter data based on filtered_ids
    obs = data[data['soma_joinid'].isin(filtered_ids)]
    relabel_df = pd.read_csv(relabel_path, sep='\t')  # Adjust the separator as needed
    # Take the first column as the join key
    join_key = relabel_df.columns[0]
    # Ensure the join_key is in both the AnnData object and the relabel DataFrame
    if join_key not in obs.columns:
        raise ValueError(f"{join_key} not found in AnnData object observations.")
    if join_key not in relabel_df.columns:
        raise ValueError(f"{join_key} not found in relabel DataFrame.")
    # Perform the left join to update the metadata
    obs = obs.merge(relabel_df, on=join_key, how='left', suffixes=(None, "_y"))
    # this line ensures that only cells in the relabel # file are kept in the obs
    # filteres out ambiguous cells
    celltypes = obs[ref_keys[0]].unique()
    final_idx = []
    for celltype in celltypes:
        celltype_ids = obs[obs[ref_keys[0]] == celltype]['soma_joinid'].values
        # Sample if there are enough observations, otherwise take all
        if len(celltype_ids) > subsample:
            subsampled_cell_idx = random.sample(list(celltype_ids), subsample)
        else:
            subsampled_cell_idx = celltype_ids.tolist()
        # Append subsampled indices to final list
        final_idx.extend(subsampled_cell_idx)

    # Return final indices
    return final_idx

def relabel(query, relabel_path, join_key=None, sep="\t"):
    # takes full adata and returns full adata
    # Read the relabel table from the file
    relabel_df = pd.read_csv(relabel_path, sep=sep)  # Adjust the separator as needed
    # Take the first column as the join key
    if join_key is None:
        join_key = relabel_df.columns[0]
    # Ensure the join_key is in both the AnnData object and the relabel DataFrame
    if join_key not in query.obs.columns:
        raise ValueError(f"{join_key} not found in AnnData object observations.")
    if join_key not in relabel_df.columns:
        raise ValueError(f"{join_key} not found in relabel DataFrame.")
    # Left join = only cell types in relabel file are kept
    query.obs = query.obs.merge(relabel_df, on=join_key, how='left', suffixes=(None, "_y"))
    columns_to_drop = [col for col in query.obs.columns if col.endswith('_y')]
    query.obs.drop(columns=columns_to_drop, inplace=True)
    return query



def aggregate_labels(query: pd.DataFrame, mapping_df: pd.DataFrame, ref_keys: list, predicted=False):
    """
    Aggregate subclass labels or predicted labels to higher levels using mapping_df.
    If predicted=True, operates on columns like 'predicted_subclass', otherwise on obs columns.
    """
    query.index = query.index.astype(int)
    for i in range(1, len(ref_keys)):
        lower_key = ref_keys[i-1]
        higher_level_key = ref_keys[i]
        mapping = mapping_df.set_index(lower_key)[higher_level_key].to_dict()
        if predicted:
            query["predicted_" + higher_level_key] = query["predicted_" + lower_key].map(mapping)
            query["predicted_" + higher_level_key] = query["predicted_" + higher_level_key].fillna(query["predicted_" + lower_key])
        else:
            query[higher_level_key] = query[lower_key].map(mapping)
            query[higher_level_key] = query[higher_level_key].fillna(query[lower_key])
    return query


def extract_data(data, filtered_ids, subsample=10, organism=None, census=None, 
    cell_columns=None, dataset_info=None, dims=20, relabel_path="biof501_proj/meta/relabel/census_map_human.tsv'", 
    ref_keys=["rachel_subclass","rachel_class","rachel_family"], original_celltypes=None, seed=42):
    
    subsampled_ids = subsample_cells(data, filtered_ids, subsample, relabel_path=relabel_path, ref_keys=ref_keys, seed=seed)
    adata = cellxgene_census.get_anndata(
        census=census,
        organism=organism,
        obs_column_names=cell_columns,
        obs_coords=subsampled_ids,
        var_value_filter = "nnz > 10",
        obs_embeddings=["scvi"])
   
    # Filter the AnnData object to include only the specified cell types 
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.filter_genes(adata, min_counts=200)

    print("Subsampling successful.")
    newmeta = adata.obs.merge(dataset_info, on="dataset_id", suffixes=(None,"y"))
    adata.obs = newmeta
    
    # before relabeling, need to map back to author types using new_observation_joinid
    # huge pain in the ass
    # only do this if original celltypes passed
    if isinstance(original_celltypes, pd.DataFrame) and not original_celltypes.empty: 
        adata.obs = map_author_labels(adata.obs, original_celltypes)
    # Assuming relabel_wrapper is defined
    adata = relabel(adata, relabel_path=relabel_path, sep='\t')
    # Convert all columns in adata.obs to factors (categorical type in pandas)
    return adata

def split_and_extract_data(data, split_column, subsample=500, organism=None, census=None, 
                           cell_columns=None, dataset_info=None, dims=20, relabel_path="/biof501_proj/meta/relabel/census_map_human.tsv",
                           ref_keys=["rachel_subclass","rachel_class","rachel_family"], seed=42, original_celltypes=None):
    # Get unique split values from the specified column
    unique_values = data[split_column].unique()
    refs = {}
    # don't use Human Variation Study splits
    
    for split_value in unique_values:
        # Filter the brain observations based on the split value
        
        filtered_ids = data[data[split_column] == split_value]['soma_joinid'].values
        collection_names = data[data[split_column] == split_value]['collection_name'].unique()
        if split_column == "dataset_id" and "HVS: Human variation study" in collection_names:
            # send stderr message
            print(f"Skipping {split_value} due to collection name 'HVS: Human variation study'.")
            continue
        if len(filtered_ids) < 1000:
            print(f"Skipping {split_value} due to insufficient cells ({len(filtered_ids)} < 1000).")
            continue
       # obs_filter = f"{split_column} == '{split_value}'"
       # don't need this since we're subsampling the extact cell coordinates
        print(f"Processing {split_value} with {len(filtered_ids)} cells.")
        adata = extract_data(data, filtered_ids, subsample, organism, census, 
                             cell_columns, dataset_info, dims=dims, relabel_path=relabel_path, 
                             ref_keys=ref_keys, 
                             original_celltypes = original_celltypes, 
                             seed=seed)
        dataset_titles = adata.obs['dataset_title'].unique()

        if split_column == "tissue": 
            name_to_use = split_value
        elif split_column == "dataset_id":
            name_to_use = dataset_titles[0]
        else:
            name_to_use = split_value

        refs[name_to_use] = adata

    return refs


def get_original_celltypes(columns_file="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/author_cell_annotations/2025-01-30/original_celltype_columns.tsv", 
                           author_annotations_path="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/author_cell_annotations/2025-01-30"):
    original_celltype_columns = pd.read_csv(columns_file, sep="\t", low_memory=False)

    original_celltypes = {}
    for file in os.listdir(author_annotations_path):
        if "obs.tsv" in file:
            dataset_title = file.split(".")[0]
            og_obs = pd.read_csv(os.path.join(author_annotations_path, file), sep="\t", low_memory=False)
            # check if all observation_joinid are unique
            assert og_obs["observation_joinid"].nunique() == og_obs.shape[0]
            og_column = original_celltype_columns[original_celltype_columns["dataset_title"] == dataset_title]["author_cell_type"].values[0]
            og_obs["author_cell_type"] = og_obs[og_column]
            original_celltypes[dataset_title] = og_obs
            original_celltypes[dataset_title]["new_dataset_title"] = dataset_title
            
    for dataset_title, obs in original_celltypes.items():
        original_celltypes[dataset_title]["new_observation_joinid"] = original_celltypes[dataset_title]["observation_joinid"].apply(lambda x: f"{dataset_title}_{x}")
    
        # concat all original_celltypes
    aggregate_obs = pd.concat([original_celltypes[ref_name] for ref_name in original_celltypes.keys()])
    # prevent duplicate observation_joinid in aggregate_obs
    assert aggregate_obs["new_observation_joinid"].nunique() == aggregate_obs.shape[0]
    
    return aggregate_obs

def map_author_labels(obs, original_celltypes):
    obs["new_dataset_title"] = obs["dataset_title"].apply(lambda x: x.replace(" ", "_")
                                                                .replace("\\/", "_")
                                                                .replace("(", "")
                                                                .replace(")", "")
                                                                .replace("\\", "")
                                                                .replace("'", "")
                                                                .replace(":", "")
                                                                .replace(";", "")
                                                                .replace("&", ""))

    obs["new_observation_joinid"] = obs["new_dataset_title"].astype(str) + "_" + obs["observation_joinid"].astype(str)
    
    mapping = dict(zip(original_celltypes["new_observation_joinid"], original_celltypes["author_cell_type"]))
    obs["author_cell_type"] = obs["new_observation_joinid"].map(mapping)

    return obs


def get_filtered_obs(census, organism, organ="brain", is_primary=True, disease="normal"):
    value_filter = (
        f"tissue_general == '{organ}' and "
        f"is_primary_data == {str(is_primary)} and "
        f"disease == '{disease}'"
    )
    return cellxgene_census.get_obs(census, organism, value_filter=value_filter)

def get_census(census_version="2024-07-01", organism="homo_sapiens", subsample=5, split_column="dataset_id", dims=50, organ="brain",
               ref_collections=["Transcriptomic cytoarchitecture reveals principles of human neocortex organization"," SEA-AD: Seattle Alzheimer’s Disease Brain Cell Atlas"],
               relabel_path="../meta/census_map_human.tsv", 
               seed=42, 
               ref_keys=["rachel_subclass","rachel_class","rachel_family"],
               original_celltypes=None):

    census = cellxgene_census.open_soma(census_version=census_version)
    dataset_info = census.get("census_info").get("datasets").read().concat().to_pandas()
    cellxgene_obs = get_filtered_obs(census, organism, organ=organ, is_primary=True, disease="normal")
    
    cellxgene_obs = cellxgene_obs.merge(dataset_info, on="dataset_id", suffixes=(None,"_y"))
    cellxgene_obs.drop(columns=['soma_joinid_y'], inplace=True)
    cellxgene_obs_filtered = cellxgene_obs[cellxgene_obs['collection_name'].isin(ref_collections)] 

    # Adjust organism naming for compatibility
    organism_name_mapping = {
        "homo_sapiens": "Homo sapiens",
        "mus_musculus": "Mus musculus"
    }
    organism = organism_name_mapping.get(organism, organism)

    cell_columns = [
        "assay", "cell_type", "tissue",
        "tissue_general", "suspension_type",
        "disease", "dataset_id", "development_stage",
        "soma_joinid","observation_joinid"
    ]
    
    # add author_cell_type to obs
    # this will enable proper relabeling and subsampling
    # need to add it back in after getting ids
    if isinstance(original_celltypes, pd.DataFrame) and not original_celltypes.empty:
        cellxgene_obs_filtered = map_author_labels(cellxgene_obs_filtered, original_celltypes)
    # Get individual datasets and embeddings
    refs = split_and_extract_data(
        cellxgene_obs_filtered, split_column=split_column,
        subsample=subsample, organism=organism,
        census=census, cell_columns=cell_columns,
        dataset_info=dataset_info, dims=dims,
        relabel_path=relabel_path,
        ref_keys=ref_keys, seed=seed, 
        original_celltypes=original_celltypes
    )
    # Get embeddings for all data together
    filtered_ids = cellxgene_obs_filtered['soma_joinid'].values
    adata = extract_data(
        cellxgene_obs_filtered, filtered_ids,
        subsample=subsample, organism=organism,
        census=census,
        cell_columns=cell_columns, 
        dataset_info=dataset_info, dims=dims,
        relabel_path=relabel_path, 
        ref_keys=ref_keys, seed = seed, 
        original_celltypes=original_celltypes
    )
    refs["whole cortex"] = adata
    for name, ref in refs.items():
        dataset_title = name.replace(" ", "_")
        for col in ref.obs.columns:
            if ref.obs[col].dtype.name =='category':
    # Convert to Categorical and remove unused categories
                ref.obs[col] = pd.Categorical(ref.obs[col].cat.remove_unused_categories())
     
    return refs



def process_query(query, model_file_path, batch_key="sample", seed=42):
    scvi.settings.seed = seed # For `scvi`
    # Ensure the input AnnData object is valid
    if not isinstance(query, ad.AnnData):
        raise ValueError("Input must be an AnnData object.")

    if "feature_id" in query.var.columns:
        query.var.set_index("feature_id", inplace=True)

    query.obs["n_counts"] = query.X.sum(axis=1)
    query.obs["joinid"] = list(range(query.n_obs))
    query.obs["batch"] = query.obs[batch_key]

    # Prepare the query AnnData for scVI
    scvi.model.SCVI.prepare_query_anndata(query, model_file_path)
    vae_q = scvi.model.SCVI.load_query_data(query, model_file_path)

    # Set the model to trained and get latent representation
    vae_q.is_trained = True
    latent = vae_q.get_latent_representation()
    query.obsm["scvi"] = latent

    return query


def map_valid_labels(query, ref_keys, mapping_df):
    # deal with differing levels of granularity
    # Loop over each reference label (e.g., "celltype", "subclass", "class", etc.)
    for key in ref_keys:
        # Get all unique labels in the query at this level (e.g., all unique 'subclass' values)
        original=query[key].unique()
        for og in original: # For each original label at this level
            # Get the highest level in the hierarchy 
            # Search for columns in mapping_df where this label appears in *any* column
            matching_cols = mapping_df.columns[mapping_df.apply(lambda col: og in col.values, axis=0)]
            print(f"Matching columns for {og}: {matching_cols}")
            if len(matching_cols) == 0:
                continue  # likely "unknown", skip
            else:
                level = matching_cols[-1]
                # Check if level is above key in the hierarchy
                if mapping_df.columns.get_loc(level) > mapping_df.columns.get_loc(key):
                    print(f"Level {level} is above level {key} in the hierarchy.")        
                    og_index = query.index[query[key] == og]
                    # Replace the value in "predicted_" column with corresponding predicted value at `level`
                    for idx in og_index:
                        # Find the replacement value from `mapping_df` for this level
                        replacement = query.loc[idx, "predicted_" + level]
                        print(f"Replacing predictions for {og} with {replacement} to match {level}")
                        # replace predicted id with appropriate level
                        query["predicted_" + key] = query["predicted_" + key].astype("object")
                        query.loc[idx, "predicted_" + key] = replacement#.iloc[0]
                        query["predicted_" + key] = query["predicted_" + key].astype("category")

    return query            


def rfc_pred(ref, query, ref_keys, seed):
    """
    Fit a RandomForestClassifier at the most granular level and aggregate probabilities for higher levels.
    
    Parameters:
    - ref: Reference data with labels.
    - query: Query data for prediction.
    - ref_keys: List of ordered keys from most granular to highest level (e.g., ["rachel_subclass", "rachel_class", "rachel_family"]).
    - tree: Dictionary representing the hierarchy of classes.
    
    Returns:
    - probabilities: Dictionary with probabilities for each level of the hierarchy.
    """
    probabilities = {}
    
    # The most granular key is the first in the ordered list
    granular_key = ref_keys[0]
    
    # Initialize and fit the RandomForestClassifier at the most granular level
    rfc = RandomForestClassifier(class_weight='balanced', random_state=seed)
    rfc.fit(ref.obsm["scvi"], ref.obs[granular_key].values)
    # Predict probabilities at e most granular level
    probs_granular = rfc.predict_proba(query.obsm["scvi"])
    class_labels_granular = rfc.classes_
    query.obs[granular_key] = query.obs[granular_key].astype(str)
    base_score = rfc.score(query.obsm["scvi"], query.obs[granular_key].values)

    # Store granular level probabilities
    probabilities[granular_key] = {
        "probabilities": probs_granular,
        "class_labels": class_labels_granular,
        "accuracy": base_score
    }
    
    return probabilities


def roc_analysis(probabilities, query, key):
    optimal_thresholds = {}
    metrics={}
  #  for key in ref_keys:
       # print(key) 
    #probs = probabilities[key]["probabilities"]
    probs = np.array(probabilities)
    class_labels = probabilities.columns.values
    optimal_thresholds[key] = {}
        
    # Binarize the class labels for multiclass ROC computation
    true_labels = label_binarize(query[key].values, classes=class_labels)
        
    # Find the optimal threshold for each class
    metrics[key] = {}
    for i, class_label in enumerate(class_labels):
        optimal_thresholds[key][class_label] = {}
        # check for positive samples
        # usually positive samples are 0 when a ref label is
        # replaced with a parent label
        # since it is not in the original query labels
        # but it is being annotated during the label transfer
        # these should not be evaluated ?
        # incorrect classifications will be reflected in the AUC and F1 of the og label
        # eg. ET is not in query so it is changed to "deep layer non-IT"
        # but these cells are CT or NP in the ref, so they're still incorrect
        # not entirely sure how to deal with this
        positive_samples = np.sum(true_labels[:, i] == 1)
        if positive_samples == 0:
            print(f"Warning: No positive samples for class {class_label}, skipping eval and setting threshold to 0.5")
            optimal_thresholds[key][class_label] = 0.5
        elif positive_samples > 0:
            metrics[key][class_label]={}
            # True label one hot encoding at class label index = 
            # vector of all cells which are either 1 = label or 0 = not label
            # probs = probability vector for all cells given class label
            fpr, tpr, thresholds = roc_curve(true_labels[:, i], probs[:, i])
            roc_auc = auc(fpr, tpr) # change to roc_auc_score, ovo, average= macro, labels               
            optimal_idx = np.argmax(tpr - fpr)
            optimal_threshold = thresholds[optimal_idx]
            if optimal_threshold == float('inf'):
                optimal_threshold = 0 
            optimal_thresholds[key][class_label]=optimal_threshold
            metrics[key][class_label]["tpr"] = tpr
            metrics[key][class_label]["fpr"] = fpr
            metrics[key][class_label]["auc"] = roc_auc
            metrics[key][class_label]["optimal_threshold"] = optimal_threshold
    return metrics


def pr_analysis(prob_df, query, key, mapping_df):

    metrics = defaultdict(dict)

    probs = np.array(prob_df)
    class_labels = prob_df.columns.values         
    original_labels = query[key].unique()
    
    # map back to correct granularity
    #for label in class_labels:
        #if label not in original_labels:
            
            #parent_label = mapping_df.loc[mapping_df[key] == label, key].values[0]
            #if parent_label:
                #probs[:, class_labels == label] = probs[:, class_labels == parent
                #_label]
    
    # Binarize the class labels for multiclass ROC computation 
    
    true_labels = label_binarize(query[key].values, classes=class_labels)
    
    # Find the optimal threshold for each class

    for i, class_label in enumerate(class_labels):
        positive_samples = np.sum(true_labels[:, i] == 1)
        if positive_samples == 0:
            print(f"Warning: No positive samples for class {class_label}, skipping eval")
        elif positive_samples > 0:
            precision, recall, thresholds = precision_recall_curve(true_labels[:, i], probs[:, i])
            avg_precision = average_precision_score(true_labels[:, i], probs[:, i])

            # Find optimal threshold as the one maximizing F1 score
            f1_scores = 2 * (precision * recall) / (precision + recall + 1e-10)
            optimal_idx = np.argmax(f1_scores)

            # Store metrics
            metrics["precision"][class_label] = precision
            metrics["recall"][class_label] = recall
            metrics["average_precision"][class_label] = avg_precision
            metrics["optimal_threshold"][class_label] = thresholds[optimal_idx]
    return metrics

def check_column_ties(probabilities, class_labels):
    """
    Checks for column ties (multiple columns with the same maximum value) in each row.

    Parameters:
    - probabilities (numpy.ndarray): 2D array where rows represent samples and columns represent classes.
    - class_labels (list): List of class labels corresponding to the columns.

    Returns:
    - tie_rows (list): List of row indices where ties occur.
    - tie_columns (dict): Dictionary where the key is the row index and the value is a list of tied column indices and their class labels.
    """
    # Find the maximum probability for each row
    max_probs = np.max(probabilities, axis=1)

    # Check if there are ties (multiple columns with the maximum value in the row)
    ties = np.sum(probabilities == max_probs[:, None], axis=1) > 1

    # Get the indices of rows with ties
    tie_rows = np.where(ties)[0]

    # Find the columns where the tie occurs and associated class labels
    tie_columns = {}
    for row in tie_rows:
        tied_columns = np.where(probabilities[row] == max_probs[row])[0]
        tie_columns[row] = [(col, class_labels[col]) for col in tied_columns]
    
    return tie_rows, tie_columns


def compute_confidence_gap_cutoff(class_probs, quantile=0.1):
    """
    Compute confidence gap cutoff at given quantile.

    Args:
        probs (np.array): shape (n_cells, n_classes), normalized probabilities.
        quantile (float): quantile for cutoff.

    Returns:
        float: confidence gap cutoff.
    """
    sorted_probs = -np.sort(-class_probs, axis=1)
    gaps = sorted_probs[:, 0] - sorted_probs[:, 1]
    cutoff = np.quantile(gaps, quantile)
    return cutoff



def classify_cells(query: pd.DataFrame, ref_keys: list, cutoff: float, probabilities: pd.DataFrame, mapping_df: pd.DataFrame, use_gap: bool = False):
    """
    Wrapper function to classify cells based on probabilities. Can use gap analysis or raw probabilities.
    """
    class_metrics = {}
    
    # Only use the first ref_key
    key = ref_keys[0]
    class_metrics[key] = {}

    # Extract the class labels and probabilities (array structure)
    class_labels = probabilities.columns.values  # Class labels are the column names
    class_probs = probabilities.values  # Probabilities as a numpy array
    
    if use_gap:
        predictions = classify_by_gap(class_probs, class_labels, cutoff=cutoff)
    else:
        predictions = classify_raw(class_probs, class_labels, cutoff=cutoff)
    
    # Store predictions and confidence in `query`
    query["predicted_" + key] = predictions
    query["confidence"] = np.max(class_probs, axis=1)  # Store max probability as confidence
    
    # Aggregate predictions 
    query = aggregate_labels(query=query, ref_keys=ref_keys, mapping_df=mapping_df, predicted=True)
    
    return query


def classify_raw(class_probs, class_labels, cutoff=0):
  """
  Classify cells based on raw probabilities without gap analysis.
  """ 
  predictions = []
  if cutoff > 0:
    max_class_indices = np.argmax(class_probs, axis=1)  # Get the index of the max probability
    max_class_probs = np.max(class_probs, axis=1)  # Get the max probability

    # Set predicted classes to "unknown" if the max probability does not meet the threshold
    predictions = [
        class_labels[i] if prob > cutoff else "unknown"
        for i, prob in zip(max_class_indices, max_class_probs)
    ]
  else:
    # Direct prediction without threshold filtering
    predictions = class_labels[np.argmax(class_probs, axis=1)]

  return predictions

    
def classify_by_gap(class_probs, class_labels, cutoff=0.1): # gap must be 0.1 or higher
    # potentially normalize probabilities ? already between 0 and 1 but the distribution will differ based on type of classifier
    """
    Classify cells as a predicted type or 'unknown' based on confidence gap.

    Args:
        class_probs (np.array): shape (n_cells, n_classes), normalized probabilities.
        class_labels (list): list of class labels.
        cutoff (float): confidence gap cutoff.

    Returns:
        np.array: array of predicted classes ('known' or 'unknown').
    """
    predictions = []
    sorted_probs = -np.sort(-class_probs, axis=1)
    gaps = sorted_probs[:, 0] - sorted_probs[:, 1]
    
    # Classify based on gap
    predictions = np.where(gaps < cutoff, "unknown", class_labels[np.argmax(class_probs, axis=1)])
    
    return predictions

def aggregate_preds(query, ref_keys, mapping_df):
    preds = np.array(query["predicted_" + ref_keys[0]])
    query.index = query.index.astype(int)

    # Reorder ref_keys based on the order in tree_df columns (ignoring "cell_type")
    ref_keys = [col for col in mapping_df.columns if col in ref_keys]

    for higher_level_key in ref_keys[1:]:  # Skip the first (most granular) level
        # Get mapping from subclass to higher-level class
        mapping = mapping_df.set_index(ref_keys[0])[higher_level_key].to_dict()

        # Assign higher-level labels based on mapping
        query["predicted_" + higher_level_key] = query["predicted_" + ref_keys[0]].map(mapping)

        # Fill NA values with original subclass labels
        query["predicted_" + higher_level_key] = query["predicted_" + higher_level_key].fillna(query["predicted_" + ref_keys[0]])

    return query

def evaluate_sample_predictions(query, ref_keys, mapping_df):
    class_metrics = defaultdict(lambda: defaultdict(dict))

    for key in ref_keys:     
        true_labels = query[key].astype(str)
        predicted_labels = query[f"predicted_{key}"].astype(str)
        labels = list(set(true_labels).union(set(predicted_labels)))

        # Overall accuracy
        class_metrics[key]["overall_accuracy"] = accuracy_score(true_labels, predicted_labels)

        # Confusion matrix
        class_metrics[key]["confusion"] = {
            "matrix": confusion_matrix(true_labels, predicted_labels, labels=labels),
            "labels": labels
        }

        # NMI, ARI between cell type labels
        class_metrics[key]["nmi"] = normalized_mutual_info_score(true_labels, predicted_labels, average_method='arithmetic')
        class_metrics[key]["ari"] = adjusted_rand_score(true_labels, predicted_labels)

        # Per-label metrics
        precision, recall, f1, support = precision_recall_fscore_support(
            true_labels, predicted_labels, labels=labels, zero_division=np.nan
        )
        support_proportions = support / support.sum() if support.sum() > 0 else np.zeros_like(support)

        class_metrics[key]["label_metrics"] = {
            label: {
                "precision": precision[i],
                "recall": recall[i],
                "f1_score": f1[i],
                "support": support_proportions[i],
                "accuracy": accuracy_score(true_labels == label, predicted_labels == label)
            }
            for i, label in enumerate(labels)
        }

        # Weighted averages
        avg_p, avg_r, avg_f, _ = precision_recall_fscore_support(
            true_labels, predicted_labels, average="weighted", zero_division=np.nan
        )
        class_metrics[key]["weighted_metrics"] = {
            "precision": avg_p,
            "recall": avg_r,
            "f1_score": avg_f,
        }
        
        # Macro averages
        macro_p, macro_r, macro_f, _ = precision_recall_fscore_support(
            true_labels, predicted_labels, average="macro", zero_division=np.nan
        )
        class_metrics[key]["macro_metrics"] = {
            "precision": macro_p,
            "recall": macro_r,
            "f1_score": macro_f,
        }

    return class_metrics


def update_classification_report(class_metrics, ref_keys):
    #for query_name, query in class_metrics.items():
       # for ref_name, ref in query.items():
    for key in ref_keys:
        for key, val in class_metrics[key]["classification_report"].items():
            if isinstance(val, dict):
                if val['support'] == 0.0: 
                    val['f1-score'] = "nan" 
    return class_metrics

def plot_confusion_matrix(query_name, ref_name, key, confusion_data, output_dir):

    #new_query_name = query_name.replace(" ", "_").replace("/", "_").replace("(","").replace(")","")
   # new_ref_name = ref_name.replace(" ", "_").replace("/", "_").replace("(","").replace(")","")
 #   output_dir = os.path.join(new_query_name, new_ref_name, output_dir)
    os.makedirs(output_dir, exist_ok=True) 
    # Extract confusion matrix and labels from the confusion data
    conf_matrix = confusion_data["matrix"]
    labels = confusion_data["labels"]

    # Plot the confusion matrix
    plt.figure(figsize=(28, 12))
    sns.heatmap(conf_matrix, annot=True, fmt='g', cmap='Reds', xticklabels=labels, yticklabels=labels, annot_kws={"size": 20})
    plt.xlabel('Predicted', fontsize =20)
    plt.ylabel('True', fontsize= 20)
    plt.title(f'Confusion Matrix: {query_name} vs {ref_name} - {key}', fontsize=17)
    # Adjust tick parameters for both axes
    #plt.tick_params(axis='both', which='major', labelsize=15, width=1)  # Increase tick label size and make ticks thicker

    # Rotate both x and y tick labels by 90 degrees
    plt.xticks(rotation=45, fontsize=15)  # Rotate x-axis labels by 90 degrees
    plt.yticks(rotation=45, fontsize=15)  # Rotate y-axis labels by 90 degrees
    
    #os.makedirs(os.path.join(output_dir, new_query_name, new_ref_name), exist_ok=True)  # Create the directory if it doesn't exist
    plt.savefig(os.path.join(output_dir,f"{key}_confusion.png"))
    plt.close() 
     
  
def map_genes(query, gene_mapping):
    # Drop rows with missing values in the relevant columns
    gene_mapping = gene_mapping.dropna(subset=["ENSEMBL_ID", "OFFICIAL_SYMBOL"])
    # Set the index of gene_mapping to "ENSEMBL_ID" and ensure it's unique
    gene_mapping = gene_mapping.drop_duplicates(subset="ENSEMBL_ID")
    gene_mapping.set_index("ENSEMBL_ID", inplace=True)
        # Check if the "feature_name" column exists in query.var
    if "feature_name" not in query.var.columns:
        # Merge gene_mapping with query.var based on the index
        query.var = query.var.merge(gene_mapping["OFFICIAL_SYMBOL"], left_index=True, right_index=True, how="left")
        # Rename the merged column to "feature_name"
        query.var.rename(columns={"OFFICIAL_SYMBOL": "feature_name"}, inplace=True)
    return query


def is_outlier(query, metric: str, nmads=3):
    M = query.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


def get_lm(query, nmads=5, scale="normal"):
    # Assume dataset is an AnnData object
    # Fit linear model: log10(n_genes_per_cell) ~ log10(counts_per_cell)
    lm_model = ols(formula='log1p_n_genes_by_counts ~ log1p_total_counts', data=query.obs).fit()
    # Calculate residuals
    residuals = lm_model.resid
    # If data is normally distributed, this is similar to std 
    mad_residuals = median_abs_deviation(residuals, scale=scale)
    # Intercept adjustment (add for upper bound, subtract for lower bound)
    intercept_adjustment = np.median(residuals) + nmads * mad_residuals
    return {
        "model": lm_model,
        "intercept_adjustment": intercept_adjustment
    }
    
     
def get_qc_metrics(query, nmads):
    query.var["mito"] = query.var["feature_name"].str.startswith(("MT", "mt", "Mt"))
    query.var["ribo"] = query.var["feature_name"].str.startswith(("RP", "Rp", "rp"))
    query.var["hb"] = query.var["feature_name"].str.startswith(("HB", "Hb","hb"))
    # fill NaN values with False
    query.var["mito"].fillna(False, inplace=True)
    query.var["ribo"].fillna(False, inplace=True)
    query.var["hb"].fillna(False, inplace=True) 

    sc.pp.calculate_qc_metrics(query, qc_vars=["mito", "ribo", "hb"], log1p=True, inplace=True, percent_top=[20])

    metrics = {
        "log1p_total_counts": "umi_outlier",
        "log1p_n_genes_by_counts": "genes_outlier",
        "pct_counts_mito": "mito_outlier",
        "pct_counts_ribo": "ribo_outlier",
        "pct_counts_hb": "hb_outlier",
    }
    
    for metric, col_name in metrics.items():
        query.obs[col_name] = is_outlier(query, metric, nmads)

    lm_dict = get_lm(query, nmads=nmads)
    intercept = lm_dict["model"].params[0]
    slope = lm_dict["model"].params[1]
    

    query.obs["counts_outlier"] = (
        query.obs["log1p_n_genes_by_counts"] < (query.obs["log1p_total_counts"] * slope + (intercept - lm_dict["intercept_adjustment"]))
        ) | (
        query.obs["log1p_n_genes_by_counts"] > (query.obs["log1p_total_counts"] * slope + (intercept + lm_dict["intercept_adjustment"]))
        ) #| (
       # query.obs["umi_outlier"] ) | (query.obs["genes_outlier"])
        

    query.obs["total_outlier"] = (
        query.obs["counts_outlier"] | query.obs["mito_outlier"] | query.obs["ribo_outlier"] | query.obs["hb_outlier"] 
        | query.obs["predicted_doublet"] | query.obs["genes_outlier"] | query.obs["umi_outlier"]
    )
    
    query.obs["non_outlier"] = ~query.obs["total_outlier"]

    return query

def qc_preprocess(query, target_sum=1e4, n_top_genes=2000, n_neighbors=10, n_pcs=30, resolution=0.3):
    sc.pp.normalize_total(query, target_sum=target_sum)
    sc.pp.log1p(query)
    sc.pp.highly_variable_genes(query, n_top_genes=n_top_genes, subset=False)
    sc.pp.pca(query)
    sc.pp.neighbors(query, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(query)
    sc.tl.leiden(query, resolution=resolution)
    
    return query

def map_celltype_hierarchy(query, markers_file):
    # Load the markers table
    df = pd.read_csv(markers_file, sep=None, header=0)
    df.drop(columns="markers", inplace=True)
    query.obs = query.obs.merge(df, left_on="subclass", right_on="subclass", how="left", suffixes=("", "_y"))
    return query

def get_gene_to_celltype_map(df, organism="mus_musculus"):
    # Read the marker file
    df = df[df["markers"].notnull()]
    gene_to_celltype = {}

    for _, row in df.iterrows():
        cell_type = row["shortname"]
        genes = [gene.strip() for gene in row["markers"].split(",")]
        for gene in genes:
            if organism == "mus_musculus":
                gene = gene.lower().capitalize()
                # Handle multiple cell types mapping to the same gene
            if gene not in gene_to_celltype:
                gene_to_celltype[gene] = []
            gene_to_celltype[gene].append(cell_type)

    # Join multiple cell types into one label if needed
    gene_ct_dict = {
        gene: f"{'_'.join(set(celltypes))}: {gene}"
        for gene, celltypes in gene_to_celltype.items()
    }
    return gene_ct_dict


def make_celltype_matrices(query, markers_file, organism="mus_musculus", outdir=""):
    # Drop vars with NaN feature names
    query = query[:, ~query.var["feature_name"].isnull()]
    query.var_names = query.var["feature_name"]
    
    markers_df = pd.read_csv(markers_file, sep="\t")
    markers_df = markers_df[markers_df["organism"] == organism]
    ontology_mapping = markers_df.set_index("cell_type")["shortname"].to_dict()
    query.raw.var.index = query.raw.var["feature_name"]

    # Read marker genes
    gene_ct_dict = get_gene_to_celltype_map(markers_df, organism=organism)
    # Collect all unique markers across all families/classes/cell types
    all_markers = list(gene_ct_dict.keys())
    valid_markers = [gene for gene in all_markers if gene in query.var_names]
    removed_markers = [gene for gene in all_markers if gene not in query.var_names]
    # Write removed markers to a text file, one per line
    with open("removed_markers.txt", "w") as f:
        for gene in removed_markers:
            f.write(f"{gene}\n")
    # Filter raw expression matrix to match query.var_names
    expr_matrix = query.raw.X.toarray()
    expr_matrix = pd.DataFrame(expr_matrix, index=query.obs.index, columns=query.raw.var.index)
    
    avg_expr = expr_matrix.groupby(query.obs["predicted_subclass"]).mean()
    avg_expr = avg_expr.loc[:, valid_markers]
    
    # Scale expression across genes
    scaled_expr = (avg_expr - avg_expr.mean()) / avg_expr.std()
    scaled_expr = scaled_expr.loc[:, valid_markers]
    scaled_expr.fillna(0, inplace=True)

    # Rename columns: gene -> gene (celltype)
    scaled_expr.rename(columns=gene_ct_dict, inplace=True)
    sorted_columns = sorted(scaled_expr.columns, key=lambda x: x.split(":")[0])  
    
    # Sort by the first part of the column name
    scaled_expr = scaled_expr[sorted_columns]

    # Save matrix
    os.makedirs(outdir, exist_ok=True)
    scaled_expr.to_csv(f"{outdir}/heatmap_mqc.tsv", sep="\t")


    
def classify_by_entropy(class_probs, class_labels, quantile):
    
    if quantile == 0:
        # If quantile is 0, return the most probable class for each cell
        predictions = class_labels[np.argmax(class_probs, axis=1)]
        return predictions
    # class_probs: (n_cells, n_classes), rows sum to 1
    entropies = entropy(class_probs.T, base=2)  # or base=np.e
    # Use a cutoff, e.g., 90th percentile of low entropy = confident
    cutoff = np.quantile(entropies, 1-quantile) # top 10% entropy are unknown
   # entropies = entropy(class_probs.T, base=2)  # or base=np.e
    predictions = np.where(entropies <= cutoff, class_labels[np.argmax(class_probs, axis=1)], "unknown")
    return predictions

def plot_entropy_distribution(class_probs, class_labels, entropies, predictions, cutoff):
    df = pd.DataFrame({
        "predicted_class": predictions,
        "max_prob": np.max(class_probs, axis=1),
        "entropy": entropies,
    })
    
    # histogram of max probabilities
    # facet by predicted class
    sns.histplot(
        data=df,
        x="entropy",
        hue="predicted_class",
        multiple="layer",
        bins=50,
        kde=False,
        palette="tab20"  # optional for better class separation
    )

    plt.axvline(cutoff, color='red', linestyle='--', label=f'Cutoff: {cutoff:.3f}')
    plt.xlabel("Entropy of Probabilities")
    plt.ylabel("Cell Count")
    plt.tight_layout()
    


def plot_max_prob_hist_by_gap(class_probs, cutoff, quantile=0.1):
    """
    Plot max predicted probabilities histogram colored by known/unknown classification.

    Args:
        probs (np.array): (n_cells, n_classes), normalized probabilities.
        cutoff (float): confidence gap cutoff.
    """
    sorted_probs = -np.sort(-class_probs, axis=1)
    gaps = sorted_probs[:, 0] - sorted_probs[:, 1]
    max_probs = sorted_probs[:, 0]

    unknown_mask = gaps < cutoff
    known_mask = ~unknown_mask

    plt.figure(figsize=(8, 5))
    plt.hist(max_probs[known_mask], bins=50, alpha=0.7, label=f'Known (gap < {cutoff:.3f})', color='blue', density=False)
    plt.hist(max_probs[unknown_mask], bins=50, alpha=0.7, label=f'Unknown (gap > {cutoff:.3f})', color='red', density=False)
    plt.text(0.05, 0.95, f'Quantile: {quantile}', fontsize=12, ha='left', va='top', transform=plt.gca().transAxes)
    plt.xlabel('Max predicted probability')
    plt.ylabel('Count')
    plt.title('Max Probability Distribution by Confidence Gap Classification')
    plt.legend()
    plt.savefig('max_prob_hist_by_gap.png')
    plt.close()
    
