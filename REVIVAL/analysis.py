"""
A scripts to run all analyses.
"""

import numpy as np
import pandas as pd

from sklearn.metrics import ndcg_score


def ndcg_scale(y_true: np.ndarray, y_pred: np.ndarray):
    """Calculate the ndcg_score with neg correction"""

    if min(y_true) < 0:
        y_true = y_true - min(y_true)
    return ndcg_score(y_true[None, :], y_pred[None, :])


# https://github.com/OATML-Markslab/ProteinGym/blob/e1d603d28ed8c6a27959c993de30312a83203a16/proteingym/baselines/unirep/utils/metric_utils.py#L18
def ndcg(y_pred: np.array, y_true: np.array) -> float:
    """
    Compute the Normalized Discounted Cumulative Gain (NDCG) score.

    Args:
        y_pred (np.array): Predicted scores.
        y_true (np.array): True scores.

    Returns:
        float: NDCG score.
    """
    y_true_minmax = (y_true - y_true.min()) / (y_true.max() - y_true.min())  # Scale between 0 and 1
    return ndcg_score(y_true_minmax.reshape(1, -1), y_pred.reshape(1, -1))


def minmax(x):
    return (x - np.min(x)) / (np.max(x) - np.min(x))


# https://github.com/OATML-Markslab/ProteinGym/blob/main/proteingym/performance_DMS_benchmarks.py
def custom_ndcg(y_true, y_score, **kwargs):
    """
    Inputs:
        y_true: an array of the true scores where higher score is better
        y_score: an array of the predicted scores where higher score is better
    Options:
        quantile: If True, uses the top k quantile of the distribution
        top: under the quantile setting this is the top quantile to
            keep in the gains calc. This is a PERCENTAGE (i.e input 10 for top 10%)
    Notes:
        Currently we're calculating NDCG on the continuous value of the DMS
        I tried it on the binary value as well and the metrics seemed mostly
        the same.
    """
    if "quantile" not in kwargs:
        kwargs["quantile"] = True
    if "top" not in kwargs:
        kwargs["top"] = 10
    if kwargs["quantile"]:
        k = np.floor(y_true.shape[0] * (kwargs["top"] / 100)).astype(int)
    else:
        k = kwargs["top"]
    if isinstance(y_true, pd.Series):
        y_true = y_true.values
    if isinstance(y_score, pd.Series):
        y_score = y_score.values
    gains = minmax(y_true)
    ranks = np.argsort(np.argsort(-y_score)) + 1

    if k == "all":
        k = len(ranks)
    # sub to top k
    ranks_k = ranks[ranks <= k]
    gains_k = gains[ranks <= k]
    # all terms with a gain of 0 go to 0
    ranks_fil = ranks_k[gains_k != 0]
    gains_fil = gains_k[gains_k != 0]

    # if none of the ranks made it return 0
    if len(ranks_fil) == 0:
        return 0

    # discounted cumulative gains
    dcg = np.sum([g / np.log2(r + 1) for r, g in zip(ranks_fil, gains_fil)])

    # ideal dcg - calculated based on the top k actual gains
    ideal_ranks = np.argsort(np.argsort(-gains)) + 1
    ideal_ranks_k = ideal_ranks[ideal_ranks <= k]
    ideal_gains_k = gains[ideal_ranks <= k]
    ideal_ranks_fil = ideal_ranks_k[ideal_gains_k != 0]
    ideal_gains_fil = ideal_gains_k[ideal_gains_k != 0]
    idcg = np.sum(
        [g / np.log2(r + 1) for r, g in zip(ideal_ranks_fil, ideal_gains_fil)]
    )

    # normalize
    return dcg / idcg


def calc_top_n_percent_recall(y_true, y_score, top_n=10):
    """
    Compute top-N% recall where the top N% is used for both true values and predicted values.

    Args:
        y_true (array-like): Ground truth scores.
        y_score (array-like): Model predicted scores.
        top_n (int): Percentage of top-ranked items to consider.

    Returns:
        float: Recall of the model's top-N% predictions.
    """
    # Get the threshold value for the top N% of samples
    threshold_true = np.percentile(y_true, 100 - top_n)
    threshold_model = np.percentile(y_score, 100 - top_n)

    # Select the top N% based on true scores
    top_true = y_true >= threshold_true

    # Select the top N% based on model scores
    top_model = y_score >= threshold_model

    # Calculate True Positives (TP): How many correct top-N% predictions?
    true_positives = (top_true) & (top_model)

    # Compute recall: TP / Total Actual Positives (Top-N% True Samples)
    recall = true_positives.sum() / (top_true.sum()) if top_true.sum() > 0 else 0

    return recall