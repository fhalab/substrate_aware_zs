"""
A script to run all analyses.
"""
# Linear combo from Chenghao and Lukas

import numpy as np
import pandas as pd

from sklearn.metrics import ndcg_score
from sklearn.preprocessing import StandardScaler

import numpy as np
import pandas as pd
import os
import json
from glob import glob
from tqdm import tqdm

from copy import deepcopy

from sklearn.linear_model import LinearRegression
from scipy.optimize import minimize
from scipy.stats import spearmanr

import seaborn as sns
import matplotlib.pyplot as plt


from REVIVAL.util import checkNgen_folder, get_file_name


COMMON_COLS = [
    "lib",
    "hd",
    "ev_score",
    "esm_score",
    "esmif_score_apo-score",
    "coves_score_apo_clean-output-100_processed",
    "Triad_score_score-frompdb-cleanup",
    "ligandmpnn_score",
    "flowsite_score",
    "dH",  # "complexscore",
    "vina_apo-score-substrate_cofactor-docked",
    "var_vol",
]
COMMON_HEME_COLS = COMMON_COLS + [
    "chain_iptm_BA_avg_score_seperate_chai",
    "chain_pae_min_CA_avg_score_seperate_af3",
    "pocket-subcofcentroid-hw_avg - substrate-logp_af3-struct_separate",
    "num_hydrogen_bond_avg_af3-score_seperate",
    # "num_interactions_avg_af3-score_seperate"
]


COMMON_COL_DICT = {
    "lib": "Library",
    "hd": "Hamming distance",
    "ev_score": "EVmutation",
    "esm_score": "ESM2",
    "esmif_score_apo-score": "ESM-IF",
    "coves_score_apo_clean-output-100_processed": "CoVES",
    "Triad_score_score-frompdb-cleanup": r"ΔΔ$G_f$",
    "ligandmpnn_score": "LigandMPNN",
    "flowsite_score": "FlowSite",
    "dH": "GALigandDock",
    # "complexscore": "GALigandock",
    "vina_apo-score-substrate_cofactor-docked": "Vina",
    "var_vol": "Active-site volume",
}
COMMON_HEME_COL_dict = {
    **COMMON_COL_DICT,
    "chain_iptm_BA_avg_score_seperate_chai": "Chai-1",
    "chain_pae_min_CA_avg_score_seperate_af3": "AF3",
    "pocket-subcofcentroid-hw_avg - substrate-logp_af3-struct_separate": "Hydrophobicity",
    "num_hydrogen_bond_avg_af3-score_seperate": "Hydrogen bonds",
    # "num_interactions_avg_af3-score_seperate": "PLIP",
}

TRPB_COLS = COMMON_COLS + [
    "chain_iptm_AB_avg_score_joint_chai",
    "chain_pae_min_BA_avg_score_joint_af3",
    "2:GLU-NH_2_avg_af3-struct_joint",
    "pocket-subcofcentroid-hw_avg - substrate-logp_af3-struct_joint",
    "num_hydrogen_bond_avg_af3-score_joint",
    # "num_interactions_avg_af3-score_joint"
]

TRPB_COL_DICT = {
    **COMMON_COL_DICT,
    "chain_iptm_AB_avg_score_joint_chai": "Chai-1",
    "chain_pae_min_BA_avg_score_joint_af3": "AF3",
    "2:GLU-NH_2_avg_af3-struct_joint": "Bond distance",
    "pocket-subcofcentroid-hw_avg - substrate-logp_af3-struct_joint": "Hydrophobicity",
    "num_hydrogen_bond_avg_af3-score_joint": "Hydrogen bonds",
    # "num_interactions_avg_af3-score_joint": "PLIP",
}

PARLQ_COLS = COMMON_HEME_COLS + ["0:C-C_1_avg_af3-struct_seperate"]

PARLQ_COL_DICT = {
    **COMMON_HEME_COL_dict,
    "0:C-C_1_avg_af3-struct_seperate": "Bond distance",
}

CB_COLS = COMMON_HEME_COLS + ["0:C-B_avg_af3-struct_seperate"]

CB_COL_DICT = {**COMMON_HEME_COL_dict, "0:C-B_avg_af3-struct_seperate": "Bond distance"}

CSI_COLS = COMMON_HEME_COLS + ["0:C-Si_avg_af3-struct_seperate"]
CSI_COL_DICT = {
    **COMMON_HEME_COL_dict,
    "0:C-Si_avg_af3-struct_seperate": "Bond distance",
}

FINAL_COL_ORDER = [
    "Library",
    "Hamming distance",
    "EVmutation",
    "ESM2",
    "ESM-IF",
    "CoVES",
    r"ΔΔ$G_f$",  # "ΔΔG",
    "Vina",
    "GALigandDock",
    "AF3",
    "Chai-1",
    "LigandMPNN",
    "FlowSite",
    "Bond distance",
    "Hydrogen bonds",
    "Hydrophobicity",
    "Active-site volume",
]

LIB_ORDER = [
    "PfTrpB-7iodo",
    "PfTrpB-7methyl",
    "PfTrpB-7bromo",
    "PfTrpB-5iodo",
    "PfTrpB-5bromo",
    "PfTrpB-5chloro",
    "PfTrpB-4bromo",
    "PfTrpB-6chloro",
    "PfTrpB-5cyano",
    "PfTrpB-4cyano",
    "PfTrpB-56chloro",
    "Rma-CB",
    "Rma-CSi",
    "ParLQ-a",
    "ParLQ-b",
    "ParLQ-c",
    "ParLQ-d",
    "ParLQ-e",
    "ParLQ-f",
    "ParLQ-g",
    "ParLQ-h",
    "ParLQ-i",
]


# Define metric categories
METRICS = ["rho", "ndcg", "ndcg10", "ndcg20", "ndcg25", "top10", "top20", "top25"]
METRICS_DICT = {
    "rho": "Spearman's ρ",
    "ndcg": "NDCG",
    "ndcg10": "NDCG@10%",
    "ndcg20": "NDCG@20%",
    "ndcg25": "NDCG@25%",
    "top10": "Top 10% recall",
    "top20": "Top 20% recall",
    "top25": "Top 25% recall",
}
FITSELE_DICT = {
    "fit": "activity",
    "sele": "selectivity",
}
FIT_METRICS = {metric: [] for metric in METRICS}
SELE_METRICS = {metric: [] for metric in METRICS}
ALL_METRICS_OPTS = [f"{t}_{metric}" for t in ["fit", "sele"] for metric in METRICS]

METRIC_COLOR_THRESHOLD = {
    "rho": [-0.25, 0.55],
    "ndcg": [0.6, 0.9],
    "ndcg10": [0.1, 0.75],
    "ndcg20": [0.15, 0.75],
    "ndcg25": [0.15, 0.75],
    "top10": [0, 0.6],
    "top20": [0, 0.6],
    "top25": [0, 0.6],
}


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
    y_true_minmax = (y_true - y_true.min()) / (
        y_true.max() - y_true.min()
    )  # Scale between 0 and 1
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


def fit_linear_model(X, y):
    """
    Case 1: Pure linear regression using sklearn.
    X: shape (m, n) - m data points, n features
    y: shape (m,)   - target values
    Returns:
        w0, w1, ..., wn as a dictionary
    """
    reg = LinearRegression(fit_intercept=True)
    reg.fit(X, y)
    intercept = reg.intercept_
    coefs = reg.coef_  # array of shape (n,)
    return {"intercept": intercept, "weights": coefs}


def inference_linear_model(X, w_0, w):
    y = np.dot(X, w) + w_0
    return y


def piecewise_transform(x, alpha1, alpha2):
    """
    Single-feature piecewise-linear transform phi(x; alpha1, alpha2).
    """
    if x < alpha1:
        return 0.0
    elif x >= alpha2:
        return 1.0
    else:
        # linear transition
        return (x - alpha1) / (alpha2 - alpha1)


def inference_piecewise_model(X, params):
    """
    Computes y_hat for all data points under the piecewise model.
    params structure:
        [w0, w1, ..., w_n, p1_1, p1_2, p2_1, p2_2, ..., pn_1, pn_2]
    where p_j_1 and p_j_2 define alpha_{j1} and alpha_{j2} via reparameterization:
        alpha_{j1} = params[ n+1 + 2*(j-1) ]
        alpha_{j2} = alpha_{j1} + exp( params[ n+1 + 2*(j-1) + 1] )
    """
    try:
        m, n = X.shape
    except ValueError:
        print("X is not 2D")

    # Parse out w0, w1..w_n
    w0 = params[0]
    w = params[1 : n + 1]  # shape (n,)

    # Parse out alpha params in reparameterized form
    alpha_list = []
    for j in range(n):
        base_index = n + 1 + 2 * j
        alpha1 = params[base_index]
        alpha2 = alpha1 + np.exp(params[base_index + 1])  # ensures alpha2 > alpha1
        alpha_list.append((alpha1, alpha2))

    y_hat = np.zeros(m)
    for i in range(m):
        # Start with intercept
        pred_i = w0
        # Add contribution from each feature
        for j in range(n):
            alpha1, alpha2 = alpha_list[j]
            pred_i += w[j] * piecewise_transform(X[i, j], alpha1, alpha2)
        y_hat[i] = pred_i

    return y_hat


def objective_piecewise(params, X, y):
    """
    Sum of squared errors for the piecewise model.
    """
    y_hat = inference_piecewise_model(X, params)
    return np.sum((y - y_hat) ** 2)


def fit_piecewise_model(X, y, max_iter=1000, random_seed=42):
    """
    Case 2: Piecewise-linear model with thresholds for ALL features.
    We do a nonlinear optimization over w0, w_j, alpha_{j1}, alpha_{j2}.

    X: shape (m, n)
    y: shape (m,)

    Returns:
        Dictionary of fitted parameters:
            {'intercept': w0, 'weights': array(w1..wn),
             'alphas': [(alpha11, alpha12), ..., (alpha_n1, alpha_n2)]}
    """
    np.random.seed(random_seed)
    m, n = X.shape

    # We'll do a random initialization for the parameters
    # Parameter layout: [w0, w1..w_n, p1_1, p1_2, p2_1, p2_2, ..., pn_1, pn_2]
    # length = (n+1) + 2*n = 3n + 1
    # reparameterize alpha2_j as alpha1_j + exp(...) to ensure alpha2_j > alpha1_j
    init = np.zeros(3 * n + 1)

    # Some random spread for w_j
    init[0] = np.mean(y)  # intercept guess
    init[1 : n + 1] = 0.01 * np.random.randn(n)

    # For alpha:
    # alpha1_j = random in [min(x_j), median(x_j)]
    # alpha2_j = alpha1_j + exp( something )
    for j in range(n):
        base_index = n + 1 + 2 * j
        col_j = X[:, j]
        alpha1_guess = np.percentile(col_j, 25)
        alpha2_guess = np.percentile(col_j, 75)
        if alpha2_guess <= alpha1_guess:
            alpha2_guess = alpha1_guess + 1.0  # fallback

        # store alpha1 directly
        init[base_index] = alpha1_guess
        # store log(alpha2 - alpha1)
        init[base_index + 1] = np.log(alpha2_guess - alpha1_guess)

    # Use 'L-BFGS-B' with no explicit bounds but can set them if needed
    res = minimize(
        fun=objective_piecewise,
        x0=init,
        args=(X, y),
        method="L-BFGS-B",
        options={"maxiter": max_iter},
    )

    if not res.success:
        print("Warning: Optimization did not converge successfully.")

    # Extract final parameters
    final_params = res.x
    w0 = final_params[0]
    w = final_params[1 : n + 1]
    alpha_pairs = []
    for j in range(n):
        base_index = n + 1 + 2 * j
        a1 = final_params[base_index]
        a2 = a1 + np.exp(final_params[base_index + 1])
        alpha_pairs.append((a1, a2))

    return {
        "intercept": w0,
        "weights": w,
        "alphas": alpha_pairs,
        "objective_value": res.fun,
    }


# --------------------------------------------------------
#  Case 3: Logistic (Sigmoid) for each feature
# --------------------------------------------------------
def logistic_transform(x, a, b):
    """
    Logistic function for a single feature:
    sigma(x; a, b) = 1 / (1 + exp(-b(x - a)))
    """
    return 1.0 / (1.0 + np.exp(-b * (x - a)))


def model_logistic(params, X):
    """
    Returns predicted y_hat for logistic transformations on each feature.
    params layout:
      [w0, w1..w_n, a1, b1, a2, b2, ..., an, bn]
    => length = (n+1) + 2n = 3n+1
    """
    m, n = X.shape
    w0 = params[0]
    w = params[1 : n + 1]

    ab_pairs = []
    for j in range(n):
        base_index = n + 1 + 2 * j
        a_j = params[base_index]
        b_j = params[base_index + 1]
        ab_pairs.append((a_j, b_j))

    y_hat = np.zeros(m)
    for i in range(m):
        pred_i = w0
        for j in range(n):
            a_j, b_j = ab_pairs[j]
            pred_i += w[j] * logistic_transform(X[i, j], a_j, b_j)
        y_hat[i] = pred_i
    return y_hat


def objective_logistic(params, X, y):
    y_hat = model_logistic(params, X)
    return np.sum((y - y_hat) ** 2)


def fit_logistic_model(X, y, max_iter=1000, random_seed=123):
    """
    Fits a model:
        y = w0 + sum_j w_j * logistic(x_j; a_j, b_j)
    for each feature j.
    """
    np.random.seed(random_seed)
    m, n = X.shape

    # total params = (n+1) + 2n = 3n+1
    init = np.zeros(3 * n + 1)

    # Initialize intercept
    init[0] = np.mean(y)

    # Initialize w_j
    init[1 : n + 1] = 0.01 * np.random.randn(n)

    # Initialize (a_j, b_j)
    for j in range(n):
        base_index = n + 1 + 2 * j
        median_j = np.median(X[:, j])
        init[base_index] = median_j  # location
        init[base_index + 1] = 0.01 * np.random.randn()  # small random slope

    res = minimize(
        fun=objective_logistic,
        x0=init,
        args=(X, y),
        method="L-BFGS-B",
        options={"maxiter": max_iter},
    )

    if not res.success:
        print("Warning: Logistic optimization did not converge.")

    final = res.x
    w0 = final[0]
    w = final[1 : n + 1]
    ab_pairs = []
    for j in range(n):
        base_index = n + 1 + 2 * j
        a_j = final[base_index]
        b_j = final[base_index + 1]
        ab_pairs.append((a_j, b_j))

    return {
        "intercept": w0,
        "weights": w,
        "ab_pairs": ab_pairs,
        "objective_value": res.fun,
    }


def predict_logistic(X_new, logistic_params):
    """
    Predict y_hat for new data X_new using the learned logistic model.
    logistic_params: dict with 'intercept', 'weights', 'ab_pairs'
    """
    w0 = logistic_params["intercept"]
    w = logistic_params["weights"]
    ab_pairs = logistic_params["ab_pairs"]

    m, n = X_new.shape
    y_hat = np.zeros(m)
    for i in range(m):
        pred_i = w0
        for j in range(n):
            a_j, b_j = ab_pairs[j]
            pred_i += w[j] * logistic_transform(X_new[i, j], a_j, b_j)
        y_hat[i] = pred_i
    return y_hat


# def fit_sigmoidal_model(X_train, y_train):
#     pass
#     return None


# def inference_sigmoida_model(X):
#     pass
#     return None


def generate_X_y(campaign, x_cols: list = [], y_name: str = None):
    df = pd.read_csv(campaign)
    X = df[x_cols].to_numpy()
    scaler = StandardScaler()
    X = scaler.fit_transform(X)

    y = df[y_name].to_numpy()

    assert X.shape[0] == y.shape[0], print(
        "ERROR: X and y have different amount of instances."
    )
    return X, y


def train_test_all(
    pattern="/disk2/fli/REVIVAL2/zs/comb/minimal/*.csv", output_dir="zs/lincomb"
):
    """ """
    pattern = [f for f in glob(pattern) if "_scope" not in f]

    # define
    y_name = "fitness"

    results_cols = FINAL_COL_ORDER[1:]

    if isinstance(pattern, str):
        lib_list = sorted(glob(pattern))
    else:
        lib_list = deepcopy(pattern)

    # extract ZS data from results csv
    campaigns = {}
    for campaign in lib_list:
        X, y = generate_X_y(campaign, results_cols, y_name)
        campaigns[get_file_name(campaign)] = (X, y)

        # handle NaNs
        if np.isnan(X).any():
            print(
                f"Error: found NaN values in X at {np.where(np.isnan(X))}. Setting them to zero"
            )
            X[np.where(np.isnan(X))] = 0
        elif np.isnan(y).any():
            print(f"Error: found NaN values in y at {np.where(np.isnan(y))}")
            y[np.where(np.isnan(y))] = 0

    # define validation folds
    lin_rho_pairwise = np.zeros(
        (len(campaigns.keys()), len(campaigns.keys())),
    )
    lin_piece_rho_pairwise = np.zeros(
        (len(campaigns.keys()), len(campaigns.keys())),
    )
    logistic_rho_pairwise = np.zeros(
        (len(campaigns.keys()), len(campaigns.keys())),
    )

    lin_params = {}
    lin_piece_params = {}
    logistic_params = {}

    for i, train_fold in tqdm(enumerate(campaigns.keys())):
        # define train fold
        X_train, y_train = campaigns[train_fold]

        # fit linear model
        params_lin = fit_linear_model(X_train, y_train)
        w_0_lin = params_lin["intercept"]
        w_lin = params_lin["weights"]

        # fit piecewise lin model
        params_piece = fit_piecewise_model(X_train, y_train)
        w_0_piece = params_piece["intercept"]
        w_piece = params_piece["weights"]
        alphas_tup_piece = params_piece["alphas"]
        alphas_piece = np.array([item for tup in alphas_tup_piece for item in tup])
        piecewise_params = np.concatenate(
            [np.array([w_0_piece]), w_piece, alphas_piece]
        )

        # fit logistic model
        logistic_params = fit_logistic_model(X_train, y_train)

        for j, test_fold in enumerate(campaigns.keys()):
            # define test_fold
            X_test, y_test = campaigns[test_fold]
            print(f"xtest{X_test[:10, :10]} ytest{y_test[:10]}")
            # lin model inference
            y_hat_lin = inference_linear_model(X_test, w_0_lin, w_lin)

            # lin model scoring
            lin_corr, _ = spearmanr(y_hat_lin, y_test)

            # piecewise lin model inference
            y_hat_piece = inference_piecewise_model(X=X_test, params=piecewise_params)

            # model scoring
            piece_lin_corr, _ = spearmanr(y_hat_piece, y_test)

            # logistic model scoring
            y_hat_logistic = predict_logistic(X_test, logistic_params)

            logistic_corr, _ = spearmanr(y_hat_logistic, y_test)

            # writing to corr arrays and paramweight json
            print(
                f"INFO: Fitted linear model and obbtained spearmanrho of {lin_corr} by fitting on:\n{train_fold}\nand scoring on:\n{test_fold}."
            )
            lin_rho_pairwise[i, j] = lin_corr
            lin_params[f"train_fold_{train_fold}_test_fold_{test_fold}"] = (
                w_0_lin,
                w_lin,
            )

            print(
                f"INFO: Fitted piecewise linear model and obbtained spearmanrho of {piece_lin_corr} by fitting on:\n{train_fold}\n and scoring on:\n{test_fold}."
            )
            lin_piece_rho_pairwise[i, j] = piece_lin_corr
            lin_piece_params[
                f"train_fold_{train_fold}_test_fold_{test_fold}"
            ] = piecewise_params

            logistic_rho_pairwise[i, j] = logistic_corr
            logistic_params[
                f"train_fold_{train_fold}_test_fold_{test_fold}"
            ] = logistic_params

    checkNgen_folder(output_dir)

    lin_rho_pairwise_df = pd.DataFrame(
        lin_rho_pairwise, columns=campaigns.keys(), index=campaigns.keys()
    )
    lin_rho_pairwise_df.to_csv(os.path.join(output_dir, "lin_rho_pairwise_df.csv"))
    np.savez(os.path.join(output_dir, "lin_params.npz"), **lin_params)

    lin_piece_rho_pairwise_df = pd.DataFrame(
        lin_piece_rho_pairwise, columns=campaigns.keys(), index=campaigns.keys()
    )
    lin_piece_rho_pairwise_df.to_csv(
        os.path.join(output_dir, "lin_piece_rho_pairwise_df.csv")
    )
    np.savez(os.path.join(output_dir, "lin_piece_params.npz"), **lin_piece_params)

    logistic_rho_pairwise_df = pd.DataFrame(
        logistic_rho_pairwise, columns=campaigns.keys(), index=campaigns.keys()
    )
    logistic_rho_pairwise_df.to_csv(
        os.path.join(output_dir, "logistic_rho_pairwise_df.csv")
    )
    np.savez(os.path.join(output_dir, "logistic_params.npz"), **logistic_params)


# clean up and save minimal comb
def clean_comb(in_path, lib):

    out_path = in_path.replace("comb", "comb/minimal")

    checkNgen_folder(out_path)

    df = pd.read_csv(in_path)

    # flip the sign
    for c in df.columns[1:]:
        if (
            "Triad_score" in c
            or "chain_pae_min" in c
            or "sum_" in c
            or "naive_score" in c
            or c == "complexscore"
            or c == "dH"
            or c == "ligscore"
            or c == "recscore"
            or c == "score"
            or c == "total_score"
            or "vina_" in c
            or "0:" in c
            or "1:" in c
            or "2:" in c
            or "_vol" in c
        ):
            df[c] = -df[c]
            print(f"flipped {c}")
    
    if "TrpB" in in_path:
        append_col = ["fitness"]
        df = df[TRPB_COLS[1:] + append_col].rename(columns=TRPB_COL_DICT)
    elif "ParLQ" in in_path:
        append_col = ["fitness", "selectivity"]
        df = df[PARLQ_COLS[1:] + append_col].rename(columns=PARLQ_COL_DICT)
    elif "Rma-CB" in in_path:
        append_col = ["fitness", "selectivity"]
        df = df[CB_COLS[1:] + append_col].rename(columns=CB_COL_DICT)
    elif "Rma-CSi" in in_path:
        append_col = ["fitness", "selectivity"]
        df = df[CSI_COLS[1:] + append_col].rename(columns=CSI_COL_DICT)


    return df[FINAL_COL_ORDER[1:]+append_col].to_csv(out_path, index=False)


def get_minimal_comb(pattern="/disk2/fli/REVIVAL2/zs/comb/*.csv"):
    """
    Reduce the comb to just the critical columns
    """
    for f in tqdm(sorted(glob(pattern))):
        if "_scope" in f:
            continue  # Skip files with "_scope"
        clean_comb(f, get_file_name(f))


def process_and_save_metrics(
    input_dir="/disk2/fli/REVIVAL2/zs/comb", output_dir="/disk2/fli/REVIVAL2/zs/metrics"
):
    """
    Processes CSV files to compute ranking metrics (NDCG, Spearman, Top-N recall)
    for fitness and selectivity and saves the results as CSV files.

    Parameters:
    - input_dir: str, path to the directory containing the CSV files.
    - output_dir: str, path to the directory where the results will be saved.
    """

    # Process each CSV file
    for f in sorted(glob(os.path.join(input_dir, "*.csv"))):
        if "_scope" in f:
            continue  # Skip files with "_scope"

        df_original = pd.read_csv(f)

        # Convert fitness column to numeric if it's a string
        if df_original["fitness"].dtype == "O":
            df_original["fitness"] = pd.to_numeric(
                df_original["fitness"].str.replace(",", ""), errors="coerce"
            )

        # Keep only numerical columns
        df_original = df_original[df_original.select_dtypes(include=["number"]).columns]

        df_name = os.path.basename(f).replace(".csv", "")  # Extract filename
        print(f"Processing: {df_name}")

        # Initialize metric dictionaries
        metric_dicts = {
            metric: dict.fromkeys(df_original.columns, np.nan) for metric in METRICS
        }
        for metric in metric_dicts.values():
            metric["lib"] = df_name

        add_selectivity = "selectivity" in df_original.columns
        sele_metric_dicts = (
            {metric: dict.fromkeys(df_original.columns, np.nan) for metric in METRICS}
            if add_selectivity
            else {}
        )

        for metric in sele_metric_dicts.values():
            metric["lib"] = df_name

        for c in df_original.columns:
            if c in ["fitness"] or "_rank" in c:
                continue

            # Remove NaN rows
            df = df_original[df_original[c].notna()]
            if df.empty:
                continue

            y_true = df["fitness"].values
            y_score = df[c].values

            # Flip sign based on conditions
            if (
                "Triad_score" in c
                or "chain_pae_min" in c
                or "sum_" in c
                or "naive_score" in c
                or c == "complexscore"
                or c == "dH"
                or c == "ligscore"
                or c == "recscore"
                or c == "score"
                or c == "total_score"
                or "vina_" in c
                or "0:" in c
                or "1:" in c
                or "2:" in c
                or "_vol" in c
            ):
                y_score = -1 * y_score
            elif " - " in c:
                y_score = np.abs(y_score)

            # Compute metrics
            metric_dicts["rho"][c] = (
                spearmanr(y_true, y_score).correlation
                if np.any(y_score != y_score[0])
                else np.nan
            )
            metric_dicts["ndcg"][c] = ndcg_scale(y_true=y_true, y_pred=y_score)
            metric_dicts["ndcg10"][c] = custom_ndcg(
                y_true=y_true, y_score=y_score, top=10
            )
            metric_dicts["ndcg20"][c] = custom_ndcg(
                y_true=y_true, y_score=y_score, top=20
            )
            metric_dicts["ndcg25"][c] = custom_ndcg(
                y_true=y_true, y_score=y_score, top=25
            )
            metric_dicts["top10"][c] = calc_top_n_percent_recall(
                y_true=y_true, y_score=y_score, top_n=10
            )
            metric_dicts["top20"][c] = calc_top_n_percent_recall(
                y_true=y_true, y_score=y_score, top_n=20
            )
            metric_dicts["top25"][c] = calc_top_n_percent_recall(
                y_true=y_true, y_score=y_score, top_n=25
            )

            # Selectivity calculations
            if add_selectivity and c != "selectivity":
                y_true_sele = df["selectivity"].values
                y_score_sele = y_score

                sele_metric_dicts["rho"][c] = (
                    spearmanr(y_true_sele, y_score_sele).correlation
                    if np.any(y_score_sele != y_score_sele[0])
                    else np.nan
                )
                sele_metric_dicts["ndcg"][c] = ndcg_scale(
                    y_true=y_true_sele, y_pred=y_score_sele
                )
                sele_metric_dicts["ndcg10"][c] = custom_ndcg(
                    y_true=y_true_sele, y_score=y_score_sele, top=10
                )
                sele_metric_dicts["ndcg20"][c] = custom_ndcg(
                    y_true=y_true_sele, y_score=y_score_sele, top=20
                )
                sele_metric_dicts["ndcg25"][c] = custom_ndcg(
                    y_true=y_true_sele, y_score=y_score_sele, top=25
                )
                sele_metric_dicts["top10"][c] = calc_top_n_percent_recall(
                    y_true=y_true_sele, y_score=y_score_sele, top_n=10
                )
                sele_metric_dicts["top20"][c] = calc_top_n_percent_recall(
                    y_true=y_true_sele, y_score=y_score_sele, top_n=20
                )
                sele_metric_dicts["top25"][c] = calc_top_n_percent_recall(
                    y_true=y_true_sele, y_score=y_score_sele, top_n=25
                )

        # Append results
        for metric in METRICS:
            FIT_METRICS[metric].append(metric_dicts[metric])
            if add_selectivity:
                SELE_METRICS[metric].append(sele_metric_dicts[metric])

    # Save results as DataFrames
    os.makedirs(output_dir, exist_ok=True)

    for metric, metric_list in FIT_METRICS.items():
        df = pd.DataFrame(metric_list)
        # rename parlq to parlq-a
        df["lib"] = df["lib"].apply(lambda x: "ParLQ-a" if x == "ParLQ" else x)
        df.to_csv(f"{output_dir}/fit_{metric}.csv", index=False)

    for metric, metric_list in SELE_METRICS.items():
        df = pd.DataFrame(metric_list)
        df["lib"] = df["lib"].apply(lambda x: "ParLQ-a" if x == "ParLQ" else x)
        df.to_csv(f"{output_dir}/sele_{metric}.csv", index=False)

    print("All metrics saved successfully!")


def process_metrics2plot(df, lib_order, group_name, col_dict, cols):
    """
    Filters, renames, and processes a subset of a DataFrame.

    Parameters:
    - df (pd.DataFrame): The input DataFrame.
    - lib_order (list): List of libraries to filter.
    - group_name (str): Name of the new average group (e.g., "PfTrpB-avg").
    - col_dict (dict): Mapping of column names for renaming.
    - cols (list): List of relevant columns to keep.

    Returns:
    - processed_df (pd.DataFrame): The processed DataFrame.
    - group_avg_df (pd.DataFrame): A DataFrame containing the average values for the group.
    """
    processed_df = df[df["lib"].isin(lib_order)].copy()
    print(processed_df.head())
    processed_df = processed_df[cols].rename(columns=col_dict).reset_index(drop=True)

    # Compute mean values
    group_avg_df = processed_df.set_index("Library").mean().to_frame().T
    group_avg_df["Library"] = group_name
    group_avg_df = group_avg_df[
        ["Library"] + group_avg_df.columns.tolist()[:-1]
    ].reset_index(drop=True)

    return processed_df, group_avg_df


def plot_metrics_heatmap(
    data,
    output_path,
    lib_order,
    colorbar_label,
    vmin,
    vmax,
    figsize=(8, 1.6),
    cbar_shrink=1.0,
):
    """
    Creates and saves a heatmap based on the processed data.

    Parameters:
    - data (pd.DataFrame): Data to plot.
    - output_path (str): Path to save the figure.
    - lib_order (list): Order for reindexing.
    - colorbar_label (str): Label for the colorbar.
    - figsize (tuple): Figure size.
    - cbar_shrink (float): Shrink factor for the colorbar.
    """
    plt.figure(figsize=figsize)
    sns.heatmap(
        data.set_index("Library").reindex(lib_order),
        cmap="vlag",
        annot=True,
        fmt=".1f",
        linewidths=0.5,
        cbar_kws={"label": colorbar_label, "shrink": cbar_shrink},
        vmin=vmin,
        vmax=vmax,
    )
    
    plt.ylabel("")
    plt.savefig(output_path, format="svg", dpi=300, bbox_inches="tight")
    print(f"Saved heatmap: {output_path}")


def processnplot_metrics(
    df,
    lib_order,
    trpb_cols,
    parlq_cols,
    cb_cols,
    csi_cols,
    trpb_col_dict,
    parlq_col_dict,
    cb_col_dict,
    csi_col_dict,
    final_col_order,
    metric_name,
    output_dir,
):
    """
    Processes data and generates heatmaps for a given metric.

    Parameters:
    - df (pd.DataFrame): The input DataFrame.
    - lib_order (list): Custom sorting order.
    - trpb_cols, parlq_cols, cb_cols, csi_cols (list): Relevant columns for each group.
    - trpb_col_dict, parlq_col_dict, cb_col_dict, csi_col_dict (dict): Column renaming mappings.
    - final_col_order (list): Final order of columns for the heatmap.
    - metric_name (str): Metric name (e.g., "fit_rho").
    - output_dir (str): Directory to save outputs.
    """
    # get colorbar_label
    fit2sele, metric = metric_name.split("_")

    # Process each group
    parlq_rho_df, parlqmean = process_metrics2plot(
        df, lib_order[-9:], "ParLQ-avg", parlq_col_dict, parlq_cols
    )

    cb_df = (
        df[df["lib"] == "Rma-CB"]
        .copy()[cb_cols]
        .rename(columns=cb_col_dict)
        .reset_index(drop=True)
    )
    csi_df = (
        df[df["lib"] == "Rma-CSi"]
        .copy()[csi_cols]
        .rename(columns=csi_col_dict)
        .reset_index(drop=True)
    )

    if fit2sele == "fit":
        trpb_rho_df, trpbmean = process_metrics2plot(
            df, lib_order[:11], "PfTrpB-avg", trpb_col_dict, trpb_cols
        )

        # Small summary heatmap
        summary_data = pd.concat([trpbmean, cb_df, csi_df, parlqmean])[final_col_order]
        summary_order = ["PfTrpB-avg", "Rma-CB", "Rma-CSi", "ParLQ-avg"]
        summary_size = (8, 1.6)

        # Full heatmap for all entries
        full_data = pd.concat([trpb_rho_df, cb_df, csi_df, parlq_rho_df])[
            final_col_order
        ]
        full_order = lib_order
        full_size = (7.2, 8)

    # no trpb
    else:
        # Small summary heatmap
        summary_data = pd.concat([cb_df, csi_df, parlqmean])[final_col_order]
        summary_order = ["Rma-CB", "Rma-CSi", "ParLQ-avg"]
        summary_size = (8, 1.24)

        # Full heatmap for all entries
        full_data = pd.concat([cb_df, csi_df, parlq_rho_df])[final_col_order]
        full_order = lib_order[-11:]
        full_size = (8, 4.2)

    colorbar_label = f"{METRICS_DICT[metric]}\n({FITSELE_DICT[fit2sele]})"
    vmin, vmax = METRIC_COLOR_THRESHOLD[metric]

    # Plot heatmaps
    plot_metrics_heatmap(
        data=summary_data,
        output_path=os.path.join(output_dir, f"zs_sum_{metric_name}.svg"),
        lib_order=summary_order,
        colorbar_label=colorbar_label,
        vmin=vmin,
        vmax=vmax,
        figsize=summary_size,
    )

    plot_metrics_heatmap(
        data=full_data,
        output_path=os.path.join(output_dir, f"zs_sum_{metric_name}_all.svg"),
        lib_order=full_order,
        colorbar_label=colorbar_label,
        vmin=vmin,
        vmax=vmax,
        figsize=full_size,
        cbar_shrink=0.5,
    )


def plot_all_metrics(input_dir="zs/metrics", output_dir="figs/metrics"):

    checkNgen_folder(output_dir)

    common_params = {
        "lib_order": LIB_ORDER,
        "trpb_cols": TRPB_COLS,
        "parlq_cols": PARLQ_COLS,
        "cb_cols": CB_COLS,
        "csi_cols": CSI_COLS,
        "trpb_col_dict": TRPB_COL_DICT,
        "parlq_col_dict": PARLQ_COL_DICT,
        "cb_col_dict": CB_COL_DICT,
        "csi_col_dict": CSI_COL_DICT,
        "final_col_order": FINAL_COL_ORDER,
    }

    for m in tqdm(ALL_METRICS_OPTS):
        fit2sele, metric = m.split("_")
        df = pd.read_csv(os.path.join(input_dir, f"{m}.csv"))

        # Process and plot for Spearman's ρ (activity)
        processnplot_metrics(df, **common_params, metric_name=m, output_dir=output_dir)