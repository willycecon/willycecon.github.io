"""
get_duh.py
----------
Compute the DUH (Dominance-Unevenness-Heterogeneity) index and related
heterogeneity measures for compositional data.

Dependencies
------------
    pip install pandas numpy mpmath

Usage
-----
    from get_duh import get_duh
    result = get_duh(df, groups="g_", evepar=2)
"""

import numpy as np
import pandas as pd
from mpmath import mp, mpf, log, fabs


def get_duh(
    data: pd.DataFrame,
    groups,
    evepar: float = 2.0,
    round_digits: int = 4,
    show_groups: bool = True,
    prefix: str = None,
    precision: int = 50,
) -> pd.DataFrame:
    """
    Compute DUH and related heterogeneity measures.

    Parameters
    ----------
    data : pd.DataFrame
        Input data. Each row is an observation.
    groups : str or list of str
        If str, used as a column-name prefix to identify group proportion columns.
        If list, used directly as column names.
        Proportions should sum to 1 per row; if not, automatic normalisation
        is applied.
    evepar : float
        Evenness parameter. Must be strictly > 1. Default 2.
    round_digits : int
        Decimal places for output columns. Default 4.
    show_groups : bool
        Print the identified group column names. Default True.
    prefix : str or None
        Optional string prepended to the output column names DUH, GSI, SE.
    precision : int
        Decimal-digit precision for mpmath high-precision arithmetic used
        in the sigma_tilde / evenness calculation. Default 50.

    Returns
    -------
    pd.DataFrame
        Original data with five columns appended:
        ``sigma_1``, ``evenness``, ``[prefix]DUH``, ``[prefix]GSI``, ``[prefix]SE``.

    Notes
    -----
    High-precision arithmetic (mpmath) is applied exclusively to the
    ``sigma_tilde`` and evenness (``psi``) calculation to avoid catastrophic
    cancellation when minor-group proportions are near-equal.

    Examples
    --------
    >>> import pandas as pd
    >>> from get_duh import get_duh
    >>> df = pd.DataFrame({"g_a": [0.5, 0.8], "g_b": [0.3, 0.1], "g_c": [0.2, 0.1]})
    >>> get_duh(df, groups="g_", evepar=2)
    """
    if evepar <= 1:
        raise ValueError("evepar must be strictly greater than 1 for convexity")

    mp.dps = precision  # set mpmath decimal-place precision

    # ── Identify group columns ──────────────────────────────────────────────
    if isinstance(groups, list):
        g2    = groups
        abs_s = len(groups)
    else:
        g2    = [c for c in data.columns if c.startswith(groups)]
        abs_s = len(g2)

    if abs_s < 2:
        raise ValueError(
            f"At least 2 group columns are required. "
            f"Found {abs_s} column(s) matching '{groups}'."
        )

    subdata = data[g2].copy().astype(float)

    # ── Check / normalise proportions ───────────────────────────────────────
    row_sums = subdata.sum(axis=1)
    if not np.allclose(row_sums, 1.0, atol=1e-6):
        print("Total Group Proportions Do Not Add Up To 1, automatic scaling applied")
        subdata = subdata.div(row_sums, axis=0)
    else:
        print("Check: Group Proportions Add Up to 1")

    # ── High-precision evenness (psi) for a single observation ─────────────
    def compute_psi(minor_p_vec: np.ndarray) -> float:
        """
        Compute evenness using mpmath to avoid catastrophic cancellation.

        sigma_tilde_i = log(evepar) * p_i / sum(minor_p)
        psi_terms_i   = |sigma_tilde_i - log(evepar) / (abs_s - 1)|^evepar
        psi           = 1 - sum(psi_terms)^(1/evepar) / log(evepar)
        """
        mp_minor_sum = sum(mpf(x) for x in minor_p_vec)
        if mp_minor_sum == 0:
            return 0.0  # edge case: all weight on sigma_1

        mp_evepar   = mpf(evepar)
        mp_factor   = log(mp_evepar)
        ideal       = mp_factor / (abs_s - 1)

        sigma_tilde = [mp_factor * (mpf(p) / mp_minor_sum) for p in minor_p_vec]
        psi_terms   = [fabs(s - ideal) ** mp_evepar for s in sigma_tilde]
        psi_val     = float(1 - sum(psi_terms) ** (mpf(1) / mp_evepar) / mp_factor)
        return psi_val

    # ── Row-wise calculation ────────────────────────────────────────────────
    records = []
    for _, row in subdata.iterrows():
        p_vals   = row.to_numpy()
        sorted_p = np.sort(p_vals)[::-1]
        sigma_1  = float(sorted_p[0])
        minor_p  = sorted_p[1:]

        # Evenness (high precision)
        psi = compute_psi(minor_p)

        # DUH: log(sigma_1) / log(abs_s) scaling
        # When sigma_1 == 1 (monopoly), log(sigma_1) == 0 → DUH = 0
        if sigma_1 <= 0.0 or sigma_1 >= 1.0:
            duh = 0.0
        else:
            duh = -(np.log(sigma_1) / np.log(abs_s)) * psi
            if np.isnan(duh):
                duh = 0.0

        # HHI / Gini-Simpson Index
        hhi = float(np.sum(p_vals ** 2))
        gsi = 1.0 - hhi

        # Shannon Entropy
        with np.errstate(divide="ignore", invalid="ignore"):
            plp = np.where(p_vals > 0, -p_vals * np.log(p_vals), 0.0)
        se = float(np.sum(plp))

        records.append(
            {
                "sigma_1":  round(sigma_1, round_digits),
                "evenness": round(psi,     round_digits),
                "_DUH":     round(duh,     round_digits),
                "_GSI":     round(gsi,     round_digits),
                "_SE":      round(se,      round_digits),
            }
        )

    results = pd.DataFrame(records, index=subdata.index)

    # ── Assemble output ─────────────────────────────────────────────────────
    duh_col = f"{prefix}DUH" if prefix else "DUH"
    gsi_col = f"{prefix}GSI" if prefix else "GSI"
    se_col  = f"{prefix}SE"  if prefix else "SE"

    out = data.copy()
    out["sigma_1"]  = results["sigma_1"]
    out["evenness"] = results["evenness"]
    out[duh_col]    = results["_DUH"]
    out[gsi_col]    = results["_GSI"]
    out[se_col]     = results["_SE"]

    print(f"Universe has {abs_s} categories")
    if show_groups:
        print(", ".join(g2))
    print(f"Using {evepar}-metric")

    return out
