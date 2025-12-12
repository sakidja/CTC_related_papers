# ============================================================
# SPARC rotation analysis: run on ALL galaxies in Rotmod_LTG
# - computes indices for each galaxy
# - saves one PNG figure per galaxy
# - writes a text summary (one line per galaxy)
#   NOW ALSO INCLUDES:
#   - R_max_kpc          (outermost measured radius)
#   - R_outer_start_kpc  (0.5 * R_max, "outer half" for surplus)
#   - R_core_kpc         (end of inner rigid / steep-rise zone)
#   - R_turnover_kpc     (radius of max V_obs)
# ============================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# Folder where you extracted Rotmod_LTG.zip
ROT_DIR = Path("Rotmod_LTG")


# ------------------------------------------------------------
# 1. Loader
# ------------------------------------------------------------
def load_sparc_rc(galaxy_name: str, rotmod_dir: Path = ROT_DIR) -> pd.DataFrame:
    """
    Load a SPARC rotation curve file from the extracted Rotmod_LTG archive.

    It finds files whose names start with the given galaxy_name, e.g.:
      NGC2403_rotmod.dat, DDO161_rotmod.dat, etc.
    """

    rotmod_dir = Path(rotmod_dir)

    candidates = sorted(rotmod_dir.glob(f"{galaxy_name}*"))
    if not candidates:
        raise FileNotFoundError(
            f"No file starting with '{galaxy_name}' found in {rotmod_dir}."
        )

    path = candidates[0]
    print(f"✔ Using file: {path.name}")

    df = pd.read_csv(
        path,
        comment="#",
        delim_whitespace=True,
        header=None,
        engine="python",
    )

    print(f"  Raw shape: {df.shape}")

    # Allow up to 10 columns safely
    colnames = [
        "R", "V_obs", "e_V_obs",
        "V_gas", "V_disk", "V_bulge",
        "V_halo", "V_tot", "col9", "col10",
    ]
    df.columns = colnames[: df.shape[1]]

    # Ensure baryonic components exist
    for cname in ["V_gas", "V_disk", "V_bulge"]:
        if cname not in df.columns:
            df[cname] = 0.0

    # Keep only positive radii
    df = df[df["R"] > 0].reset_index(drop=True)
    return df


# ------------------------------------------------------------
# 2. Metrics
# ------------------------------------------------------------
def compute_baryonic_velocity(df):
    V_bar2 = df["V_gas"] ** 2 + df["V_disk"] ** 2 + df["V_bulge"] ** 2
    return np.sqrt(np.clip(V_bar2, 0, None))


def compute_fractional_uplift(df, V_bar, min_vbar=1e-3):
    V_obs2 = df["V_obs"].values ** 2
    V_bar2 = V_bar ** 2

    U = np.full_like(V_bar, np.nan, dtype=float)
    mask = V_bar > min_vbar
    U[mask] = (V_obs2[mask] - V_bar2[mask]) / V_bar2[mask]
    return U


def compute_curl_index(df, V_bar):
    R = df["R"].values
    V_obs2 = df["V_obs"].values ** 2
    V_bar2 = V_bar ** 2

    tau_obs = R * V_obs2
    tau_bar = R * V_bar2

    dR = np.diff(R)
    dtau_obs = np.diff(tau_obs) / dR
    dtau_bar = np.diff(tau_bar) / dR

    num = np.trapz(np.abs(dtau_obs - dtau_bar), R[1:])
    den = np.trapz(np.abs(dtau_bar),            R[1:])
    if den <= 0:
        return np.nan
    return num / den


def compute_surplus_index(
    df,
    V_bar,
    inner_fraction: float = 0.5,
    return_radii: bool = False,
):
    """
    SurplusIndex as before, but can optionally return:
      R_max, R_outer_start = inner_fraction * R_max
    so you can compare with your model radii.
    """
    R = df["R"].values
    V_obs2 = df["V_obs"].values ** 2
    V_bar2 = V_bar ** 2

    E_extra = np.clip(V_obs2 - V_bar2, 0, None)

    R_max = R.max()
    r0 = inner_fraction * R_max
    mask_outer = R >= r0

    num = np.trapz(E_extra[mask_outer], R[mask_outer])
    den = np.trapz(V_obs2, R)
    if den <= 0:
        S = np.nan
    else:
        S = num / den

    if return_radii:
        return S, R_max, r0
    else:
        return S


def compute_structure_radii(df):
    """
    Compute:
      - R_turnover: radius where V_obs is maximal
      - R_core: radius where the rising inner slope dV/dR
                first drops below a fraction of its initial value.
    Returns (R_core, R_turnover).
    """
    R = df["R"].values
    V = df["V_obs"].values

    if len(R) < 3:
        return np.nan, np.nan

    # --- R_turnover: radius of max V_obs ---
    i_max = int(np.nanargmax(V))
    R_turnover = float(R[i_max])

    # --- R_core: slope-based definition ---
    # finite-difference slope on midpoints
    dR = np.diff(R)
    dV = np.diff(V)
    slope = dV / dR          # length N-1
    R_mid = 0.5 * (R[1:] + R[:-1])

    # estimate initial slope using first few points with positive slope
    mask_pos = slope > 0
    if not np.any(mask_pos):
        return np.nan, R_turnover

    # use median of first up-to-3 positive slopes as "inner slope"
    positive_indices = np.where(mask_pos)[0]
    first_pos = positive_indices[:3]
    s0 = np.median(slope[first_pos])

    if s0 <= 0:
        return np.nan, R_turnover

    frac_drop = 0.3  # core ends when slope < 0.3 * initial slope
    drop_mask = slope < frac_drop * s0

    # R_core = first radius where slope has dropped enough
    R_core = np.nan
    if np.any(drop_mask):
        i_core = np.where(drop_mask)[0][0]
        R_core = float(R_mid[i_core])
    else:
        # if never drops, leave as NaN
        R_core = np.nan

    return R_core, R_turnover


# ------------------------------------------------------------
# 3. Single-galaxy analysis (with optional figure saving)
# ------------------------------------------------------------
def analyze_galaxy(
    galaxy_name: str,
    rotmod_dir: Path = ROT_DIR,
    show_plots: bool = False,
    savefig_path: Path = None,
):
    df = load_sparc_rc(galaxy_name, rotmod_dir)

    R = df["R"].values
    V_obs = df["V_obs"].values
    V_bar = compute_baryonic_velocity(df)
    U = compute_fractional_uplift(df, V_bar)

    CI = compute_curl_index(df, V_bar)

    # SurplusIndex plus radii
    S, R_max, R_outer_start = compute_surplus_index(
        df,
        V_bar,
        inner_fraction=0.5,
        return_radii=True,
    )

    # Structural radii
    R_core, R_turnover = compute_structure_radii(df)

    print(
        f"  CurlIndex = {CI:.4f}, "
        f"SurplusIndex = {S:.4f}, "
        f"mean U = {np.nanmean(U):.4f}"
    )
    print(
        f"  R_max = {R_max:.3f} kpc,  "
        f"R_outer_start (0.5 R_max) = {R_outer_start:.3f} kpc"
    )
    print(
        f"  R_core ≈ {R_core:.3f} kpc,  "
        f"R_turnover ≈ {R_turnover:.3f} kpc"
    )

    # Make plots only if requested
    if show_plots or (savefig_path is not None):
        fig, axes = plt.subplots(3, 1, figsize=(6, 10), sharex=True)

        # 1. Rotation curve
        ax = axes[0]
        ax.plot(R, V_obs, "o-", label="V_obs")
        ax.plot(R, V_bar, "s--", label="V_bar (baryons)")
        if not np.isnan(R_core):
            ax.axvline(R_core, color="orange", linestyle=":", label="R_core")
        ax.axvline(R_turnover, color="red", linestyle="--", label="R_turnover")
        ax.set_ylabel("V (km/s)")
        ax.set_title(f"{galaxy_name}: Rotation Curve")
        ax.legend()
        ax.grid(alpha=0.3)

        # 2. Uplift
        ax = axes[1]
        ax.plot(R, U, "o-")
        ax.axhline(0, color="k", linewidth=0.8)
        ax.set_ylabel("U(r)")
        ax.set_title("Fractional Uplift U(r)")
        ax.grid(alpha=0.3)

        # 3. Extra rotational energy (V_obs^2 - V_bar^2)
        V_obs2 = V_obs ** 2
        V_bar2 = V_bar ** 2
        E_extra = np.clip(V_obs2 - V_bar2, 0, None)

        ax = axes[2]
        ax.plot(R, E_extra, "o-")
        ax.axvline(R_outer_start, color="green", linestyle="--",
                   label="outer start (0.5 R_max)")
        ax.set_xlabel("R (kpc)")
        ax.set_ylabel("Extra V² (km²/s²)")
        ax.set_title(f"Extra Rotational Energy (Surplus Index ≈ {S:.3f})")
        ax.grid(alpha=0.3)
        ax.legend()

        plt.tight_layout()

        if savefig_path is not None:
            savefig_path = Path(savefig_path)
            savefig_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(savefig_path, dpi=150, bbox_inches="tight")
            print(f"  Saved figure to {savefig_path}")

        if show_plots:
            plt.show()
        else:
            plt.close(fig)

    return {
        "galaxy": galaxy_name,
        "CurlIndex": float(CI),
        "SurplusIndex": float(S),
        "MeanUplift": float(np.nanmean(U)),
        "R_max": float(R_max),
        "R_outer_start": float(R_outer_start),
        "R_core": float(R_core) if not np.isnan(R_core) else np.nan,
        "R_turnover": float(R_turnover),
    }


# ------------------------------------------------------------
# 4. Batch processing: run on ALL galaxies
# ------------------------------------------------------------
def process_all_galaxies(
    rotmod_dir: Path = ROT_DIR,
    out_fig_dir: Path = Path("sparc_figs"),
    summary_txt: Path = Path("sparc_indices.txt"),
):
    rotmod_dir = Path(rotmod_dir)
    out_fig_dir = Path(out_fig_dir)
    out_fig_dir.mkdir(parents=True, exist_ok=True)

    results = []

    # We will scan all .dat files and infer galaxy names from the prefix
    dat_files = sorted(rotmod_dir.glob("*.dat"))
    print(f"\nFound {len(dat_files)} .dat files in {rotmod_dir}")

    for path in dat_files:
        # Example filename: NGC2403_rotmod.dat -> galaxy_name = "NGC2403"
        stem = path.stem   # "NGC2403_rotmod"
        galaxy_name = stem.split("_")[0]

        print(f"\n=== Analyzing {galaxy_name} ===")
        try:
            res = analyze_galaxy(
                galaxy_name,
                rotmod_dir=rotmod_dir,
                show_plots=False,
                savefig_path=out_fig_dir / f"{galaxy_name}.png",
            )
            results.append(res)
        except Exception as e:
            print(f"  !! ERROR for {galaxy_name}: {e}")

    # Write a simple text summary (one line per galaxy)
    summary_txt = Path(summary_txt)
    with summary_txt.open("w") as f:
        f.write(
            "# galaxy CurlIndex SurplusIndex MeanUplift "
            "R_max_kpc R_outer_start_kpc R_core_kpc R_turnover_kpc\n"
        )
        for r in results:
            f.write(
                f"{r['galaxy']}  "
                f"{r['CurlIndex']:.6g}  "
                f"{r['SurplusIndex']:.6g}  "
                f"{r['MeanUplift']:.6g}  "
                f"{r['R_max']:.6g}  "
                f"{r['R_outer_start']:.6g}  "
                f"{r['R_core']:.6g}  "
                f"{r['R_turnover']:.6g}\n"
            )

    print(f"\n✔ Wrote index summary to {summary_txt}")
    print(f"✔ Saved {len(results)} figures to {out_fig_dir}")


# ------------------------------------------------------------
# 5. Main
# ------------------------------------------------------------
if __name__ == "__main__":
    process_all_galaxies(
        rotmod_dir=ROT_DIR,
        out_fig_dir=Path("sparc_figs"),
        summary_txt=Path("sparc_indices.txt"),
    )

