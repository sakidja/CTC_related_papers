import numpy as np
import rebound
import matplotlib.pyplot as plt
import os

# ============================================================
# USER CONTROLS
# ============================================================
G        = 1.0
N        = 2000
Rmax     = 15.0
R_inner  = 0.5
M_disk   = 0.2
m_part   = M_disk / N

R_scale  = 2.5
n_bins   = 50
min_bin_count = 30

softening = 0.1
dt        = 0.001
dt_output = 2.0
t_max     = 60.0

M_central = 4.0

# Spiral forcing controls (used only in Case B)
USE_SPIRAL = True
A_spiral   = 0.02
m_spiral   = 2
Omega_p    = 0.25

# Initial perturbations (seed only)
eps_seed_axisymmetric = 0.0
eps_seed_spiral_case  = 0.05

SAVE_PARTICLE_SNAPS = True
SNAPSHOT_STRIDE     = 2

OUTDIR = "crt_runs"
os.makedirs(OUTDIR, exist_ok=True)

# ============================================================
# Radial bins
# ============================================================
Rmin = R_inner
Rb = np.logspace(np.log10(Rmin), np.log10(Rmax), n_bins)  # bin centers

bin_edges = np.zeros(n_bins + 1)
bin_edges[1:-1] = 0.5 * (Rb[1:] + Rb[:-1])
bin_edges[0]    = Rmin
bin_edges[-1]   = Rmax

bin_areas = np.pi * (bin_edges[1:]**2 - bin_edges[:-1]**2)

# ============================================================
# Diagnostics helper
# ============================================================
def compute_profiles_from_state(x, y, vx, vy):
    R = np.sqrt(x*x + y*y)
    safe_R = np.where(R > 0, R, 1e-12)
    vR   = (x*vx + y*vy) / safe_R
    vphi = (x*vy - y*vx) / safe_R

    Vphi_prof   = np.full(n_bins, np.nan)
    vR_prof     = np.full(n_bins, np.nan)
    FL_prof     = np.full(n_bins, np.nan)
    Sigma_prof  = np.full(n_bins, np.nan)
    count_prof  = np.zeros(n_bins, dtype=int)

    for i in range(n_bins):
        mask = (R >= bin_edges[i]) & (R < bin_edges[i+1])
        n_here = np.sum(mask)
        count_prof[i] = n_here
        if n_here >= min_bin_count:
            Vphi_prof[i]  = np.mean(vphi[mask])                    # computed, not plotted
            vR_prof[i]    = np.mean(vR[mask])
            FL_prof[i]    = np.mean(R[mask] * vR[mask] * vphi[mask])
            Sigma_prof[i] = (n_here * m_part) / bin_areas[i]

    return Vphi_prof, vR_prof, FL_prof, Sigma_prof, count_prof

# ============================================================
# Spiral force callback
# ============================================================
def make_spiral_force(A, m, Omega_p):
    def spiral_force(reb_sim):
        s  = reb_sim.contents
        ps = s.particles
        t  = s.t

        for i in range(1, s.N):  # skip central mass
            p = ps[i]
            x, y = p.x, p.y
            R = np.hypot(x, y)
            if R < 1e-12:
                continue
            phi = np.arctan2(y, x)

            phase = m * (phi - Omega_p * t)
            aphi  = -A * np.sin(phase)

            # e_phi = (-sin phi, cos phi)
            p.ax += aphi * (-np.sin(phi))
            p.ay += aphi * ( np.cos(phi))
    return spiral_force

# ============================================================
# Snapshot plot helper (xy only)
# ============================================================
def save_xy_snapshot_plot(pref, t_snap, x, y):
    plt.figure(figsize=(6,6))
    plt.scatter(x[1:], y[1:], s=1)               # skip central mass
    plt.scatter([x[0]], [y[0]], s=40, marker="x")
    plt.gca().set_aspect("equal", "box")
    plt.xlim(-Rmax, Rmax)
    plt.ylim(-Rmax, Rmax)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(f"{os.path.basename(pref)}: particle snapshot (t={t_snap:.1f})")
    plt.tight_layout()
    plt.savefig(pref + f"_xy_t{t_snap:.1f}.png", dpi=200)
    plt.close()

# ============================================================
# Build simulation
# ============================================================
def build_sim(eps_seed, use_spiral, case_prefix):
    sim = rebound.Simulation()
    sim.G = G
    sim.softening = softening
    sim.integrator = "whfast"
    sim.dt = dt

    # Central mass
    sim.add(m=M_central, x=0., y=0., z=0., vx=0., vy=0., vz=0.)

    rng = np.random.default_rng(1234)

    for _ in range(N):
        while True:
            R = rng.exponential(scale=R_scale)
            if R_inner <= R <= Rmax:
                break
        phi = 2.0 * np.pi * rng.random()

        x = R * np.cos(phi)
        y = R * np.sin(phi)

        M_disk_enclosed = M_disk * (1.0 - np.exp(-R/R_scale)*(1.0 + R/R_scale))
        M_enclosed = M_central + M_disk_enclosed

        v_circ = np.sqrt(G * M_enclosed / R)

        # Seed perturbation in vphi (optional)
        vphi = v_circ * (1.0 + eps_seed * np.cos(2.0 * phi))

        vx = -vphi * np.sin(phi)
        vy =  vphi * np.cos(phi)

        sim.add(m=m_part, x=x, y=y, z=0.0, vx=vx, vy=vy, vz=0.0)

    sim.move_to_com()

    if use_spiral:
        sim.additional_forces = make_spiral_force(A_spiral, m_spiral, Omega_p)

    print(f"[{case_prefix}] Built system. sim.N={sim.N} | eps_seed={eps_seed} | spiral={use_spiral}")
    return sim

# ============================================================
# Run one case
# ============================================================
def run_case(case_prefix, eps_seed, use_spiral):
    sim = build_sim(eps_seed, use_spiral, case_prefix)
    pref = os.path.join(OUTDIR, case_prefix)

    # Initial state
    parts = sim.particles
    x0  = np.array([p.x  for p in parts])
    y0  = np.array([p.y  for p in parts])
    vx0 = np.array([p.vx for p in parts])
    vy0 = np.array([p.vy for p in parts])

    Vphi_initial, vR_initial, FL_initial, Sigma_initial, _ = compute_profiles_from_state(x0, y0, vx0, vy0)

    times, Vphi_history, FL_history, vR_history, Sigma_history = [], [], [], [], []
    particle_snaps = []

    times.append(0.0)
    Vphi_history.append(Vphi_initial)
    FL_history.append(FL_initial)
    vR_history.append(vR_initial)
    Sigma_history.append(Sigma_initial)

    if SAVE_PARTICLE_SNAPS:
        particle_snaps.append({"t": 0.0, "x": x0.copy(), "y": y0.copy(), "vx": vx0.copy(), "vy": vy0.copy()})
        save_xy_snapshot_plot(pref, 0.0, x0, y0)

    t = 0.0
    step_index = 0
    while t < t_max:
        t_next = min(t + dt_output, t_max)
        sim.integrate(t_next)
        t = t_next

        parts = sim.particles
        x  = np.array([p.x  for p in parts])
        y  = np.array([p.y  for p in parts])
        vx = np.array([p.vx for p in parts])
        vy = np.array([p.vy for p in parts])

        Vphi_prof, vR_prof, FL_prof, Sigma_prof, _ = compute_profiles_from_state(x, y, vx, vy)

        times.append(t)
        Vphi_history.append(Vphi_prof)
        FL_history.append(FL_prof)
        vR_history.append(vR_prof)
        Sigma_history.append(Sigma_prof)

        if SAVE_PARTICLE_SNAPS and (step_index % SNAPSHOT_STRIDE == 0):
            particle_snaps.append({"t": t, "x": x.copy(), "y": y.copy(), "vx": vx.copy(), "vy": vy.copy()})
            save_xy_snapshot_plot(pref, t, x, y)

        step_index += 1
        print(f"[{case_prefix}] Integrated to t={t:.1f}")

    times         = np.array(times)
    Vphi_history  = np.array(Vphi_history)
    FL_history    = np.array(FL_history)
    vR_history    = np.array(vR_history)
    Sigma_history = np.array(Sigma_history)

    Vphi_final = Vphi_history[-1]
    FL_final   = FL_history[-1]

    # Surplus + Curl (computed; no V(R) plotting)
    Erot_final = 0.5 * (Rb**2) * (Vphi_final**2)

    mask_outer = (Rb > 0.5 * Rmax) & np.isfinite(Erot_final)
    E_outer = np.nansum(Erot_final[mask_outer])
    E_total = np.nansum(Erot_final[np.isfinite(Erot_final)])
    Surplus_sim = E_outer / (E_total + 1e-12)

    RFL = Rb * FL_final
    finite = np.isfinite(RFL)
    if np.sum(finite) > 1:
        dRFL_dR = np.gradient(RFL[finite], Rb[finite])
        Curl_sim = np.nanmean(np.abs(dRFL_dR))
    else:
        Curl_sim = np.nan

    print(f"[{case_prefix}] Surplus Index: {Surplus_sim}")
    print(f"[{case_prefix}] Curl Index  : {Curl_sim}")

    # ----------------- SAVE PLOTS -----------------

    # FL profile 
    plt.figure(figsize=(7,5))
    plt.axhline(0.0, color="gray", lw=0.8)
    plt.plot(Rb, FL_initial, "--k", label="Initial $F_L$")
    plt.plot(Rb, FL_final, "r", label=f"Final $F_L$ (t={t_max:.1f})")
    plt.xlabel("Radius R")
    plt.ylabel(r"$F_L(R)=\langle R v_R v_\phi\rangle$")
    plt.title(f"{case_prefix}: angular-momentum flux profile")
    plt.legend()
    plt.tight_layout()
    plt.savefig(pref + "_FL_profile.png", dpi=200)
    plt.close()

    # FL heatmap
    plt.figure(figsize=(7,5))
    FL_plot = np.where(np.isfinite(FL_history), FL_history, 0.0)
    plt.imshow(FL_plot, aspect="auto", origin="lower",
               extent=[Rb.min(), Rb.max(), times[0], times[-1]],
               cmap="RdBu_r")
    plt.colorbar(label=r"$F_L(R,t)$")
    plt.xlabel("Radius R")
    plt.ylabel("Time")
    plt.title(f"{case_prefix}: $F_L(R,t)$")
    plt.tight_layout()
    plt.savefig(pref + "_FL_heatmap.png", dpi=200)
    plt.close()

    # vR heatmap
    plt.figure(figsize=(7,5))
    vR_plot = np.where(np.isfinite(vR_history), vR_history, 0.0)
    plt.imshow(vR_plot, aspect="auto", origin="lower",
               extent=[Rb.min(), Rb.max(), times[0], times[-1]],
               cmap="RdBu_r")
    plt.colorbar(label=r"$\langle v_R\rangle(R,t)$")
    plt.xlabel("Radius R")
    plt.ylabel("Time")
    plt.title(f"{case_prefix}: radial streaming")
    plt.tight_layout()
    plt.savefig(pref + "_vR_heatmap.png", dpi=200)
    plt.close()

    # Sigma heatmap
    plt.figure(figsize=(7,5))
    S_plot = np.where(np.isfinite(Sigma_history), Sigma_history, 0.0)
    plt.imshow(S_plot, aspect="auto", origin="lower",
               extent=[Rb.min(), Rb.max(), times[0], times[-1]],
               cmap="viridis")
    plt.colorbar(label=r"$\Sigma(R,t)$")
    plt.xlabel("Radius R")
    plt.ylabel("Time")
    plt.title(f"{case_prefix}: surface density")
    plt.tight_layout()
    plt.savefig(pref + "_Sigma_heatmap.png", dpi=200)
    plt.close()

    # Surplus vs Curl
    plt.figure(figsize=(5,4))
    plt.scatter(Surplus_sim, Curl_sim, s=80)
    plt.xlabel("Surplus Index")
    plt.ylabel("Curl Index")
    plt.title(f"{case_prefix}: Surplus vs Curl")
    plt.tight_layout()
    plt.savefig(pref + "_surplus_vs_curl.png", dpi=200)
    plt.close()

    # ----------------- SAVE SNAPS -----------------
    if SAVE_PARTICLE_SNAPS:
        print(f"[{case_prefix}] Saving particle snapshots...")
        for snap in particle_snaps:
            t_snap = snap["t"]
            dat_name = pref + f"_particles_t{t_snap:.1f}.dat"
            npz_name = pref + f"_particles_t{t_snap:.1f}.npz"
            data = np.column_stack([snap["x"], snap["y"], snap["vx"], snap["vy"]])
            np.savetxt(dat_name, data, header=f"t={t_snap:.1f}\n x y vx vy", fmt="%.6f")
            np.savez(npz_name, t=t_snap, x=snap["x"], y=snap["y"], vx=snap["vx"], vy=snap["vy"])

    return Surplus_sim, Curl_sim

# ============================================================
# MAIN
# ============================================================
if __name__ == "__main__":
    S0, C0 = run_case("caseA_eps0", eps_seed_axisymmetric, use_spiral=False)
    S1, C1 = run_case("caseB_spiral", eps_seed_spiral_case, use_spiral=USE_SPIRAL)

    print("\n================= SUMMARY =================")
    print(f"Case A (eps=0, no spiral):   Surplus={S0:.6f}  Curl={C0:.6f}")
    print(f"Case B (spiral forcing on):  Surplus={S1:.6f}  Curl={C1:.6f}")
    print(f"Outputs saved in: {OUTDIR}/")

