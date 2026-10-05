# Vectorized Fully-Coupled Newton-Raphson Solver for Aestimo 1D Drift-Diffusion
# Replaces legacy Krylov solver with an exact analytic block-tridiagonal direct solver
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import logging
import warnings

logger = logging.getLogger(__name__)


class NewtonConvergenceError(RuntimeError):
    """A failed voltage step cannot be exported as a converged device result."""

    def __init__(self, voltage, diagnostics):
        self.voltage = float(voltage)
        self.diagnostics = dict(diagnostics)
        super().__init__(
            "Mode 10 failed at Va={:.6g} V: {} (iterations={}, residual={:.6g}).".format(
                self.voltage, diagnostics.get('reason', 'not converged'),
                diagnostics.get('iterations', 0), diagnostics.get('residual_norm', float('inf'))
            )
        )

def Ber(x):
    """Numerically stable Scharfetter-Gummel Bernoulli function: B(x) = x / (exp(x) - 1)."""
    ax = np.abs(x)
    res = np.empty_like(x, dtype=float)
    small = ax < 1e-4
    res[small] = 1.0 - 0.5 * x[small] + (x[small]**2) / 12.0
    pos = x > 80.0
    res[pos] = 0.0
    neg = x < -80.0
    res[neg] = -x[neg]
    norm = (~small) & (~pos) & (~neg)
    xn = x[norm]
    res[norm] = xn / np.expm1(xn)
    return res

def Ber_prime(x):
    """Derivative of the Bernoulli function B'(x) = (exp(x) - 1 - x*exp(x)) / (exp(x) - 1)^2."""
    ax = np.abs(x)
    res = np.empty_like(x, dtype=float)
    small = ax < 1e-4
    res[small] = -0.5 + x[small] / 6.0
    pos = x > 80.0
    res[pos] = 0.0
    neg = x < -80.0
    res[neg] = -1.0
    norm = (~small) & (~pos) & (~neg)
    xn = x[norm]
    ex = np.exp(xn)
    exm1 = ex - 1.0
    res[norm] = (exm1 - xn * ex) / (exm1**2)
    return res

class CoupledNewtonSolver:
    """
    Fully-Coupled Newton-Raphson Drift-Diffusion Solver for 1D Semiconductor Devices.
    Simultaneously solves non-linear Poisson, electron continuity, and hole continuity
    equations on a 1D grid with exact analytical Jacobian and sparse direct factorization.
    """
    def __init__(self, model, n_max, dx, ni, dop, Ldi, Ppz_Psp, pol_surf_char, Nc, Nv, fi_stat, n_stat=None, p_stat=None):
        self.model = model
        self.n_max = int(n_max)
        self.dx = float(dx)
        self.dx2 = self.dx * self.dx
        self.ni = ni
        self.dop = dop
        self.Ldi = Ldi
        self.Ppz_Psp = Ppz_Psp
        self.pol_surf_char = pol_surf_char
        self.Nc = Nc
        self.Nv = Nv
        
        self.q = 1.602176634e-19
        self.kb = 1.380649e-23
        self.T = float(getattr(model, 'T', 300.0))
        self.Vt = self.kb * self.T / self.q
        self.eps0 = 8.8541878128e-12
        
        # Equilibrium potential and carrier references
        self.fi_eq = fi_stat.copy()
        self.ni_ref = 1e18  # standard reference carrier density for normalization (m^-3)
        
        if n_stat is not None and p_stat is not None:
            self.n_eq = np.clip(n_stat, 1e-30, 1e30)
            self.p_eq = np.clip(p_stat, 1e-30, 1e30)
        else:
            self.n_eq = None
            self.p_eq = None
            
        self.eps = model.eps
        self.eps_mid = 2.0 / (1.0 / self.eps[:-1] + 1.0 / self.eps[1:])
        self.coef_eps = (self.Vt / (self.q * self.ni_ref * self.dx2)) * self.eps_mid
        self.scale_t = self.dx2 / self.Vt
        self.J_scale = 0.1 * (self.q * self.Vt * self.ni_ref / self.dx)  # mA/cm^2 scale
        
        self.dfi_eq = self.fi_eq[1:] - self.fi_eq[:-1]
        self.d_psi_n0 = None
        self.d_psi_p0 = None
        self.ni_ratio2 = None
        
        # Dynamic storage for terminal metrics
        self.last_Jtot = 0.0
        self.Jn = None
        self.Jp = None
        self.Jtot = None
        self.mun_mid = None
        self.mup_mid = None
        self.TAUN0 = None
        self.TAUP0 = None
        self.G_norm = None
        self.last_diagnostics = {}

    def require_convergence(self, ok, fi, n, p, Va):
        """Stop the caller before failed/non-finite states enter output arrays."""
        valid = all(np.all(np.isfinite(values)) for values in (fi, n, p))
        if not ok or not valid:
            diagnostics = dict(self.last_diagnostics)
            if not valid:
                diagnostics['reason'] = 'non-finite state'
            raise NewtonConvergenceError(Va, diagnostics)

    def _init_equilibrium_state(self, n_init, p_init):
        if self.n_eq is None:
            self.n_eq = np.clip(n_init, 1e-30, 1e30)
            self.p_eq = np.clip(p_init, 1e-30, 1e30)
        self.d_psi_n0 = np.log(self.n_eq[1:] / self.n_eq[:-1])
        self.d_psi_p0 = -np.log(self.p_eq[1:] / self.p_eq[:-1])
        
        ni_phys = getattr(self.model, 'ni_phys', np.full(self.n_max, 1.8e12))
        self.ni_ratio2 = (ni_phys / self.ni_ref)**2

    def compute_residual_only(self, fi, n, p, Va=0.0):
        N = self.n_max
        F = np.zeros(3 * N)
        
        # Left contact (x=0, Anode)
        # In Aestimo forward bias lowers the anode barrier: fi[0] = fi_eq[0] + Va/Vt
        F[0] = fi[0] - (self.fi_eq[0] + Va / self.Vt)
        F[1] = n[0] - self.n_eq[0]
        F[2] = p[0] - self.p_eq[0]
        
        # Internal field and Bernoulli flux arguments
        dfi = fi[1:] - fi[:-1]
        delta_dfi = dfi - self.dfi_eq
        psi_n = delta_dfi + self.d_psi_n0
        psi_p = delta_dfi + self.d_psi_p0
        
        Bn_pos = Ber(psi_n)
        Bn_neg = Ber(-psi_n)
        Bp_pos = Ber(psi_p)
        Bp_neg = Ber(-psi_p)
        
        Jn = self.mun_mid * (n[1:] * Bn_pos - n[:-1] * Bn_neg)
        Jp = self.mup_mid * (p[:-1] * Bp_pos - p[1:] * Bp_neg)
        
        denom = self.TAUP0 * (n + np.sqrt(self.ni_ratio2)) + self.TAUN0 * (p + np.sqrt(self.ni_ratio2))
        denom = np.maximum(denom, 1e-20)
        U_srh = (n * p - self.ni_ratio2) / denom
        R_minus_G = U_srh - self.G_norm
        
        # Vectorized Poisson residual relative to equilibrium Gauss law
        F_fi = (self.coef_eps[1:] * (dfi[1:] - self.dfi_eq[1:])
              - self.coef_eps[:-1] * (dfi[:-1] - self.dfi_eq[:-1])
              + (p[1:-1] - self.p_eq[1:-1]) - (n[1:-1] - self.n_eq[1:-1]))
        F[3::3][:-1] = F_fi
        
        # Vectorized Continuity residuals
        F_n = (Jn[1:] - Jn[:-1]) - self.scale_t * R_minus_G[1:-1]
        F[4::3][:-1] = F_n
        
        F_p = -(Jp[1:] - Jp[:-1]) - self.scale_t * R_minus_G[1:-1]
        F[5::3][:-1] = F_p
        
        # Right contact (x=L, Cathode)
        r_last = 3 * (N - 1)
        F[r_last] = fi[N-1] - self.fi_eq[N-1]
        F[r_last+1] = n[N-1] - self.n_eq[N-1]
        F[r_last+2] = p[N-1] - self.p_eq[N-1]
        
        return F

    def compute_residual_and_jacobian(self, fi, n, p, Va=0.0):
        N = self.n_max
        F = self.compute_residual_only(fi, n, p, Va=Va)
        
        rows, cols, vals = [], [], []
        def add(r, c, v):
            rows.append(r)
            cols.append(c)
            vals.append(v)
            
        # Node 0 Dirichlet BCs
        add(0, 0, 1.0)
        add(1, 1, 1.0)
        add(2, 2, 1.0)
        
        dfi = fi[1:] - fi[:-1]
        delta_dfi = dfi - self.dfi_eq
        psi_n = delta_dfi + self.d_psi_n0
        psi_p = delta_dfi + self.d_psi_p0
        
        Bn_pos = Ber(psi_n)
        Bn_neg = Ber(-psi_n)
        Bpn_pos = Ber_prime(psi_n)
        Bpn_neg = Ber_prime(-psi_n)
        
        Bp_pos = Ber(psi_p)
        Bp_neg = Ber(-psi_p)
        Bpp_pos = Ber_prime(psi_p)
        Bpp_neg = Ber_prime(-psi_p)
        
        denom = self.TAUP0 * (n + np.sqrt(self.ni_ratio2)) + self.TAUN0 * (p + np.sqrt(self.ni_ratio2))
        denom = np.maximum(denom, 1e-20)
        dU_dn = (p * denom - (n * p - self.ni_ratio2) * self.TAUP0) / (denom**2)
        dU_dp = (n * denom - (n * p - self.ni_ratio2) * self.TAUN0) / (denom**2)
        
        for i in range(1, N - 1):
            r_fi = 3 * i
            r_n = 3 * i + 1
            r_p = 3 * i + 2
            
            # Poisson derivatives
            add(r_fi, 3*(i-1), self.coef_eps[i-1])
            add(r_fi, 3*i, -(self.coef_eps[i] + self.coef_eps[i-1]))
            add(r_fi, 3*(i+1), self.coef_eps[i])
            add(r_fi, 3*i+1, -1.0)
            add(r_fi, 3*i+2, +1.0)
            
            # Electron continuity derivatives
            dJn_i_dfi_i = self.mun_mid[i] * (-n[i+1] * Bpn_pos[i] - n[i] * Bpn_neg[i])
            dJn_im1_dfi_i = self.mun_mid[i-1] * (n[i] * Bpn_pos[i-1] + n[i-1] * Bpn_neg[i-1])
            add(r_n, 3*(i-1), dJn_im1_dfi_i)
            add(r_n, 3*i, dJn_i_dfi_i - dJn_im1_dfi_i)
            add(r_n, 3*(i+1), -dJn_i_dfi_i)
            
            add(r_n, 3*(i-1)+1, self.mun_mid[i-1] * Bn_neg[i-1])
            add(r_n, 3*i+1, -self.mun_mid[i] * Bn_neg[i] - self.mun_mid[i-1] * Bn_pos[i-1] - self.scale_t * dU_dn[i])
            add(r_n, 3*(i+1)+1, self.mun_mid[i] * Bn_pos[i])
            add(r_n, 3*i+2, -self.scale_t * dU_dp[i])
            
            # Hole continuity derivatives (exact signs verified against numerical finite differences)
            dJp_i_dfi_i = self.mup_mid[i] * (-p[i] * Bpp_pos[i] - p[i+1] * Bpp_neg[i])
            dJp_im1_dfi_i = self.mup_mid[i-1] * (p[i-1] * Bpp_pos[i-1] + p[i] * Bpp_neg[i-1])
            add(r_p, 3*(i-1), -dJp_im1_dfi_i)
            add(r_p, 3*i, -(dJp_i_dfi_i - dJp_im1_dfi_i))
            add(r_p, 3*(i+1), dJp_i_dfi_i)
            
            add(r_p, 3*(i-1)+2, self.mup_mid[i-1] * Bp_pos[i-1])
            add(r_p, 3*i+2, -(self.mup_mid[i] * Bp_pos[i] + self.mup_mid[i-1] * Bp_neg[i-1]) - self.scale_t * dU_dp[i])
            add(r_p, 3*(i+1)+2, self.mup_mid[i] * Bp_neg[i])
            add(r_p, 3*i+1, -self.scale_t * dU_dn[i])
            
        # Node N-1 Dirichlet BCs
        r_last = 3 * (N - 1)
        add(r_last, r_last, 1.0)
        add(r_last+1, r_last+1, 1.0)
        add(r_last+2, r_last+2, 1.0)
        
        J_mat = sp.csc_matrix((vals, (rows, cols)), shape=(3*N, 3*N))
        return F, J_mat

    def solve_step(self, fi, n, p, Va=0.0, max_iter=25, tol=0.02):
        if not isinstance(max_iter, (int, np.integer)) or max_iter < 0:
            raise ValueError('max_iter must be a non-negative integer')
        if not np.isfinite(tol) or tol <= 0:
            raise ValueError('tol must be finite and positive')
        # A previous voltage step's currents must not survive a failed solve.
        self.Jn = self.Jp = self.Jtot = None
        self.last_Jtot = float('nan')
        self.last_diagnostics = {'iterations': 0, 'residual_norm': float('inf'),
                                 'reason': 'iteration limit', 'voltage': float(Va)}
        for it in range(max_iter + 1):
            self.last_diagnostics['iterations'] = it
            if not all(np.all(np.isfinite(values)) for values in (fi, n, p)):
                self.last_diagnostics['reason'] = 'non-finite state'
                break
            if np.any(n <= 0) or np.any(p <= 0):
                self.last_diagnostics['reason'] = 'non-positive carrier density'
                break
            F, J = self.compute_residual_and_jacobian(fi, n, p, Va=Va)
            if not np.all(np.isfinite(F)):
                self.last_diagnostics['reason'] = 'non-finite residual'
                break
            res_norm = float(np.max(np.abs(F)))
            self.last_diagnostics['residual_norm'] = res_norm
            self.last_diagnostics['poisson_residual'] = float(np.max(np.abs(F[0::3])))
            self.last_diagnostics['electron_residual'] = float(np.max(np.abs(F[1::3])))
            self.last_diagnostics['hole_residual'] = float(np.max(np.abs(F[2::3])))
            # Small potential corrections alone do not imply continuity convergence.
            if res_norm < tol:
                self.compute_currents(fi, n, p)
                if np.all(np.isfinite(self.Jtot)) and np.isfinite(self.last_Jtot):
                    self.last_diagnostics['reason'] = 'converged'
                    return fi, n, p, True
                self.Jn = self.Jp = self.Jtot = None
                self.last_Jtot = float('nan')
                self.last_diagnostics['reason'] = 'non-finite current'
                break
            if it == max_iter:
                break
            if not np.all(np.isfinite(J.data)):
                self.last_diagnostics['reason'] = 'non-finite Jacobian'
                break
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter('error', spla.MatrixRankWarning)
                    dX = spla.spsolve(J, -F)
            except (spla.MatrixRankWarning, RuntimeError, ValueError, np.linalg.LinAlgError):
                self.last_diagnostics['reason'] = 'linear solve failed'
                break
            if not np.all(np.isfinite(dX)):
                self.last_diagnostics['reason'] = 'non-finite Newton correction'
                break
            d_fi = dX[0::3]
            d_n = dX[1::3]
            d_p = dX[2::3]
            
            max_dfi = np.max(np.abs(d_fi))
            alpha = 1.0
            # Damping on significant majority/minority carriers to prevent inversion
            sig_n = (d_n < 0) & (n > 10.0)
            if np.any(sig_n):
                alpha = min(alpha, float(np.min(-0.85 * n[sig_n] / d_n[sig_n])))
            sig_p = (d_p < 0) & (p > 10.0)
            if np.any(sig_p):
                alpha = min(alpha, float(np.min(-0.85 * p[sig_p] / d_p[sig_p])))
                
            if max_dfi * alpha > 2.0:
                alpha = 2.0 / max_dfi
                
            fi = fi + alpha * d_fi
            n = np.maximum(1e-30, n + alpha * d_n)
            p = np.maximum(1e-30, p + alpha * d_p)
            
        return fi, n, p, False

    def compute_currents(self, fi, n, p):
        dfi = fi[1:] - fi[:-1]
        delta_dfi = dfi - self.dfi_eq
        psi_n = delta_dfi + self.d_psi_n0
        psi_p = delta_dfi + self.d_psi_p0
        
        Bn_pos = Ber(psi_n)
        Bn_neg = Ber(-psi_n)
        Bp_pos = Ber(psi_p)
        Bp_neg = Ber(-psi_p)
        
        Jn = self.mun_mid * (n[1:] * Bn_pos - n[:-1] * Bn_neg)
        Jp = self.mup_mid * (p[:-1] * Bp_pos - p[1:] * Bp_neg)
        
        self.Jn = Jn * self.J_scale
        self.Jp = Jp * self.J_scale
        self.Jtot = self.Jn + self.Jp
        
        mid_lo = int(0.2 * self.n_max)
        mid_hi = int(0.8 * self.n_max)
        self.last_Jtot = float(np.median(self.Jtot[mid_lo:mid_hi]))
        return self.last_Jtot

    def solve(self, fi_init, n_init, p_init, mun, mup, TAUN0, TAUP0, Cn0, Cp0, G_opt, iteration=1, Va=0.0):
        if self.d_psi_n0 is None:
            self._init_equilibrium_state(n_init, p_init)
            
        # Ensure mobilities in SI m^2/Vs
        mun_clean = mun.copy()
        mup_clean = mup.copy()
        if np.max(mun_clean) > 10.0:  # was in cm^2/Vs
            mun_clean *= 1e-4
            mup_clean *= 1e-4
            
        self.mun_mid = (mun_clean[:-1] + mun_clean[1:]) / 2.0
        self.mup_mid = (mup_clean[:-1] + mup_clean[1:]) / 2.0
        self.TAUN0 = TAUN0
        self.TAUP0 = TAUP0
        
        if np.isscalar(G_opt):
            G_opt_arr = np.full(self.n_max, float(G_opt))
        else:
            G_opt_arr = G_opt
        self.G_norm = G_opt_arr / self.ni_ref
        
        fi, n, p, ok = self.solve_step(
            fi_init, n_init, p_init, Va=Va,
            max_iter=getattr(self.model, 'dd_max_iterations', 25),
            tol=getattr(self.model, 'dd_residual_tolerance', 0.02)
        )
        return fi, n, p, ok

