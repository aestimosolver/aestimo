import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import scipy.optimize
import logging
import copy
from aeslibs.aestimo_poisson1d import Ber

logger = logging.getLogger(__name__)

class CoupledNewtonSolver:
    def __init__(self, model, n_max, dx, ni, dop, Ldi, Ppz_Psp, pol_surf_char, Nc, Nv, fi_stat):
        self.model = model
        self.n_max = n_max
        self.dx = dx
        self.dx_m = dx * 1e-9
        self.dx2 = self.dx_m**2
        self.ni = ni
        self.ni_phys = getattr(model, 'ni_phys', ni)
        self.dop = dop
        self.Ldi = Ldi
        self.Ppz_Psp = Ppz_Psp
        self.pol_surf_char = pol_surf_char
        self.Nc = Nc
        self.Nv = Nv
        self.fi_stat = fi_stat
        
        self.Vt = 1.380649e-23 * getattr(model, 'T', 300.0) / 1.602176634e-19
        self.q = 1.602176634e-19
        self.eps0 = 8.8541878128e-12
        
    def solve(self, fi_init, n_init, p_init, mun, mup, TAUN0, TAUP0, Cn0, Cp0, G_opt, iteration):
        # We pack variables into X = [fi_0, n_0, p_0, fi_1, n_1, p_1, ...]
        # fi is normalized potential
        # n and p are normalized carrier densities (divided by ni)
        
        X_init = np.zeros(3 * self.n_max)
        X_init[0::3] = fi_init
        X_init[1::3] = n_init
        X_init[2::3] = p_init
        
        # We will use scipy.optimize.root with Jacobian-Free Newton-Krylov
        # This completely avoids assembling the massive analytical Jacobian!
        print(f"--> Starting Fully-Coupled Newton-Krylov solver (iteration {iteration})...")
        logger.info(f"    Starting Fully-Coupled Newton-Krylov solver (iteration {iteration})...")
        
        # Scaling vectors to balance the variables for the solver
        # fi is O(1) - O(100)
        # n, p can be O(1e18) / ni, which is O(1e6) - O(1e12)
        # We solve for log(n) and log(p) internally to keep variables O(1) and strictly positive!
        
        # Solve for fi, n, p directly to avoid exponential scaling issues
        Y_init = X_init.copy()
        
        def residual(Y):
            fi = Y[0::3]
            n = Y[1::3]
            p = Y[2::3]
            
            F = np.zeros_like(Y)
            F_fi = F[0::3]
            F_n = F[1::3]
            F_p = F[2::3]
            
            # --- Poisson Residual ---
            Ldi2 = self.Ldi * self.Ldi
            dop_out = self.dop / self.ni
            Ppz_Psp_out = self.Ppz_Psp / self.ni if self.Ppz_Psp is not None else np.zeros(self.n_max)
            
            for i in range(1, self.n_max - 1):
                pol_charge = 0.0
                if self.Ppz_Psp is not None and i < self.n_max - 2:
                    pol_charge = (self.Ppz_Psp[i] - self.Ppz_Psp[i+1]) / self.dx_m / self.ni[i]
                if self.pol_surf_char is not None and self.pol_surf_char[i] != 0:
                    pol_charge = self.pol_surf_char[i] / self.ni[i]
                    
                # Aestimo uses Ldi2 / dx2 for the normalized Poisson equation
                coef = Ldi2[i] / self.dx2
                charge_scale = max(1.0, abs(dop_out[i]))
                F_fi[i] = (coef * (fi[i-1] - 2*fi[i] + fi[i+1]) + p[i] - n[i] + dop_out[i] + pol_charge) / charge_scale
                
            # Poisson Boundaries (Dirichlet)
            F_fi[0] = fi[0] - fi_init[0]
            F_fi[-1] = fi[-1] - fi_init[-1]
            
            # --- Continuity Residuals ---
            # Bernoulli function coefficients
            for i in range(1, self.n_max - 1):
                mun_avg_m1 = (mun[i-1] + mun[i])/2.0
                mun_avg_p1 = (mun[i] + mun[i+1])/2.0
                mup_avg_m1 = (mup[i-1] + mup[i])/2.0
                mup_avg_p1 = (mup[i] + mup[i+1])/2.0
                
                dfi_m1 = fi[i-1] - fi[i]
                dfi_p1 = fi[i+1] - fi[i]
                
                J_n_m1 = mun_avg_m1 * (Ber(dfi_m1)*n[i-1] - Ber(-dfi_m1)*n[i])
                J_n_p1 = mun_avg_p1 * (Ber(-dfi_p1)*n[i] - Ber(dfi_p1)*n[i+1])
                div_Jn = J_n_m1 - J_n_p1
                
                J_p_m1 = mup_avg_m1 * (Ber(-dfi_m1)*p[i-1] - Ber(dfi_m1)*p[i])
                J_p_p1 = mup_avg_p1 * (Ber(dfi_p1)*p[i] - Ber(-dfi_p1)*p[i+1])
                div_Jp = J_p_m1 - J_p_p1
                
                # Recombination
                ni_r = self.ni[i]
                ni_p = self.ni_phys[i]
                ni_ratio2 = (ni_p / ni_r)**2
                
                # SRH + TAT
                trap_scale = max(getattr(self.model, 'trap_density_scale', 1.0), 1e-12)
                E_field = (self.Vt / self.dx_m) * abs((fi[i+1] - fi[i-1]) / 2.0)
                gamma = 0.0
                tat_field = float(getattr(self.model, 'tat_field', 1e10))
                if E_field > 1e4 and tat_field < 1e9:
                    ratio = E_field / tat_field
                    gamma = 2.0 * np.sqrt(3.0 * np.pi) * ratio * np.exp(np.clip(ratio**2, 0, 20))
                
                gamma *= trap_scale
                n1_norm = np.maximum((ni_p / ni_r), 1e-20)
                p1_norm = np.maximum((ni_p / ni_r), 1e-20)
                
                # Use np.abs for safe density in denominator to avoid negative roots
                n_safe = np.abs(n[i])
                p_safe = np.abs(p[i])
                
                denom = TAUP0[i]*(n_safe + n1_norm) + TAUN0[i]*(p_safe + p1_norm)
                denom = max(denom / max(1.0 + gamma, 1e-12), 1e-20)
                
                U_srh = (n_safe*p_safe - ni_ratio2) / denom
                
                # Auger
                U_aug = 0.0
                if Cn0 is not None and Cp0 is not None:
                    U_aug = (Cn0[i]*n_safe + Cp0[i]*p_safe) * (ni_r**2) * (n_safe*p_safe - ni_ratio2)
                    
                # Optical Gen
                G_op = G_opt[i] / ni_r if isinstance(G_opt, np.ndarray) else G_opt / ni_r
                
                U_total = U_srh + U_aug - G_op
                
                scale = (self.dx2 / self.Vt)
                
                F_n[i] = scale * U_total - div_Jn
                F_p[i] = scale * U_total - div_Jp
                
            # Continuity Boundaries (Dirichlet)
            F_n[0] = Y[1] - Y_init[1]
            F_n[-1] = Y[-2] - Y_init[-2]
            F_p[0] = Y[2] - Y_init[2]
            F_p[-1] = Y[-1] - Y_init[-1]
            
            return F
            
        # Call Newton-Krylov solver (Jacobian-Free Newton-Krylov GMRES)
        try:
            sol = scipy.optimize.root(residual, Y_init, method='krylov', options={'fatol': 1e-2, 'maxiter': 50, 'disp': False})
            if sol is not None and sol.success:
                Y_final = sol.x
                success_flag = True
            else:
                logger.warning(f"    Newton solver fallback: maintaining initial step values.")
                Y_final = Y_init
                success_flag = False
        except Exception as e:
            logger.error(f"    Newton solver exception: {e}")
            Y_final = Y_init
            success_flag = False
            
        fi_out = Y_final[0::3]
        n_out = np.abs(Y_final[1::3])
        p_out = np.abs(Y_final[2::3])
        
        return fi_out, n_out, p_out, success_flag
