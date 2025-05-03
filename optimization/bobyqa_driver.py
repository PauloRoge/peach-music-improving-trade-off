# bobyqa_driver.py — versão estável + fallback correto (best_x)
import numpy as np
import pybobyqa

print(">>> BOBYQA driver (best_x) carregado")

def bobyqa_optimize(URA, Un, ref, wavelength,
                    x0, lower, upper,
                    tol=1e-4, maxiter=500):

    # ------------- conversões ----------------
    URA   = np.asarray(URA, dtype=float)
    Un    = np.asarray(Un,  dtype=complex)
    ref   = np.asarray(ref, dtype=float).ravel()
    x0    = np.asarray(x0,  dtype=float).ravel()
    lower = np.asarray(lower).ravel()
    upper = np.asarray(upper).ravel()
    bounds = np.column_stack((lower, upper))      # shape (n_dim,2)

    k  = 2*np.pi / float(wavelength)
    Pn = Un @ Un.conj().T                         # projeção ruído

    def cost(x):
        xv, yv = x
        d_ref = np.sqrt((xv-ref[0])**2 + (yv-ref[1])**2 + ref[2]**2)
        d_k   = np.sqrt((URA[:,0]-xv)**2 +
                        (URA[:,1]-yv)**2 +
                         URA[:,2]**2)
        phase = -k * (d_ref - d_k)
        a     = np.exp(1j*phase)
        denom = np.real(a.conj() @ (Pn @ a))
        return denom + 1e-12          # custo sempre finito > 0

    # ------------- solver --------------------
    res = pybobyqa.solve(cost,                     # x0 é POSICIONAL
                         x0,
                         bounds=bounds,
                         rhobeg=0.3*np.min(upper-lower),
                         rhoend=tol,
                         maxfun=int(maxiter))

    # -------- fallback seguro ---------------
    if res.x is None:
        print(">>> BOBYQA não convergiu – devolvendo best_x ou x0")
        cand = getattr(res, 'best_x', x0)
        return np.asarray(cand, dtype=float)

    return np.asarray(res.x, dtype=float)
