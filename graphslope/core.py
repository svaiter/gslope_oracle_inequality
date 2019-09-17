import numpy as np
import scipy.linalg as lin
from scipy.sparse.linalg import svds
from sklearn.isotonic import isotonic_regression


def prox_l1sorted(v, w):
    """Proximal operator of the OWL norm dot(w, reversed(sort(v)))
    Follows description and notation from:
    X. Zeng, M. Figueiredo,
    The ordered weighted L1 norm: Atomic formulation, dual norm,
    and projections.
    eprint http://arxiv.org/abs/1409.4271
    (From pyowl)
    XXX
    """

    v_abs = np.abs(v)
    ix = np.argsort(v_abs)[::-1]
    v_abs = v_abs[ix]
    # project to K+ (monotone non-negative decreasing cone)
    v_abs = isotonic_regression(v_abs - w, y_min=0, increasing=False)

    # undo the sorting
    inv_ix = np.zeros_like(ix)
    inv_ix[ix] = np.arange(len(v))
    v_abs = v_abs[inv_ix]

    return np.sign(v) * v_abs


def graphslope_fista(D, y, lams, operator_norm=None,
                     warm_start=None, tol=1e-4, maxiter=1000, retall=False):
    """ Solve the problem
    argmin 1/2 ||y - beta||_2^2 + SLOPE_lams(D^T beta)
    """
    hnormy = 0.5 * lin.norm(y) ** 2
    def prox_G(beta, tau):
        return beta - tau * prox_l1sorted(beta/tau, 1./tau * lams)
    def grad_F(theta):
        return D.T.dot(D.dot(theta) - y)
    def primal_energy(beta):
        return 0.5 * lin.norm(y-beta)**2 + np.sum(lams * np.sort(np.abs(D.T.dot(beta)))[::-1])
    def dual_energy(theta):
        return 0.5 * lin.norm(D.dot(theta) - y) ** 2 - hnormy

    L = operator_norm if operator_norm else svds(D,k=1)[1][0] ** 2

    if not warm_start:
        theta = np.zeros(D.shape[1])
        beta = np.zeros_like(y)
    else:
        (beta, theta) = warm_start
    primal, dual = primal_energy(beta), dual_energy(prox_G(theta,1.))
    gap = primal + dual
    z = theta.copy()

    iter = 1
    t = 1. # Initial relaxation for FISTA

    primal_list = []
    dual_list = []
    gap_list = []
    while (gap > tol) and (iter < maxiter):
        theta_new = prox_G(z - 1. / L * grad_F(z), 1. / L)
        beta = y - D.dot(theta_new)
        primal, dual = primal_energy(beta), dual_energy(theta_new)
        gap = primal + dual

        tnew = (1. + np.sqrt(1. + 4. * t ** 2)) / 2.
        z = theta_new + (t - 1.) / tnew * (theta_new - theta)
        theta = theta_new
        t = tnew
        if retall:
            primal_list.append(primal)
            dual_list.append(dual)
            gap_list.append(gap)
        iter += 1

    if iter == maxiter:
        print("Convergence not achieved. Gap: {0}. Tol: {1}".format(gap, tol))
    if retall:
        dict_res = {
            'beta': beta,
            'theta': theta,
            'final_iter': iter,
            'final_gap': gap,
            'primal_list': np.asarray(primal_list),
            'dual_list': np.asarray(dual_list),
            'gap_list': np.asarray(gap_list)
        }
        return dict_res
    else:
        return beta


class GraphSlope:
    """Graph-Slope regression model.

    The Graph-Slope regression optimization problem is defined by::

        min 1/2n ||y - beta||_2^2 + norm_slope(D^T beta, lambdas)

        where norm_slope(., lambdas) is the reorded l1-norm with weights
        lambdas.

    Argument
    --------
    incidence : ndarray
        The incidence matrix D of the underlying graph.

    Parameters
    ----------
    lambdas : ndarray or float, default: 1.0
        The weights associated to the slope norm. Ensure they
        are decaying.

    warm_start : boolean, default: False
        If True and self.coef_ and self.dual_coef_ exist, those values are
        used to start the optimization algorithm.

    incidence_norm : float, default: None
        If None, the matrix norm of the incidence matrix is computed,
        otherwise used as a value in the optimization algorithm.

    max_iter : integer, default:1000
        Maximum number of iterations.

    tol : float, default:0.0001
        Dual gap tolerance to stop the optimization.

    Attributes
    ----------
    coef_ : ndarray
        parameter (primal) vector (beta in the cost function).

    dual_coef_ : ndarray
        parameter (dual) vector.

    Examples
    --------
    XXX

    References
    ----------
    A sharp oracle inequality for Graph-Slope,
    Pierre Bellec, Joseph Salmon and Samuel Vaiter
    2017.

    The ordered weighted L1 norm: Atomic formulation, dual norm,
    and projections.
    X. Zeng, M. Figueiredo.
    2014.

    SLOPE - Adaptive Variable Selection via Convex Optimization.
    Malgorzata Bogdan, Ewout van den Berg, Weijie Su, Emmanuel Candes.
    2014.
    """

    def __init__(self, incidence, lambdas=1.0, incidence_norm=None, warm_start=False, max_iter=1000, tol=0.0001):
        self.incidence = incidence
        self.lambdas = lambdas
        self.warm_start = warm_start
        self.max_iter = max_iter
        self.tol = tol
        self.incidence_norm = incidence_norm

        self.coef_ = None
        self.dual_coef_ = None

    def fit(self, y):

        self.y_ = y

        if not self.incidence_norm:
            self.incidence_norm = svds(self.incidence,k=1)[1][0] ** 2
        if np.isscalar(self.lambdas):
            self.lambdas = self.lambdas * np.ones(self.incidence.shape[1])

        if self.warm_start and self.coef_ is not None and self.dual_coef_ is not None:
            warmed_up = (self.coef_, self.dual_coef_)
        else:
            warmed_up = None

        # Beware, multiply lambdas by n to go from 1/2n ||.||^2 to 1/2 ||.||^2
        res = graphslope_fista(self.incidence, self.y_, y.shape[0] * self.lambdas,
                               operator_norm=self.incidence_norm, warm_start=warmed_up,
                               tol=self.tol, maxiter=self.max_iter, retall=True)

        self.coef_ = res['beta']
        self.dual_coef_ = res['theta']
        self.dual_gaps_ = res['gap_list']
        self.final_iter_ = res['final_iter']
        self.final_gap_ = res['final_gap']
        self.raw_output_ = res

        return self
