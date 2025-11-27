import numpy as np
import pandas as pd
from scipy.stats import multivariate_normal
from numpy.linalg import multi_dot
import scipy.linalg
import os
from joblib import Parallel, delayed
import multiprocessing

p = 100
n_min = 6
n_max = 30
r_new_list = [0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
np.random.seed(1234)


def get_err(cov_est, cov_t):
    error = np.sum((cov_est - cov_t) ** 2) / (100 ** 2)
    return error


for r in r_new_list[6:11]:
    corr = np.zeros((p, p))
    for i in range(p):
        for j in range(p):
            corr[i, j] = r ** abs(i - j)

    sd_s = 1
    sd_b = np.sqrt(10)
    sd = np.diag(np.concatenate((np.repeat(sd_s, p // 2), np.repeat(sd_b, p // 2))))
    cov = multi_dot([sd, corr, sd])

    def simulation(n):
        def sim_err(iter):
            err = pd.DataFrame(np.zeros((1, 8)),
                               columns=['oracle_c', 'oracle_v', 'lw', 'rblw', 'oas_c', 'oas_v', 'corr_oracle',
                                        'corr_oas'])
            sample = multivariate_normal.rvs(np.repeat(0, p), cov, size=n)

            scov = np.cov(sample, rowvar=False)
            f1 = (np.sum(np.diag(scov)) / p) * np.eye(p)
            f2 = np.diag(np.diag(scov))
            tr_s = np.sum(np.diag(scov))
            tr_s2 = np.sum(np.diag(np.dot(scov.T, scov)))
            tr_t = np.sum(np.diag(cov))
            tr_t2 = np.sum(np.diag(np.dot(cov.T, cov)))
            tr_tdiag2 = np.sum(np.diag(np.diag(np.diag(cov)) * np.diag(np.diag(cov))))
            tr_sdiag2 = np.sum(np.diag(np.diag(np.diag(scov)) * np.diag(np.diag(scov))))
            rou_o_constant = ((1 - 2 / p) * tr_t2 + tr_t ** 2) / ((n + 1 - 2 / p) * tr_t2 + (1 - n / p) * tr_t ** 2)
            scov_o_constant = (1 - rou_o_constant) * scov + rou_o_constant * f1
            rou_o_variant = ((1 / n) * tr_t2 - (2 / n) * tr_tdiag2 + (1 / n) * tr_t ** 2) / (
                    ((n + 1) / n) * tr_t2 + (1 / n) * tr_t ** 2 - ((n + 2) / n) * tr_tdiag2)
            scov_o_variant = (1 - rou_o_variant) * scov + rou_o_variant * f2

            num = 0
            for i in range(n):
                fill = np.sum((np.outer(sample[i, :], sample[i, :]) - scov) ** 2)
                num = num + fill

            den = (n ** 2) * (tr_s2 - (tr_s ** 2) / p)
            rou_lw = num / den
            rou_lw = min(rou_lw, 1)
            scov_lw = (1 - rou_lw) * scov + rou_lw * f1

            rou_rblw = (((n - 2) / n) * tr_s2 + tr_s ** 2) / ((n + 2) * (tr_s2 - (tr_s ** 2) / p))
            rou_rblw = min(rou_rblw, 1)
            scov_rblw = (1 - rou_rblw) * scov + rou_rblw * f1

            phi = (tr_s2 - (tr_s ** 2) / p) / (tr_s ** 2 + (1 - 2 / p) * tr_s2)
            rou_oas_constant = min(1 / ((n + 1 - 2 / p) * phi), 1)
            scov_oas_constant = (1 - rou_oas_constant) * scov + rou_oas_constant * f1

            phi1 = (tr_s2 - tr_sdiag2) / (tr_s2 + tr_s ** 2 - 2 * tr_sdiag2)
            rou_oas_variant = min(1 / ((n + 1) * phi1), 1)
            scov_oas_variant = (1 - rou_oas_variant) * scov + rou_oas_variant * f2

            scorr = scipy.linalg.cov2cor(scov)
            corr = scipy.linalg.cov2cor(cov)
            tr_t2_cor = np.sum(np.diag(np.dot(corr.T, corr)))
            tr_t_cor = np.sum(np.diag(corr))
            tr_s_cor = np.sum(np.diag(scorr))
            tr_s2_cor = np.sum(np.diag(np.dot(scorr.T, scorr)))
            f1_cor = (np.sum(np.diag(scorr)) / p) * np.eye(p)

            rou_o_constant_cor = ((1 - 2 / p) * tr_t2_cor + tr_t_cor ** 2) / (
                    (n + 1 - 2 / p) * tr_t2_cor + (1 - n / p) * tr_t_cor ** 2)
            scorr_o_constant_cor = (1 - rou_o_constant_cor) * scorr + rou_o_constant_cor * f1_cor
            scov_o_constant_cor = np.dot(np.dot(np.diag(np.sqrt(np.diag(scov))), scorr_o_constant_cor),
                                          np.diag(np.sqrt(np.diag(scov))))

            phi_cor = (tr_s2_cor - (tr_s_cor ** 2) / p) / (tr_s_cor ** 2 + (1 - 2 / p) * tr_s2_cor)
            rou_oas_constant_cor = min(1 / ((n + 1 - 2 / p) * phi_cor), 1)
            scorr_oas_constant = (1 - rou_oas_constant_cor) * scorr + rou_oas_constant_cor * f1_cor
            scov_oas_constant_cor = np.dot(np.dot(np.diag(np.sqrt(np.diag(scov))), scorr_oas_constant),
                                            np.diag(np.sqrt(np.diag(scov))))

            err.loc[0, 'oracle_c'] = get_err(scov_o_constant, cov)
            err.loc[0, 'oracle_v'] = get_err(scov_o_variant, cov)
            err.loc[0, 'lw'] = get_err(scov_lw, cov)
            err.loc[0, 'rblw'] = get_err(scov_rblw, cov)
            err.loc[0, 'oas_c'] = get_err(scov_oas_constant, cov)
            err.loc[0, 'oas_v'] = get_err(scov_oas_variant, cov)
            err.loc[0, 'corr_oracle'] = get_err(scov_o_constant_cor, cov)
            err.loc[0, 'corr_oas'] = get_err(scov_oas_constant_cor, cov)

            return err

        num_cores = multiprocessing.cpu_count()
        res = Parallel(n_jobs=num_cores)(
            delayed(sim_err)(i) for i in range(1, 1001))
        res = pd.concat(res).assign(sample_size=n)
        return res

    res_sample_size = pd.concat([simulation(n) for n in range(n_min, n_max + 1)])
    res_sample_size.to_pickle(f'r={r}.pkl')
