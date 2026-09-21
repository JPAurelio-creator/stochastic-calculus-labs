"""
Verificação numérica do lema de Itô.

Acompanha https://quanttrainer.com.br/lema-de-ito.html

São três afirmações, na ordem em que aparecem na página.

1. A variação quadrática do movimento browniano converge para T, enquanto a
   variação total diverge. É daí que vem o termo de segunda ordem que a regra
   da cadeia usual não tem.

2. Para o movimento browniano geométrico, o lema de Itô dá
   E[S_T] = S_0 * exp(mu*T), enquanto a mediana é S_0 * exp((mu - sigma^2/2)*T).
   A diferença entre as duas é o termo de Itô, e ela cresce com a volatilidade.

3. O drift da log-trajetória é mu - sigma^2/2, não mu. Conferido por simulação.
"""

import numpy as np

T = 1.0
SEED = 20260828


def variacoes(n, rng):
    """Variação quadrática e variação total de uma trajetória com n passos."""
    dt = T / n
    incrementos = rng.standard_normal(n) * np.sqrt(dt)
    return np.sum(incrementos**2), np.sum(np.abs(incrementos))


def parte_1():
    print("1. Variação quadrática contra variação total, uma trajetória")
    print(f"{'passos':>10}  {'var. quadrática':>16}  {'var. total':>12}  "
          f"{'sqrt(2n/pi)':>12}")
    print("-" * 56)
    rng = np.random.default_rng(SEED)
    for n in [100, 1_000, 10_000, 100_000, 1_000_000]:
        vq, vt = variacoes(n, rng)
        # A variação total de n incrementos gaussianos tem valor esperado
        # n * E|N(0, dt)| = n * sqrt(2*dt/pi) = sqrt(2n/pi).
        print(f"{n:>10,}  {vq:>16.3f}  {vt:>12.1f}  {np.sqrt(2*n/np.pi):>12.1f}"
              .replace(",", " "))
    print("-" * 56)
    print(f"A coluna do meio converge para T = {T:g}. A da direita cresce "
          f"como sqrt(n) e não converge.\n")
    assert abs(vq - T) < 0.01, "variação quadrática não convergiu para T"


def parte_2():
    print("2. Média contra mediana do movimento browniano geométrico")
    print(f"{'mu':>5}  {'sigma':>6}  {'E[S_T]/S_0':>12}  {'exp(mu T)':>10}  "
          f"{'mediana/S_0':>12}  {'exp((mu-s^2/2)T)':>17}")
    print("-" * 72)
    rng = np.random.default_rng(SEED)
    n_sim = 4_000_000
    for mu, sigma in [(0.08, 0.20), (0.10, 0.30), (0.15, 0.50)]:
        w_t = rng.standard_normal(n_sim) * np.sqrt(T)
        s_t = np.exp((mu - 0.5 * sigma**2) * T + sigma * w_t)  # com S_0 = 1
        media_teorica = np.exp(mu * T)
        mediana_teorica = np.exp((mu - 0.5 * sigma**2) * T)
        print(f"{mu:>5.0%}  {sigma:>6.0%}  {s_t.mean():>12.3f}  "
              f"{media_teorica:>10.3f}  {np.median(s_t):>12.3f}  "
              f"{mediana_teorica:>17.3f}")
        assert abs(s_t.mean() - media_teorica) < 0.01
        assert abs(np.median(s_t) - mediana_teorica) < 0.01
    print("-" * 72)
    print("A média cresce com mu. A mediana cresce com mu - sigma^2/2.")
    print("O buraco entre as duas é o termo de Itô, e ele abre com sigma.\n")


def parte_3():
    print("3. Drift da log-trajetória, simulado contra teórico")
    mu, sigma, n_sim = 0.10, 0.30, 4_000_000
    rng = np.random.default_rng(SEED)
    w_t = rng.standard_normal(n_sim) * np.sqrt(T)
    log_s = (mu - 0.5 * sigma**2) * T + sigma * w_t

    drift_simulado = log_s.mean() / T
    erro_padrao = log_s.std(ddof=1) / np.sqrt(n_sim) / T
    drift_ito = mu - 0.5 * sigma**2

    print(f"  mu = {mu:.0%}, sigma = {sigma:.0%}")
    print(f"  E[log(S_T/S_0)]/T simulado ... {drift_simulado:.5f} "
          f"(erro padrão {erro_padrao:.5f})")
    print(f"  mu - sigma^2/2 .............. {drift_ito:.5f}")
    print(f"  mu, a resposta errada ....... {mu:.5f}")
    z_ito = abs(drift_simulado - drift_ito) / erro_padrao
    z_mu = abs(drift_simulado - mu) / erro_padrao
    print(f"  Afastamento de mu - sigma^2/2, em erros padrão: {z_ito:.2f}")
    print(f"  Afastamento de mu, em erros padrão: {z_mu:.2f}")
    assert z_ito < 4.0, "simulação discorda do drift de Itô"
    assert z_mu > 10.0, "simulação não consegue separar as duas hipóteses"
    print("\nOK, a simulação concorda com Itô e rejeita mu com folga.")


if __name__ == "__main__":
    parte_1()
    parte_2()
    parte_3()
