"""
Verificação numérica da redução de Black-Scholes à equação do calor.

Acompanha https://quanttrainer.com.br/black-scholes-equacao-do-calor.html

A afirmação testada é que a mudança de variáveis da página não é uma analogia,
é uma identidade. Resolvendo a equação do calor pura por diferenças finitas e
desfazendo a substituição, o número que sai tem que ser o mesmo da fórmula
fechada de Black-Scholes, dentro do erro de discretização.

A cadeia de substituições é a clássica.

    x = ln(S/K),   tau = sigma^2 (T - t) / 2,   k = 2r / sigma^2
    V(S, t) = K * v(x, tau)
    v(x, tau) = exp(-(k-1)x/2 - (k+1)^2 tau/4) * u(x, tau)

e sobra u_tau = u_xx, que é a equação do calor sem nenhum parâmetro. Toda a
estrutura financeira, taxa de juros, volatilidade e prazo, fica dentro das
substituições e nenhuma dentro da equação.

A integração usa Crank-Nicolson, que é incondicionalmente estável e de segunda
ordem no tempo e no espaço.
"""

import numpy as np
from scipy.linalg import solve_banded
from scipy.stats import norm

K = 100.0          # strike, fixo, porque o erro é medido como % do strike
X_MIN, X_MAX = -5.0, 5.0
N_X = 2000
N_TAU = 2000


def black_scholes_call(s, k, r, sigma, t):
    d1 = (np.log(s / k) + (r + 0.5 * sigma**2) * t) / (sigma * np.sqrt(t))
    d2 = d1 - sigma * np.sqrt(t)
    return s * norm.cdf(d1) - k * np.exp(-r * t) * norm.cdf(d2)


def u_contorno(x, tau, k):
    """Valor assintótico de u para x grande, vindo da própria substituição."""
    return (np.exp(0.5 * (k + 1) * x + 0.25 * (k + 1) ** 2 * tau)
            - np.exp(0.5 * (k - 1) * x + 0.25 * (k - 1) ** 2 * tau))


def resolve_equacao_do_calor(r, sigma, t_venc):
    """Integra u_tau = u_xx por Crank-Nicolson e devolve (x, u_final, k)."""
    k = 2.0 * r / sigma**2
    tau_max = 0.5 * sigma**2 * t_venc

    x = np.linspace(X_MIN, X_MAX, N_X + 1)
    dx = x[1] - x[0]
    dtau = tau_max / N_TAU
    lam = dtau / dx**2

    # Condição inicial, que é o payoff da call já transformado.
    u = np.maximum(np.exp(0.5 * (k + 1) * x) - np.exp(0.5 * (k - 1) * x), 0.0)

    # Matriz tridiagonal de Crank-Nicolson, em formato de banda.
    a = 0.5 * lam
    banda = np.zeros((3, N_X - 1))
    banda[0, 1:] = -a
    banda[1, :] = 1.0 + 2.0 * a
    banda[2, :-1] = -a

    for n in range(N_TAU):
        tau_atual, tau_novo = n * dtau, (n + 1) * dtau
        interior = u[1:-1]
        lado_direito = (interior
                        + a * (u[2:] - 2.0 * interior + u[:-2]))

        # Contorno em x_min, onde a opção não vale nada, e em x_max, onde ela
        # vira o ativo menos o strike descontado.
        u_esq_novo, u_esq_velho = 0.0, 0.0
        u_dir_novo = u_contorno(X_MAX, tau_novo, k)
        u_dir_velho = u_contorno(X_MAX, tau_atual, k)

        lado_direito[0] += a * (u_esq_novo - u_esq_velho) + a * u_esq_novo
        lado_direito[-1] += a * (u_dir_novo - u_dir_velho) + a * u_dir_novo

        u = np.concatenate(([u_esq_novo],
                            solve_banded((1, 1), banda, lado_direito),
                            [u_dir_novo]))

    return x, u, k, tau_max


def preco_por_diferencas_finitas(s_alvos, r, sigma, t_venc):
    """Desfaz as substituições e devolve o preço em unidades de dinheiro."""
    x, u, k, tau_max = resolve_equacao_do_calor(r, sigma, t_venc)
    v = np.exp(-0.5 * (k - 1) * x - 0.25 * (k + 1) ** 2 * tau_max) * u
    precos = K * v
    x_alvos = np.log(np.asarray(s_alvos) / K)
    return np.interp(x_alvos, x, precos)


def main():
    print("Equação do calor por diferenças finitas contra Black-Scholes")
    print(f"Grade de {N_X} pontos em x e {N_TAU} passos em tau, "
          f"Crank-Nicolson\n")
    print(f"{'r':>5}  {'sigma':>6}  {'T':>5}  {'maior desvio':>13}  "
          f"{'% do strike':>12}")
    print("-" * 50)

    s_alvos = [80.0, 90.0, 100.0, 110.0, 120.0]
    pior_absoluto, pior_relativo = 0.0, 0.0

    for r in [0.0, 0.05, 0.10, 0.15]:
        for sigma in [0.20, 0.30, 0.45]:
            for t_venc in [0.5, 1.0]:
                numerico = preco_por_diferencas_finitas(s_alvos, r, sigma,
                                                        t_venc)
                fechado = black_scholes_call(np.array(s_alvos), K, r, sigma,
                                             t_venc)
                desvio = np.max(np.abs(numerico - fechado))
                pior_absoluto = max(pior_absoluto, desvio)
                pior_relativo = max(pior_relativo, desvio / K)
                print(f"{r:>5.0%}  {sigma:>6.0%}  {t_venc:>5.1f}  "
                      f"{desvio:>13.5f}  {desvio / K:>12.4%}")

    print("-" * 50)
    print(f"Maior desvio em todo o varrido: {pior_absoluto:.5f}, "
          f"ou {pior_relativo:.4%} do strike")
    assert pior_relativo < 0.0001, "desvio acima do erro de discretização"
    print("\nOK, as duas rotas chegam ao mesmo número.")


if __name__ == "__main__":
    main()
