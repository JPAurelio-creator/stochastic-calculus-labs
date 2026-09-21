"""
Verificação numérica do teorema de Girsanov.

Acompanha https://quanttrainer.com.br/teorema-de-girsanov.html

A afirmação testada é a seguinte. O preço de uma opção europeia não depende do
drift real mu do ativo. Simulando sob quatro mundos diferentes e reponderando
cada trajetória pela derivada de Radon-Nikodym

    dQ/dP = exp(-lambda * W_T - 0.5 * lambda^2 * T),    lambda = (mu - r) / sigma,

o preço reponderado tem que bater com a fórmula fechada de Black-Scholes nos
quatro casos.

Também confere a outra face da mesma moeda, que o valor esperado reponderado do
ativo é S_0 * exp(r*T) qualquer que seja o mundo de origem.
"""

import numpy as np
from scipy.stats import norm

S0, K, R, SIGMA, T = 100.0, 105.0, 0.04, 0.25, 1.0
N_SIM = 6_000_000
DRIFTS = [-0.05, 0.04, 0.12, 0.25]
SEED = 20260910


def black_scholes_call(s0, k, r, sigma, t):
    d1 = (np.log(s0 / k) + (r + 0.5 * sigma**2) * t) / (sigma * np.sqrt(t))
    d2 = d1 - sigma * np.sqrt(t)
    return s0 * norm.cdf(d1) - k * np.exp(-r * t) * norm.cdf(d2)


def preco_reponderado(mu, rng):
    """Simula sob P (drift real mu) e reponderá para Q."""
    w_t = rng.standard_normal(N_SIM) * np.sqrt(T)
    s_t = S0 * np.exp((mu - 0.5 * SIGMA**2) * T + SIGMA * w_t)

    lam = (mu - R) / SIGMA
    # densidade de Radon-Nikodym dQ/dP avaliada em cada trajetória
    rn = np.exp(-lam * w_t - 0.5 * lam**2 * T)

    payoff = np.maximum(s_t - K, 0.0)
    amostra = np.exp(-R * T) * payoff * rn
    preco = np.mean(amostra)
    # Erro padrão da própria amostra. Ele não é o mesmo em todas as linhas,
    # e é por isso que ele precisa ser calculado e não chutado. Ver nota
    # sobre variância no fim do arquivo.
    erro_padrao = amostra.std(ddof=1) / np.sqrt(N_SIM)
    esperado_q = np.mean(s_t * rn)
    return lam, preco, erro_padrao, esperado_q


def main():
    referencia = black_scholes_call(S0, K, R, SIGMA, T)
    alvo_esperanca = S0 * np.exp(R * T)

    print(f"S0={S0:g}  K={K:g}  r={R:.0%}  sigma={SIGMA:.0%}  T={T:g} ano")
    print(f"{N_SIM:,} simulações por linha, semente {SEED}\n".replace(",", " "))
    print(f"{'mu':>6}  {'lambda':>8}  {'MC reponderado':>15}  {'erro pad.':>10}  "
          f"{'desvio/erro':>12}  {'E_Q[S_T]':>10}")
    print("-" * 74)

    rng = np.random.default_rng(SEED)
    pior_z = 0.0
    for mu in DRIFTS:
        lam, preco, erro_padrao, esperado_q = preco_reponderado(mu, rng)
        z = (preco - referencia) / erro_padrao
        pior_z = max(pior_z, abs(z))
        print(f"{mu:>6.0%}  {lam:>8.3f}  {preco:>15.4f}  {erro_padrao:>10.4f}  "
              f"{z:>12.2f}  {esperado_q:>10.2f}")

    print("-" * 74)
    print(f"Fórmula fechada de Black-Scholes: {referencia:.4f}")
    print(f"E_Q[S_T] tem que dar S0*exp(r*T) = {alvo_esperanca:.3f}")
    print(f"Maior afastamento em erros padrão: {pior_z:.2f}")

    # O teste certo não é o desvio absoluto, é o desvio medido em erros padrão.
    # Abaixo de 4 desvios o resultado é ruído de Monte Carlo, acima disso é
    # sinal de que alguma coisa está errada na conta.
    assert pior_z < 4.0, "afastamento grande demais para ser ruído de simulação"
    print("\nOK, os quatro mundos concordam com a fórmula fechada.")


if __name__ == "__main__":
    main()

# Nota sobre a variância, que é o detalhe que Girsanov não apaga.
#
# A reponderação é exata e o estimador é não viesado nos quatro mundos. O que
# ela não preserva é a precisão, e a coluna de erro padrão mostra isso.
#
# O padrão observado não é "quanto mais longe do mundo neutro ao risco, pior".
# A linha de mu = -5% tem o maior erro padrão e a de mu = 25% tem o menor,
# embora seja ela a mais distante em lambda. A razão é a interação entre o
# reponderador e o payoff. Com drift baixo, pouquíssimas trajetórias terminam
# acima do strike, e essas poucas carregam peso alto, então a média se apoia em
# um punhado de pontos. Com drift alto, muitas trajetórias terminam no dinheiro
# e o payoff fica bem amostrado, e a reponderação apenas corrige o excesso.
#
# Esse é exatamente o princípio da amostragem por importância, que usa a
# liberdade dada por Girsanov na direção contrária. Em vez de sofrer a troca de
# medida, escolhe-se de propósito simular sob um drift que põe massa onde o
# payoff é diferente de zero, e depois se corrige com dQ/dP. O preço sai o
# mesmo, com menos trajetórias. Trocar a medida não custa acurácia, custa ou
# devolve precisão, dependendo da direção da troca.
