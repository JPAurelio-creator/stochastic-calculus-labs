# Stochastic Calculus Labs

Interactive stochastic calculus simulators, plus the scripts that check every
result numerically. Everything runs in the browser, with no install and no account.

[English](#english) · [Português](#português)

---

## English

### Simulators

Self-contained HTML files with no external dependencies. Open one in any
browser, or run it live on the site.

| Lab | What it does | Run it live | Source |
|---|---|---|---|
| Heat equation | Solves the Black-Scholes PDE by finite differences in your browser and overlays the closed-form price. Pick a payoff (call, put, digital, straddle, butterfly, or draw your own), drag time to expiry, and read price, delta and gamma. | [heat-equation-lab](https://quanttrainer.com.br/heat-equation-lab.html) | [`heat-equation-lab.html`](heat-equation-lab.html) |
| Brownian motion | Draws Brownian sample paths inside the ±2√t envelope, to show the square-root-of-time scaling behind every volatility calculation. | [brownian-motion-lab](https://quanttrainer.com.br/brownian-motion-lab.html) | [`brownian-motion-lab.html`](brownian-motion-lab.html) |

### Numerical checks

```bash
pip install -r requirements.txt
python verificacoes/ito_variacao_quadratica.py
python verificacoes/black_scholes_equacao_do_calor.py
python verificacoes/girsanov_mudanca_de_medida.py
```

Each script runs on its own, prints a table and ends with assertions. If an
assertion fails, a claim on the matching page is wrong. File names and code
comments are in Portuguese, and the numbers read the same in any language.

- **`ito_variacao_quadratica.py`** shows quadratic variation converging to `T`
  while total variation diverges like `sqrt(n)`, which is where the second-order
  term in Itô's lemma comes from. It then measures the drift of the log-path of
  geometric Brownian motion, which comes out as `mu - sigma^2/2` and not `mu`,
  rejecting the wrong hypothesis by about 300 standard errors.
- **`black_scholes_equacao_do_calor.py`** solves the parameter-free heat equation
  `u_tau = u_xx` with Crank-Nicolson, undoes the change of variables and compares
  with the closed form across rates from 0% to 15%, vols from 20% to 45% and
  maturities from six months to a year. Largest deviation 0.0016% of strike.
- **`girsanov_mudanca_de_medida.py`** simulates under four different real-world
  drifts, reweights each path by the Radon-Nikodym derivative and shows all four
  prices landing on the same Black-Scholes value. The standard errors differ
  between rows, and the comment at the end of the file explains why. That detail
  is the principle behind importance sampling.

### Formula sheet and derivations

A one-page [stochastic calculus formula sheet](https://quanttrainer.com.br/formulario-calculo-estocastico.html)
covers Itô's lemma, the product rule, quadratic variation, GBM, Ornstein-Uhlenbeck,
CIR, the Brownian bridge, Girsanov, Feynman-Kac, Black-Scholes and the Greeks. Each
entry lists the assumptions it needs and the case where it breaks. The sheet and the
full derivations are written in Portuguese.

---

## Português

Simuladores interativos de cálculo estocástico e os scripts que conferem cada
resultado numericamente. Os dois laboratórios existem em português,
[`lab-equacao-do-calor.html`](lab-equacao-do-calor.html) e
[`lab-movimento-browniano.html`](lab-movimento-browniano.html), e rodam ao vivo no
[site](https://quanttrainer.com.br/lab-equacao-do-calor.html).

Os scripts de verificação existem porque nenhuma afirmação matemática das páginas
foi publicada sem conferência numérica antes. As deduções completas estão nos
artigos sobre o [lema de Itô](https://quanttrainer.com.br/lema-de-ito.html),
[Black-Scholes como equação do calor](https://quanttrainer.com.br/black-scholes-equacao-do-calor.html),
[Girsanov](https://quanttrainer.com.br/teorema-de-girsanov.html) e
[Feynman-Kac](https://quanttrainer.com.br/feynman-kac.html), e as fórmulas estão
reunidas no [formulário de cálculo estocástico](https://quanttrainer.com.br/formulario-calculo-estocastico.html).

---

## License

MIT. Use it however you like, including in class.
