# Stochastic Calculus Labs

Simuladores interativos de cálculo estocástico e os scripts que verificam
numericamente cada resultado, em português.

O material nasceu de uma trilha de matemática para entrevistas de quant
trading. Os simuladores rodam no navegador, sem instalação e sem cadastro, e os
scripts de verificação existem porque nenhuma afirmação matemática das páginas
foi publicada sem ser conferida por simulação antes.

## Simuladores

Arquivos HTML autocontidos, sem dependência externa. Basta abrir no navegador.

| Lab | O que faz | Página |
|---|---|---|
| [`lab-movimento-browniano.html`](lab-movimento-browniano.html) | Trajetórias de movimento browniano, variação quadrática contra variação total | [Lema de Itô](https://quanttrainer.com.br/lema-de-ito.html) |
| [`lab-equacao-do-calor.html`](lab-equacao-do-calor.html) | Difusão por diferenças finitas precificando uma opção, com Black-Scholes sobreposto | [Black-Scholes e a equação do calor](https://quanttrainer.com.br/black-scholes-equacao-do-calor.html) |

## Verificações numéricas

```bash
pip install -r requirements.txt
python verificacoes/ito_variacao_quadratica.py
python verificacoes/black_scholes_equacao_do_calor.py
python verificacoes/girsanov_mudanca_de_medida.py
```

Cada script é executável sozinho, imprime uma tabela e termina com asserções.
Se uma asserção falhar, alguma afirmação da página correspondente está errada.

**[`ito_variacao_quadratica.py`](verificacoes/ito_variacao_quadratica.py)**
Mostra a variação quadrática convergindo para `T` enquanto a variação total
diverge como `sqrt(n)`, que é a origem do termo de segunda ordem. Depois separa
média de mediana no movimento browniano geométrico e mede o drift da
log-trajetória, que dá `mu - sigma^2/2` e não `mu`, com a hipótese errada
rejeitada por cerca de 300 erros padrão.

**[`black_scholes_equacao_do_calor.py`](verificacoes/black_scholes_equacao_do_calor.py)**
Resolve `u_tau = u_xx`, a equação do calor sem nenhum parâmetro, por
Crank-Nicolson, desfaz a mudança de variáveis e compara com a fórmula fechada.
Varre taxas de 0% a 15%, volatilidades de 20% a 45% e prazos de seis meses a um
ano. Maior desvio observado de 0,0016% do strike, que é erro de discretização.

**[`girsanov_mudanca_de_medida.py`](verificacoes/girsanov_mudanca_de_medida.py)**
Simula sob quatro drifts reais diferentes, de −5% a 25%, reponderá cada
trajetória pela derivada de Radon-Nikodym e mostra os quatro preços caindo
sobre o mesmo valor de Black-Scholes. O script também reporta o erro padrão de
cada linha, que não é igual entre elas, e o comentário no fim do arquivo explica
por quê. Esse detalhe é o princípio da amostragem por importância.

## Por que em português

Quase todo o material bom de cálculo estocástico está em inglês, e uma parte
grande do que existe em português é tradução automática ou resumo sem dedução.
As páginas fazem a conta inteira, passo a passo, com as passagens que os livros
costumam pular.

## Licença

MIT. Use como quiser, inclusive em aula.
