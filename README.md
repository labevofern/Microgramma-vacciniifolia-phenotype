
dados_lineares <- readxl::read_excel(
  here::here("data", "linear", "Dados lineares.xlsx")
)
```

Para tornar o fluxo totalmente automatizado, os quatro scripts podem ser posteriormente divididos em funções e coordenados por um arquivo principal ou pelo pacote [`targets`](https://books.ropensci.org/targets/).

## Arquivos intermediários

Os scripts usam nomes diferentes para algumas matrizes. Antes das análises integrativas, padronize ou copie os arquivos da seguinte forma:

| Arquivo produzido | Nome esperado na integração |
|---|---|
| `Dados lineares_Matriz_Distancia_Euclidiana.csv` | `dist_linear.csv` |
| `Dados outline_Matriz_Distancia_Euclidiana.csv` | `dist_outline.csv` |
| `dist_ambiental.csv` | `dist_ambiental.csv` |

Confirme que todas as matrizes são quadradas, simétricas e possuem os mesmos identificadores de espécimes.

## Principais saídas

### Morfometria linear

- diagnóstico da imputação por `missForest`;
- resultados de Kruskal–Wallis e Dunn;
- gráficos de violino, PCA e LDA;
- matriz de confusão e acurácia da LDA;
- matriz e heatmap de distâncias euclidianas;
- relações morfologia–ambiente com correção FDR de Benjamini–Hochberg;
- figura integrada com gráfico de efeitos, matriz de relações e rede;
- tabela e matriz gráfica de alometria.

### Morfometria de contornos

- visualização dos contornos e diagnóstico do número de harmônicos;
- scores e gráficos de PCA e LDA;
- formas médias por domínio;
- PERMANOVA dos coeficientes de Fourier;
- matriz e heatmap de distâncias;
- mapa com exemplares representativos.

### Ambiente e distribuição espacial

- tabela de ocorrências com variáveis ambientais;
- seleção de variáveis por VIF;
- scores e gráfico de NMDS;
- PERMANOVA e testes de dispersão multivariada;
- perfis de silhueta e agrupamentos PAM;
- testes de Mantel;
- mapas e tabelas de frequência dos agrupamentos.

### Análises integrativas

- gráficos de dbRDA para medidas lineares e contornos;
- testes globais e por variável ambiental;
- comparação de Procrustes;
- figura consolidada para publicação.

## Parâmetros analíticos importantes

- semente aleatória: `123` nas principais análises estocásticas;
- remoção de outliers: regra de `1,5 × IQR` dentro de cada variável e domínio;
- ajuste do teste de Dunn: Bonferroni;
- seleção ambiental: `VIF < 7`;
- exclusão inicial de variáveis ambientais: mais de 15% de valores ausentes;
- PERMANOVA, betadisper, dbRDA e Mantel: `9.999` permutações nos blocos principais;
- Procrustes: `999` permutações;
- dbRDA: variáveis com `p < 0,1` são mantidas no modelo reduzido;
- GLMs: variáveis padronizadas, distribuição gaussiana e correção FDR de Benjamini–Hochberg;
- Fourier elíptico: oito harmônicos no bloco atual;
- agrupamento PAM: `k = 2` no bloco atual, após avaliação por silhueta;
- NMDS ambiental: cinco dimensões e distância euclidiana.

Esses parâmetros devem ser revisados se o conjunto de dados ou a pergunta biológica forem alterados.

## Ajustes necessários antes de uma execução completa

O código contém alguns pontos dependentes da sessão de trabalho original:

1. **Caminhos absolutos:** todos os `setwd()` devem ser convertidos em caminhos relativos.
2. **Nomes das matrizes:** padronize os nomes descritos na seção “Arquivos intermediários”.
3. **Identificador `Herbario`:** remova espaços extras e confira duplicatas antes das junções.
4. **`Linear_.R`:** as ocorrências de `ggreRHD` e do argumento `reRHD` parecem ser erros de digitação; o pacote e o argumento usados nos demais blocos são `ggrepel` e `repel`.
5. **`Outline_.R`:** o objeto `lf_fou1` utilizado na chamada da PCA não é criado anteriormente. Confirme se o objeto correto é `lf_fou` ou se o fator deve ser obtido de outro objeto.
6. **`Ambiental_.R`:** o renomeio das variáveis pressupõe exatamente 39 camadas depois do filtro de valores ausentes. Verifique a quantidade e a ordem das camadas antes de renomeá-las.
7. **Colunas por posição:** alguns blocos selecionam variáveis por número de coluna. Prefira nomes explícitos para evitar mudanças silenciosas quando a planilha for alterada.
8. **Buffer espacial:** valide a unidade de `width` depois de definir ou reprojetar o CRS.
9. **Execução por blocos:** objetos como `df_media_final`, `lf_fou_avg` e as matrizes de distância devem existir antes dos blocos integrativos.

## Reprodutibilidade

Para registrar as versões dos pacotes utilizadas em uma análise concluída:

```r
sessionInfo()
```

Para maior reprodutibilidade, recomenda-se inicializar o projeto com [`renv`](https://rstudio.github.io/renv/):

```r
install.packages("renv")
renv::init()
renv::snapshot()
```

## Dados e privacidade

Antes de publicar os dados, confirme se coordenadas de coleta, códigos de herbário e informações sobre espécimes podem ser disponibilizados publicamente. Caso necessário, deposite os dados em um repositório apropriado e mantenha no GitHub apenas metadados, scripts e instruções de acesso.

## Citação

Caso utilize ou adapte este fluxo, cite o repositório e a publicação associada. Um modelo pode ser preenchido quando a URL e o ano da versão pública estiverem definidos:

```text
Aljahara J (ANO). Análises integradas de macromorfologia foliar e ambiente.
Repositório GitHub. URL: https://github.com/USUARIO/REPOSITORIO
```
