# ==============================================================
# Análise morfométrica de medidas lineares
# ==============================================================
# Autor: Juliana Aljahara
# Programa: Doutoranda em Biologia Vegetal
# Instituição: Universidade Federal de Pernambuco (UFPE)
# Orientadora: Dra. Thaís E. Almeida
# ==============================================================================
# 1. PACOTES E PREPARAÇÃO
# ==============================================================================

library(readxl)        # Importação de planilhas do Excel (.xls, .xlsx) direto para o R
library(ggplot2)       # Criação de gráficos estáticos de alta qualidade baseados na gramática de gráficos
library(MASS)          # Funções estatísticas avançadas e suporte para modelos lineares e distribuições
library(dplyr)         # Manipulação, filtragem e organização de tabelas de dados (data.frames)
library(tidyr)         # Organização e reestruturação de dados (converter entre formatos largo e longo)
library(vegan)         # Análises de ecologia de comunidades (diversidade, ordenação e estatística multivariada)
library(factoextra)    # Visualização e extração de resultados de análises multivariadas (como PCA e Cluster)
library(corrplot)      # Visualização gráfica de matrizes de correlação
library(missForest)    # Imputação de dados faltantes (NA) usando algoritmos de Random Forest
library(purrr)         # Programação funcional para aplicar funções a listas ou colunas de forma eficiente
library(FSA)           # Ferramentas específicas para análise de dados de fisheries (ciências pesqueiras)
library(ggpubr)        # Facilita a criação de gráficos prontos para publicação integrados ao ggplot2
library(pheatmap)      # Criação de mapas de calor (heatmaps) altamente customizáveis e com agrupamento automático
library(RColorBrewer)  # Fornece paletas de cores prontas e perceptualmente uniformes para mapas e gráficos
library(ggbeeswarm)    # Criação de gráficos de dispersão (stripcharts) que evitam a sobreposição de pontos
library(tidyverse)     # Coleção de pacotes (dplyr, ggplot2, tidyr, etc.) para ciência de dados e fluxo de trabalho integrado
library(rlang)         # Fornece ferramentas para programação genérica e manipulação de expressões e ambientes no R
library(vegan)         # Para a PERMANOVA (adonis2)
library(ggbeeswarm)    # Para geom_quasirandom
library(ggrepel)       # Para formatações de texto/rótulos se necessário
library(rstatix)       # Para testes de Kruskal-Wallis e Dunn estruturados
library(writexl)       # Para salvar arquivos em formato .xlsx

# Cores do artigo
my_colors <- c("#6A0572", "#AB83A1", "#F67280", "#F8B195")

# Diretório base
setwd("C:/Users/nikso/OneDrive/_JULIANA ALJAHARA/1_Análise de macromorfologica/Diretório_Linear")

dir.create("2_Tabelas_e_Imagens", showWarnings = FALSE)

# ==============================================================================
# 2. LEITURA E LIMPEZA (Sincronizado com o print do seu Excel)
# ==============================================================================

dados <- read_excel("2_Tabelas_e_Imagens/Dados lineares.xlsx")

# Renomeando as colunas conforme a ordem do seu print:
# 1:Indivíduos, 2:Numero, 3:Herbário, 4:Domínio, 5:Lat, 6:Log, 7-13:Medidas
df_base <- dados %>%
  rename(
    id = 1,          # Mudamos para 'id' para manter compatibilidade com o pivot_longer
    numero = 2,
    Herbario = 3,
    `Dominio fitogeografico`= 4,       # Domínio fitogeográfico
    Lat = 5,
    Log = 6,
    AAN = 7, BAN = 8, LBL = 9, TOL = 10, PEL = 11, PED = 12, RHD = 13, LFA = 14
  )

# Converter medidas para numérico (vírgula por ponto) e bioma para fator
df_base <- df_base %>%
  mutate(across(7:13, ~as.numeric(gsub(",", ".", .)))) %>%
  mutate(`Dominio fitogeografico` = as.factor(`Dominio fitogeografico`))

# Preparar dados para imputação
medidas_com_na <- df_base %>% dplyr::select(AAN, BAN, LBL, TOL, PEL, PED, RHD, LFA)

cat("\n=== Rodando Imputação de NAs ===\n")
set.seed(123)
imputacao <- missForest(as.matrix(medidas_com_na))
print(imputacao$OOBerror)

# Criar a tabela de performance global
tabela_performance <- data.frame(
  Metrica = c("Algoritmo", "Erro de Imputação (NRMSE)", "Variáveis Analisadas", "Total de Observações"),
  Valor = c(
    "missForest (Random Forest)", 
    round(imputacao$OOBerror, 6), 
    paste(ncol(medidas_com_na), "variáveis"),
    nrow(df_base)
  )
)

# Salvar
write.csv(tabela_performance, "Tabela_Suplementar_Qualidade_Imputacao_linear.csv", row.names = FALSE)

print(tabela_performance)

# CRIAR DF_LIMPO (Aqui garantimos que a coluna 'id' e 'bioma' existam)
df_limpo <- data.frame(
  id = df_base$id,
  bioma = df_base$`Dominio fitogeografico`,
  Herbario = df_base$Herbario,
  Lat = df_base$Lat,
  Log = df_base$Log,
  imputacao$ximp
)

df_media_final <- df_limpo %>%
  # 1. Agrupamos APENAS RHDos identificadores principais
  dplyr::group_by(id, Herbario) %>%
  
  # 2. Fazemos o resumo de tudo de uma vez só
  dplyr::summarise(
    # Mantém o primeiro valor encontrado para as colunas categóricas/coordenadas
    bioma = first(bioma),
    Herbario = first(Herbario),
    Lat = mean(Lat, na.rm = TRUE),
    Log = mean(Log, na.rm = TRUE),
    
    # Calcula a média para as colunas morfológicas/numéricas da imputação
    across(c(AAN, BAN, LBL, TOL, PEL, PED, RHD, LFA), mean, na.rm = TRUE),
    
    .groups = "drop"
  )

# ==============================================================================
# 3. PREPARAÇÃO DOS DADOS E REMOÇÃO DE OUTLIERS
# ==============================================================================

df_processado <- df_media_final %>%
  pivot_longer(
    cols = -c(id, bioma, Herbario, Lat, Log), 
    names_to = "Variavel", 
    values_to = "Valor"
  ) %>%
  group_by(Variavel, bioma) %>%
  mutate(
    Q1 = quantile(Valor, 0.25, na.rm = TRUE),
    Q3 = quantile(Valor, 0.75, na.rm = TRUE),
    IQR = Q3 - Q1,
    limite_inf = Q1 - 1.5 * IQR,
    limite_sup = Q3 + 1.5 * IQR,
    is_outlier = Valor < limite_inf | Valor > limite_sup
  ) %>%
  # Mantém apenas as observações que NÃO são outliers
  filter(!is_outlier) %>% 
  ungroup() %>%
  # REMOVE AS COLUNAS REPETIDAS QUE JÁ CUMPRIRAM SEU PARHD:
  select(-c(Q1, Q3, IQR, limite_inf, limite_sup, is_outlier))

# ------------------------------------------------------------------------------
# ANÁLISE ESTATÍSTICA (DADOS LIMPOS) - INCLUINDO TODOS OS PARES
# ------------------------------------------------------------------------------

# 2.1. Teste Global de Kruskal-Wallis
res_kruskal <- df_processado %>%
  group_by(Variavel) %>%
  kruskal_test(Valor ~ bioma)

# 2.2. Teste Post-hoc de Dunn com ajuste de Bonferroni 
# Calcula as coordenadas para ABSOLUTAMENTE TODAS as combinações de biomas
res_dunn <- df_processado %>%
  group_by(Variavel) %>%
  dunn_test(Valor ~ bioma, p.adjust.method = "bonferroni") %>%
  add_xy_position(x = "bioma", formula = Valor ~ bioma, scales = "free_y")

# [ALTERADO]: Não filtramos mais por p.adj < 0.05. Vamos usar o 'res_dunn' completo.

# Visualização completa dos Resultados no Console
print("--- Resultado Kruskal-Wallis ---")
print(res_kruskal)

print("--- Resultado Post-hoc de Dunn (Todos os pares) ---")
print(res_dunn)

# ------------------------------------------------------------------------------
# GRÁFICO DE VIOLINO / BOXPLOT COM BARRAS DE SIGNIFICÂNCIA PAREADA
# ------------------------------------------------------------------------------

p_visibilidade_total <- df_processado %>%
  ggplot(aes(x = bioma, y = Valor, fill = bioma)) +
  
  # Violino (Refletindo a distribuição sem os extremos)
  geom_violin(alpha = 0.2, trim = FALSE, color = "gray80") +
  
  # Pontos (Beeswarm) - Distribuição interna dos dados
  geom_quasirandom(aes(color = bioma), size = 1.2, alpha = 0.3, width = 0.25) +
  
  # Boxplot Interno (outliers ocultados pois já foram removidos do DF)
  geom_boxplot(width = 0.1, fill = "white", color = "black", alpha = 0.7, 
               outlier.shape = NA) + 
  
  # Adiciona o P-valor Global (Kruskal-Wallis) fixado no canto superior esquerdo
  stat_compare_means(
    method = "kruskal.test", 
    aes(label = paste0("italic(p)[Kruskal] == ", ..p.format..)),
    parse = TRUE, label.y.npc = 0.99, label.x.npc = 0.05, hjust = 0, size = 3.8
  ) + 
  
  # ADICIONA AS BARRAS PAR A PAR DO POST-HOC DE DUNN
  # Modifique para label = "p.adj.format" se preferir o valor numérico (ex: 0.01) em vez de asteriscos (*).
  stat_pvalue_manual(
    res_dunn, 
    label = "p.adj.signif", 
    tip.length = 0.01, 
    step.increase = 0.06, 
    hide.ns = TRUE
  ) +
  
  # Facetamento por Variável com eixos livres no eixo Y
  facet_wrap(~Variavel, scales = "free_y") +
  
  # Customização Estética
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  theme_minimal(base_size = 14) +
  labs(x = "Phytogeographic domain", y = "Measured Value") +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    panel.grid.major.x = element_blank(),
    panel.spacing = unit(1.5, "lines")
  )

# Exibir o gráfico gerado
print(p_visibilidade_total)

# ------------------------------------------------------------------------------
# SALVAR ARQUIVOS E IMAGENS
# ------------------------------------------------------------------------------

# Salvar Imagem do Gráfico Finalizado
ggsave("2_Tabelas_e_Imagens/Violin_Visibilidade_Dunn_Signif.png", 
       plot = p_visibilidade_total, width = 14, height = 10, dpi = 300)

# Limpar as colunas de coordenadas gráficas do Dunn antes de exportar as tabelas
res_dunn_tabela <- res_dunn %>% 
  select(-c(.y., groups, y.position, x, xmin, xmax))

# Salvar tabelas dos resultados estatísticos
write.csv(res_kruskal, "2_Tabelas_e_Imagens/Resultado_Kruskal_Wallis_Sem_Outliers.csv", row.names = FALSE)
write.csv(res_dunn, "2_Tabelas_e_Imagens/Resultados_Dunn_PostHoc.csv", row.names = FALSE)
write_xlsx(res_dunn, path = "2_Tabelas_e_Imagens/Resultados_Dunn_PostHoc.xlsx")

# ==============================================================================
# 4. PCA - ANÁLISE DE COMPONENTES PRINCIPAIS
# ==============================================================================

library(factoextra)
library(dplyr)

# 1. Isolar as variáveis numéricas a partir da coluna 6
medidas <- df_limpo %>% 
  select(6:last_col())

# 2. Isolar o vetor de biomas do mesmo dataframe para mapear as cores
bioma_fatores <- df_limpo$bioma

# 3. Executar a PCA (com os dados padronizados)
pca <- prcomp(medidas, scale. = TRUE)

p_pca <- fviz_pca_biplot(pca,
                         geom.ind = "point",
                         pointshape = 21,
                         pointsize = 5,
                         fill.ind = bioma_fatores,
                         col.ind = "black",
                         col.var = "black",
                         reRHD = TRUE,
                         legend.title = "",
                         title = "PCA Biplot - Linear Measurements",
                         palette = my_colors) +
  theme_classic()

print(p_pca)

ggsave("2_Tabelas_e_Imagens/PCA_Biplot.png",
       plot = p_pca,
       width = 10,
       height = 8,
       dpi = 300)

# ==============================================================================
# 5. PERMANOVA
# ==============================================================================

set.seed(123)

# Seleciona da coluna 6 até a última coluna do dataframe
medidas_selecionadas <- df_media_final[, 6:ncol(df_media_final)]

# Executa o PERMANOVA aplicando o scale() apenas nas colunas corretas
permanova_res <- adonis2(scale(medidas_selecionadas) ~ bioma,
                         data = df_media_final,
                         method = "euclidean",
                         permutations = 9999)

print("=== PERMANOVA Results ===")
print(permanova_res)

# ==============================================================================
# 6. LDA - ANÁLISE DISCRIMINANTE LINEAR
# ==============================================================================

# 1. Calculando as médias das colunas numéricas agrupadas por Herbario
# Mantemos o 'bioma' no agrupamento para não perdê-lo na sumarização
df_med <- df_media_final %>%
  group_by(Herbario, bioma) %>%
  summarise(
    across(c(AAN, BAN, PEL, TOL, RHD, PED, RHD), mean, na.rm = TRUE),
    .groups = "drop"
  )

# 2. Preparando os dados para a LDA
# A variável resposta (o que queremos prever)
bioma2 <- df_med$bioma

# Selecionamos apenas as colunas numéricas (os códigos informados) para a análise
medidas2 <- df_med %>% 
  dplyr::select(AAN, BAN, PEL, TOL, RHD, PED, RHD)

# 3. Executando a LDA com Validação Cruzada (CV = TRUE)
lda_cv <- lda(bioma2 ~ ., data = data.frame(bioma2, medidas2), CV = TRUE)

# 4. Resultados
pred <- lda_cv$class

cat("\n=== Confusion Matrix (Cross-Validation) ===\n")
tab <- table(Real = bioma2, Predicted = pred)
print(tab)

acc_geral <- sum(diag(tab)) / sum(tab)

cat("\nOverall Accuracy:",
    round(acc_geral * 100, 2), "%\n")

# ==============================================================================
# 7. PLOT LDA (COM ESTÉTICA APLICADA)
# ==============================================================================

# 1. Preparação do modelo e predição
lda_model <- lda(bioma2 ~ ., data = data.frame(bioma2, medidas2))
lda_pred <- predict(lda_model)

# 2. Criar o rótulo da acurácia
label_acc <- paste0("Accuracy: ", round(acc_geral * 100, 2), "%")

if (ncol(lda_pred$x) >= 2) {
  
  lda_df <- data.frame(LD1 = lda_pred$x[,1],
                       LD2 = lda_pred$x[,2],
                       bioma = bioma2)
  
  p_lda <- ggplot(lda_df, aes(LD1, LD2, fill = bioma, color = bioma)) +
    # Estética dos pontos: shape 21 permite preenchimento (fill) e borda (color)
    geom_point(shape = 21, size = 5, color = "black", alpha = 0.8) +
    stat_ellipse(show.legend = FALSE, linetype = 2) +
    scale_fill_manual(values = my_colors) +
    scale_color_manual(values = my_colors) +
    # Configurações de dimensão de fonte e tema limpo
    theme_classic(base_size = 14) +
    labs(title = "",
         subtitle = label_acc,
         x = "LD1",
         y = "LD2",
         fill = "Bioma",
         color = "Bioma") +
    # Negrito nos títulos dos eixos e remoção de grades
    theme(legend.position = "bottom",
          panel.grid = element_blank(),
          axis.title = element_text(face = "bold"))
  
} else {
  
  lda_df <- data.frame(LD1 = lda_pred$x[,1], bioma = bioma2)
  
  p_lda <- ggplot(lda_df, aes(x = LD1, fill = bioma)) +
    geom_density(alpha = 0.7, color = "black") +
    scale_fill_manual(values = my_colors) +
    theme_classic(base_size = 14) +
    labs(title = "",
         subtitle = label_acc,
         x = "LD1",
         y = "Densidade",
         fill = "Bioma") +
    theme(legend.position = "bottom",
          axis.title = element_text(face = "bold"))
}

# Exibir
print(p_lda)

# Salvar com as dimensões específicas (10x8 pol, 300 DPI)
ggsave("2_Tabelas_e_Imagens/LDA_Plot_Customizado.png",
       plot = p_lda,
       width = 10,
       height = 8,
       dpi = 300)

# ==============================================================================
# EXTRAÇÃO DOS SCORES DA LDA (COM BIOMA E HERBÁRIO)
# ==============================================================================

# 1. Executar o modelo e predição (seu código original)
lda_model <- lda(bioma2 ~ ., data = data.frame(bioma2, medidas2))
lda_pred <- predict(lda_model)

# 2. Criar o DataFrame dos Scores da LDA
df_scores_lda <- as.data.frame(lda_pred$x)

# 3. Adicionar as colunas de identificação (Bioma e Herbário)
df_scores_lda <- df_scores_lda %>%
  mutate(
    Herbario = df_med$Herbario,  # Puxa a coluna diretamente do dataframe sumarizado
    bioma    = bioma2
  ) %>%
  # Coloca o Herbário e o Bioma nas primeiras colunas, seguidos RHDos scores (LDs)
  select(Herbario, bioma, everything()) 

# Exibir as primeiras linhas para conferir
head(df_scores_lda)

# [Opcional] Salvar a tabela na sua pasta de resultados se quiser exportar para o Excel/CSV
write.csv(df_scores_lda, 
          file = "2_Tabelas_e_Imagens/Scores_LDA_Herbario.csv", 
          row.names = FALSE)

# ==============================================================================
# 8. HEATMAP MATRIZ DE CONFUSÃO
# ==============================================================================

conf_df <- as.data.frame(tab)

conf_df <- conf_df %>%
  group_by(Real) %>%
  mutate(Percent = round(Freq / sum(Freq) * 100, 1))

conf_df <- conf_df %>%
  mutate(Label = ifelse(Freq == 0, "",
                        paste0(Percent, "% (", Freq, ")")))

p_conf <- ggplot(conf_df, aes(x = Predicted, y = Real, fill = Percent)) +
  geom_tile(color = "white", size = 0.8) +
  geom_text(aes(label = Label), size = 4, fontface = "bold") +
  scale_fill_gradient(low = "#F8F8F8", high = "#F67280") +
  theme_classic(base_size = 14) +
  labs(title = "",
       x = "Predicted Class",
       y = "Real Class") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))

print(p_conf)

ggsave("2_Tabelas_e_Imagens/Confusion_Matrix.png",
       plot = p_conf,
       width = 8,
       height = 6,
       dpi = 300)

# ==============================================================================
# PROCESSAMENTO DE DADOS E MATRIZ DE DISTÂNCIA (BIOMAS)
# ==============================================================================

# 1. Padronização (Z-score)
# Garante que variáveis com escalas diferentes (ex: RHD vs AAN) tenham o mesmo peso
medidas_scaled <- scale(medidas2)
rownames(medidas_scaled) <- df_med$Herbario

# 2. Cálculo da Matriz de Distância Euclidiana
dist_matrix <- dist(medidas_scaled, method = "euclidean")
dist_table <- as.matrix(dist_matrix)

# 3. Salvamento dos Dados (Tabelas)
write.csv(dist_table, "2_Tabelas_e_Imagens/Dados lineares_Matriz_Distancia_Euclidiana.csv")
write.csv(df_med, "2_Tabelas_e_Imagens/Dados lineares_por variavel_GLM.csv", row.names = FALSE)

# 4. Configuração da Identidade Visual (Biomas)
annotation_data <- data.frame(Bioma = df_med$bioma)
rownames(annotation_data) <- df_med$Herbario

ann_colors = list(
  Bioma = c(
    "Atlantic Forest" = "#6A0572", 
    "Cerrado"         = "#F67280", 
    "Caatinga"        = "#AB83A1", 
    "Pampa"           = "#F8B195"
  )
)

# 5. Geração do Heatmap de Alta Definição
# O pheatmap com 'filename' ajusta a imagem automaticamente ao tamanho do paRHD
pheatmap(dist_table, 
         clustering_method = "ward.D2", 
         annotation_row = annotation_data, 
         annotation_col = annotation_data, 
         annotation_colors = ann_colors,
         # Roxo: Branco (Baixa distância/Muita relação) -> Escuro (Alta distância/Pouca relação)
         color = colorRampPalette(brewer.pal(9, "Purples"))(100),
         
         # --- Estética e Simetria ---
         border_color = "white",
         lwd = 0.5,
         main = "", 
         
         # --- Textos e Legendas ---
         fontsize_row = 7, 
         fontsize_col = 7,
         angle_col = 45,
         annotation_names_row = FALSE, 
         annotation_names_col = FALSE,
         
         # --- Agrupamento e Quebras ---
         cutree_rows = 4, 
         cutree_cols = 4,
         treeheight_row = 60,
         treeheight_col = 60,
         
         # --- Configuração do Arquivo Final ---
         filename = "2_Tabelas_e_Imagens/Heatmap_Distancia_Final_Ajustado.png",
         width = 15,   # Tamanho ideal para manter a proporção
         height = 10,  
         res = 300
)

# ==========================================================
# 2. DIRETÓRIO E DADOS
# ==========================================================

setwd("C:/Users/nikso/OneDrive/_JULIANA ALJAHARA/1_Análise de macromorfologica/Diretório_Integrativa/2_Tabelas_e_Imagens")

# Matrizes de distância
dist_linear  <- read.csv("dist_linear.csv", row.names = 1)
dist_outline <- read.csv("dist_outline.csv", row.names = 1)

# Variáveis ambientais
env <- read.csv("Especimes_e_Variaveis_Selecionadas_VIF.csv", row.names = 1)

# ==============================================================================
# 9. GLM COM PADRONIZAÇÃO E UNIÃO SEGURA
# ==============================================================================

# Carregar o pacote necessário para manipulação de dados
library(dplyr)

# 2. Preparação e Limpeza rigorosa dos nomes para garantir o Merge
# (Substitua 'df_media_final' RHDo nome do seu dataframe atual, caso mudou)
df_final <- df_media_final 

# Remove espaços em branco extras que podem causar erro na união
df_final$Herbario <- trimws(as.character(df_final$Herbario))
env$Herbario      <- trimws(as.character(env$Herbario))

# 3. União dos dados (Merge seguro)
df_glm_completo <- left_join(df_final, env, by = "Herbario")

# --- CHECKPOINT DE ERRO ---
# Checa se as colunas ambientais vieram vazias (indica erro no nome do herbário)
# Ajustei para testar diretamente uma coluna que veio do arquivo 'env'
var_teste <- colnames(env)[2] # Pega a primeira variável ambiental real do seu arquivo
n_nas <- sum(is.na(df_glm_completo[[var_teste]])) 

if(n_nas > 0) {
  warning("Atenção: ", n_nas, " linhas não encontraram correspondência ambiental. Verifique os nomes dos Herbários!")
}

# ==========================================================
# 4. Definir listas de colunas para o GLM
# ==========================================================
# IMPORTANTE: Verifique se os números das colunas (índices) continuam os mesmos 
# no seu dataframe final unificado (df_glm_completo)
medidas_lista    <- colnames(df_glm_completo)[6:13]
ambientais_lista <- colnames(df_glm_completo)[22:41]

# ==========================================================
# 2. Função GLM
# ==========================================================

rodar_glm_ambiental <- function(medida, ambiental) {
  
  temp_df <- df_glm_completo %>%
    dplyr::select(all_of(c(medida, ambiental))) %>%
    # Garante que ambos sejam numéricos antes do drop_na
    mutate(across(everything(), ~as.numeric(as.character(.)))) %>% 
    tidyr::drop_na()
  
  if(nrow(temp_df) < 5) return(NULL)
  
  # O scale() agora não deve falhar
  temp_df[[medida]] <- scale(temp_df[[medida]])
  temp_df[[ambiental]] <- scale(temp_df[[ambiental]])
  
  modelo <- tryCatch(
    glm(as.formula(paste(medida, "~", ambiental)),
        data = temp_df,
        family = gaussian()),
    error = function(e) return(NULL)
  )
  
  if(!is.null(modelo)) {
    coefs <- summary(modelo)$coefficients
    
    if(nrow(coefs) > 1) {
      return(data.frame(
        Medida = medida,
        Variavel_Ambiental = ambiental,
        Estimativa = coefs[2, 1],
        Erro_Padrao = coefs[2, 2],
        p_value = coefs[2, 4]
      ))
    }
  }
}

# ==========================================================
# 3. Rodar todos os modelos
# ==========================================================

grid_testes <- expand.grid(
  medida = medidas_lista,
  ambiental = ambientais_lista
)

resultado_final_glm <- map2_df(
  grid_testes$medida,
  grid_testes$ambiental,
  rodar_glm_ambiental
)

# ==========================================================
# 4. Controle de FDR (Benjamini-Hochberg) e classificação de significância
# ==========================================================

resultado_final_glm <- resultado_final_glm %>%
  mutate(
    p_FDR = p.adjust(p_value, method = "BH"),
    Significancia = case_when(
      p_FDR < 0.001 ~ "***",
      p_FDR < 0.01  ~ "**",
      p_FDR < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )

# ==========================================================
# 5. Preparar dados para TODOS os gráficos (SEM WARNING)
# ==========================================================

dados_plot <- purrr::pmap_df(
  list(medidas_lista %>% rep(each = length(ambientais_lista)),
       ambientais_lista %>% rep(times = length(medidas_lista))),
  
  function(med, amb) {
    
    temp <- df_glm_completo %>%
      dplyr::select(all_of(c(med, amb))) %>%
      # ADICIONE ESTA LINHA ABAIXO:
      mutate(across(everything(), ~as.numeric(as.character(.)))) %>% 
      tidyr::drop_na()
    
    if(nrow(temp) < 5) return(NULL)
    
    data.frame(
      x = scale(temp[[amb]])[,1],
      y = scale(temp[[med]])[,1],
      Medida = med,
      Variavel_Ambiental = amb
    )
  }
)

# ==========================================================
# 6. Juntar p-values
# ==========================================================

dados_plot <- dados_plot %>%
  left_join(resultado_final_glm,
            by = c("Medida", "Variavel_Ambiental"))

# ==========================================================
# 7. Criar gráfico FINAL
# ==========================================================

# 1. Preparar os dados (mantendo a lógica anterior)
resultado_final_glm <- resultado_final_glm %>%
  mutate(
    sig_color = ifelse(p_FDR < 0.05, "Significant", "Non-significant"),
    p_label = paste0("italic(p)[FDR] == ", signif(p_FDR, 2))
  )

# 2. Criar o gráfico com fontes maiores
p <- ggplot(dados_plot, aes(x = x, y = y)) +
  geom_point(size = 1.8, alpha = 0.7) +
  
  geom_smooth(method = "glm", 
              method.args = list(family = "gaussian"), 
              color = "red", 
              se = TRUE, 
              linewidth = 0.6) +
  
  facet_grid(Medida ~ Variavel_Ambiental, scales = "free") +
  
  geom_text(
    data = resultado_final_glm,
    aes(label = p_label, color = sig_color),
    x = Inf, y = Inf,
    hjust = 1.1, vjust = 1.5,   # Aumentei um pouco o vjust para não encostar no topo com a fonte maior
    size = 3.5,                 # AUMENTADO: de 2.5 para 3.5 (p-valores)
    parse = TRUE,
    inherit.aes = FALSE
  ) +
  
  scale_color_manual(values = c("Significant" = "blue", "Non-significant" = "black")) +
  
  labs(
    x = "Environmental variable",
    y = "Measurement"
  ) +
  
  theme_bw(base_size = 14) +    # AUMENTADO: base_size maior ajuda no escalonamento geral
  
  theme(
    # Títulos das facetas (os nomes das variáveis em cima e ao lado)
    strip.text = element_text(size = 7.4, face = "bold"), # AUMENTADO: de 7 para 10
    
    # Números nos eixos
    axis.text = element_text(size = 9),                  # AUMENTADO: de 6 para 9
    
    # Títulos dos eixos (X e Y)
    axis.title = element_text(size = 11, face = "bold"), # AUMENTADO: de 10 para 12 e negrito
    
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black"),
    strip.background = element_rect(fill = "white"),
    legend.position = "none"
  )

print(p)

# ==========================================================
# 8. Salvar imagem
# ==========================================================

ggsave("GLM_Facets.png",
       p,
       width = 25,   # Aumentado de 18
       height = 18,  # Aumentado de 12
       units = "in",
       dpi = 300,
       scale = 0.8)  # O scale ajuda a "dar zoom" nos elementos se necessário

# ==========================================================
# 6. Juntar e Filtrar apenas Significativas
# ==========================================================

# 1. Filtramos os resultados do GLM para manter apenas relações significativas após FDR (BH)
resultado_sig <- resultado_final_glm %>%
  filter(p_FDR < 0.05) %>%
  mutate(
    sig_color = "blue",
    p_label = paste0("italic(p)[FDR] == ", signif(p_FDR, 2))
  )

# 2. Filtramos o dados_plot usando dplyr::select explicitamente
# O inner_join garante que apenas as relações significativas permaneçam no gráfico
dados_plot_sig <- dados_plot %>%
  inner_join(resultado_sig %>% 
               dplyr::select(Medida, Variavel_Ambiental, p_label, sig_color), 
             by = c("Medida", "Variavel_Ambiental"))

# ==========================================================
# 7. Criar gráfico apenas das SIGNIFICATIVAS
# ==========================================================

p_sig <- ggplot(dados_plot_sig, aes(x = x, y = y)) +
  geom_point(size = 1.8, alpha = 0.7) +
  
  geom_smooth(method = "glm", 
              method.args = list(family = "gaussian"), 
              color = "red", 
              se = TRUE, 
              linewidth = 0.6) +
  
  # A MUDANÇA ESTÁ AQUI:
  # facet_wrap organiza os gráficos em sequência, eliminando os vazios.
  # scales = "free" permite que cada gráfico tenha seu próprio eixo.
  facet_wrap(~ Medida + Variavel_Ambiental, scales = "free", ncol = 4) +
  
  geom_text(
    aes(label = p_label),
    color = "blue",
    x = Inf, y = Inf,
    hjust = 1.1, vjust = 1.5, 
    size = 3.5,
    parse = TRUE
  ) +
  
  labs(
    x = "Environmental variable",
    y = "Measurement",
    title = ""
  ) +
  
  theme_bw(base_size = 14) +
  
  theme(
    # Ajuste o tamanho da fonte do título de cada painel para caber os dois nomes
    strip.text = element_text(size = 8, face = "bold"), 
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black"),
    strip.background = element_rect(fill = "white"),
    legend.position = "none"
  )

print(p_sig)

# Ao salvar, você pode reduzir a altura, já que haverá menos linhas
ggsave("GLM_Facets_signficativas.png",
       p_sig,
       width = 16, 
       height = 12, 
       units = "in",
       dpi = 300)

#=========================================================
# FIGURA INTEGRATIVA (VERSÃO AJUSTADA PARA PUBLICAÇÃO)
#=========================================================

library(ggplot2)
library(ggreRHD)
library(dplyr)
library(patchwork)
library(ggraph)
library(igraph)
library(ggrepel)

#---------------------------
# PALETA
#---------------------------
cor_neg        <- "#16425B" # Azul Marinho Colonial
cor_pos        <- "#D9A05B" # Ouro Velho
cor_node_morph <- "#A8422B" # Vermelho Sangue de Boi
cor_node_env   <- "#CBD4C2" # Cinza Esverdeado Pálido

# Tipografia
t_titulo  <- 18
t_axis    <- 13
t_leg     <- 12
t_labels  <- 4.2

#---------------------------
# DADOS
#---------------------------
dados_fig <- resultado_final_glm %>%
  mutate(
    Direction = ifelse(Estimativa > 0,"Positive","Negative"),
    Weight = abs(Estimativa),
    p_lab = ifelse(p_FDR < 0.001,"<0.001",
                   sprintf("%.3f",p_FDR))
  )

#=========================================================
# PANEL A — VOLCANO COM LEGENDA AJUSTADA
#=========================================================
# Criando a tabela de rótulos filtrando os dados originais
rotulos <- dados_fig %>% 
  filter(p_FDR < 0.05) # Significância após controle de FDR pelo método de Benjamini-Hochberg

p1 <- ggplot(
  dados_fig,
  aes(Estimativa,-log10(p_FDR))
)+
  
  geom_vline(
    xintercept=0,
    linetype=2,
    linewidth=.5,
    alpha=.4
  )+
  
  geom_hline(
    yintercept=-log10(0.05),
    linetype=2,
    linewidth=.5,
    alpha=.4
  )+
  
  geom_point(
    aes(
      fill=Direction,
      size=Weight
    ),
    shape=21,
    color="white",
    stroke=.8,
    alpha=.95
  ) +
  
  # CORREÇÃO AQUI: Mudado de geom_text_reRHD para geom_text_repel
  geom_text_repel(
    data=rotulos,
    aes(label=Medida),
    size=4,
    fontface="bold",
    box.padding=.6,
    point.padding=.5
  ) +
  
  scale_fill_manual(
    values=c(
      Positive=cor_pos,
      Negative=cor_neg
    )
  )+
  
  # legenda das magnitudes melhor
  scale_size_continuous(
    name = "Effect size",
    range = c(4,10),                 # mesmos tamanhos do gráfico
    breaks = c(0.1,0.2,0.3,0.4,0.5),
    guide = guide_legend(
      order = 2,
      title.position="top",
      label.position="right",
      keyheight=unit(0.9,"cm"),
      keywidth=unit(1,"cm")
    )
  )+
  
  guides(
    fill = guide_legend(
      title="Direction",
      order=1
    ),
    
    size = guide_legend(
      override.aes = list(
        shape=21,
        fill="grey70",
        colour="grey30",
        stroke=.6
      ),
      keyheight=unit(1.1,"cm"),
      keywidth=unit(1.1,"cm")
    )
  )+
  
  labs(
    title="Effect Size and Significance",
    x="Standardized Estimate",
    y=expression(-log[10](italic(p)[FDR]))
  )+
  
  theme_classic(base_size=14)+
  theme(
    plot.title=
      element_text(
        face="bold",
        size=16
      ),
    
    axis.title=
      element_text(
        face="bold",
        size=12
      ),
    
    legend.position="right",
    
    legend.title=
      element_text(
        face="bold",
        size=11
      ),
    
    legend.text=
      element_text(
        size=10
      )
  )

#=========================================================
# PANEL B 
#=========================================================

dados_fig <- dados_fig %>%
  mutate(
    # Suaviza a escala para os círculos não crescerem desgovernados
    Significancia_Tam = -log10(p_FDR)
  )

p2 <- ggplot(
  dados_fig,
  aes(
    x = Medida,                  # Invertido: Trait no X
    y = Variavel_Ambiental       # Invertido: Driver no Y
  )) +
  
  # Linhas de grade discretas ao fundo para guiar a leitura estruturada
  geom_hline(yintercept = unique(dados_fig$Variavel_Ambiental), color = "grey94", linewidth = 0.4) +
  geom_vline(xintercept = unique(dados_fig$Medida), color = "grey94", linewidth = 0.4) +
  
  # Os balões (tamanhos ligeiramente mais controlados)
  geom_point(
    aes(
      fill = Estimativa,
      size = Significancia_Tam
    ),
    shape = 21,
    color = "white",
    stroke = 0.8,
    alpha = 0.9
  ) +
  
  # Rótulos dos valores de P: Centralizados para evitar qualquer corte nas bordas
  geom_text(
    aes(label = p_lab),
    size = 2.8,               # Um pouco menor para caber melhor dentro das bolhas
    fontface = "bold",        
    color = "black",          
    hjust = 0.5,              # Centraliza horizontalmente no meio do círculo
    vjust = 0.5               # Centraliza verticalmente no meio do círculo
  ) +
  
  # Paleta de cores oficial
  scale_fill_gradient2(
    low = cor_neg,
    mid = "white",
    high = cor_pos,
    midpoint = 0,
    name = "Estimate"
  ) +
  
  # Escala de tamanho calibrada e sem legenda indesejada
  scale_size_continuous(
    range = c(2.5, 7.5), 
    guide = "none" 
  ) +
  
  labs(
    title = "Global Relationship Matrix",
    x = "Morphological Traits",
    y = "Environmental Drivers"
  ) +
  
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    panel.grid = element_blank(),
    
    # Margem direita expandida para o texto do P na última coluna não cortar
    plot.margin = margin(t = 20, r = 35, b = 20, l = 20, unit = "pt"),
    
    # Textos dos eixos retos e limpos, já que não há mais aperto
    axis.text.x = element_text(
      size = 10,
      face = "bold",
      colour = "grey10"
    ),
    
    axis.text.y = element_text(
      size = 9,                 # Um tiquinho menor para listar todas as variáveis perfeitamente
      face = "bold",
      colour = "grey10"
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 12
    ),
    
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10)
  )

print(p2)

#=========================================================
# PANEL C — VOLTA AO CIRCULAR ORIGINAL (ARRUMADO)
#=========================================================

edges <- dados_fig %>%
  filter(p_FDR < 0.05) %>%
  rename(
    from=Variavel_Ambiental,
    to=Medida
  )

nodes <- data.frame(
  name=unique(
    c(edges$from,edges$to)
  )) %>%
  mutate(
    Type=
      ifelse(
        name %in% edges$to,
        "Morphology",
        "Environment"
      )
  )

g <- graph_from_data_frame(
  edges,
  vertices=nodes,
  directed=FALSE
)

p3 <- ggraph(
  g,
  layout="linear",
  circular=TRUE
)+
  
  geom_edge_arc(
    aes(
      width=Weight,
      color=Direction
    ),
    alpha=.55
  )+
  
  geom_node_point(
    aes(fill=Type),
    shape=21,
    size=8,
    stroke=1.3,
    color="white"
  ) +
  
  # CORREÇÃO AQUI: Mudado para usar o argumento repel correto
  geom_node_text(
    aes(label=name),
    size=4,
    fontface="bold",
    repel=TRUE                  # Ativa o ggrepel por dentro do ggraph
  ) +
  
  scale_fill_manual(
    values=c(
      Morphology=cor_node_morph,
      Environment=cor_node_env
    )
  )+
  
  scale_edge_color_manual(
    values=c(
      Positive=cor_pos,
      Negative=cor_neg
    )
  )+
  
  scale_edge_width(
    range=c(.6,2.6)
  )+
  
  theme_void()+
  labs(
    title="Significant Interactions Network",
    fill="Node Class",
    edge_color="Direction",
    edge_width="Effect Size"
  )+
  
  theme(
    plot.title=
      element_text(
        face="bold",
        size=16,
        hjust=.5
      ),
    
    legend.position="bottom",
    legend.box="horizontal",
    
    legend.title=
      element_text(
        face="bold",
        size=11
      ),
    
    legend.text=
      element_text(
        size=10
      )
  )+
  
  guides(
    fill=guide_legend(nrow=1),
    edge_colour=guide_legend(nrow=1),
    edge_width=guide_legend(nrow=1)
  )

print(p3)

#=========================================================
# COMPOSIÇÃO FINAL MAIS BALANCEADA
#=========================================================

fig_final <-
  (p1 | p2) /
  p3 +
  plot_layout(
    heights=c(1,1.15)   # rede não exagerada
  )

print(fig_final)

ggsave(
  "Integrated_Figure_Final.png",
  fig_final,
  width=16,
  height=12,
  dpi=600,
  bg="white"
)

ggsave(
  "Panel_C_network.png",
  plot = p3,
  width = 16,
  height = 10,
  dpi = 600,
)

# ==============================================================================
# ALOMETRIA E INTEGRAÇÃO MORFOLÓGICA
# Abordagem revisada: SMA + hipóteses dimensionais específicas
# ==============================================================================

# DETXAR A TABELA NO MODELO DA ANÁLISE
library(dplyr)
library(tidyr)

df_wide <- df_processado %>%
  pivot_wider(
    names_from = Variavel,
    values_from = Valor
  )

str(df_wide)

names(df_wide)











# ==============================================================================
# ALOMETRIA E INTEGRAÇÃO MORFOLÓGICA
# Abordagem revisada: SMA + hipóteses dimensionais específicas
# ==============================================================================

cat("\n--- Iniciando análise revisada de alometria e integração ---\n")


# ==============================================================================
# 0. PACOTES
# ==============================================================================

pacotes <- c(
  "dplyr",
  "tidyr",
  "readr",
  "ggplot2",
  "GGally",
  "smatr",
  "writexl"
)

faltantes <- pacotes[
  !pacotes %in% rownames(installed.packages())
]

if(length(faltantes) > 0){
  install.packages(faltantes)
}

invisible(
  lapply(
    pacotes,
    library,
    character.only = TRUE
  )
)


# ==============================================================================
# 1. PALETA DE CORES
# ==============================================================================

COL_NAVY <- "#16425B"   # Azul Marinho Colonial
COL_GOLD <- "#D9A05B"   # Ouro Velho
COL_RED  <- "#A8422B"   # Vermelho Sangue de Boi
COL_PALE <- "#CBD4C2"   # Cinza Esverdeado Pálido


# ==============================================================================
# 2. DEFINIR A BASE DE DADOS
# ==============================================================================
#
# Se df_wide já existir, ele será utilizado diretamente.
#
# Caso contrário, o script tenta gerar df_wide a partir de df_processado.
#
# Estrutura esperada:
#
# id | Herbario | bioma | Lat | Log |
# AAN | BAN | LBL | TOL | PEL | PED | RHD | LFA
#
# ==============================================================================

if(exists("df_wide")){
  
  df_analise <- df_wide
  
} else if(exists("df_processado")){
  
  df_analise <- df_processado %>%
    pivot_wider(
      names_from = Variavel,
      values_from = Valor
    )
  
} else {
  
  stop(
    "ERRO: nem 'df_wide' nem 'df_processado' foram encontrados no ambiente."
  )
  
}

cat("\nBase de dados definida com sucesso.\n")


# ==============================================================================
# 3. VERIFICAR VARIÁVEIS NECESSÁRIAS
# ==============================================================================

variaveis_necessarias <- c(
  "AAN",
  "BAN",
  "LBL",
  "TOL",
  "PEL",
  "PED",
  "RHD",
  "LFA"
)

variaveis_faltando <- setdiff(
  variaveis_necessarias,
  names(df_analise)
)

if(length(variaveis_faltando) > 0){
  
  stop(
    paste(
      "ERRO: faltam as seguintes variáveis:",
      paste(
        variaveis_faltando,
        collapse = ", "
      )
    )
  )
  
}

cat("Todas as variáveis morfométricas foram encontradas.\n")


# ==============================================================================
# 4. VERIFICAR SE AS VARIÁVEIS SÃO NUMÉRICAS
# ==============================================================================

nao_numericas <- variaveis_necessarias[
  !sapply(
    df_analise[variaveis_necessarias],
    is.numeric
  )
]

if(length(nao_numericas) > 0){
  
  stop(
    paste(
      "ERRO: estas variáveis não são numéricas:",
      paste(
        nao_numericas,
        collapse = ", "
      )
    )
  )
  
}

cat("Todas as variáveis morfométricas são numéricas.\n")


# ==============================================================================
# 5. DEFINIÇÃO DAS HIPÓTESES ALOMÉTRICAS
# ==============================================================================
#
# OBJETIVO
# ------------------------------------------------------------------------------
#
# A análise anterior utilizava regressão OLS e comparava automaticamente todos
# os slopes com b = 1.
#
# Isso foi alterado porque a expectativa de isometria depende da dimensão
# física das variáveis.
#
#
# RELAÇÃO ALOMÉTRICA:
#
#       Y = a * X^b
#
# Em escala logarítmica:
#
#       log(Y) = log(a) + b * log(X)
#
#
# b  = slope alométrico observado
# b0 = slope esperado sob similaridade geométrica
#
#
# ------------------------------------------------------------------------------
# COMPRIMENTO × COMPRIMENTO
# ------------------------------------------------------------------------------
#
# As duas variáveis possuem dimensão L.
#
# Portanto:
#
#       b0 = 1
#
#
# ------------------------------------------------------------------------------
# ÁREA × COMPRIMENTO
# ------------------------------------------------------------------------------
#
# Área possui dimensão L².
# Comprimento possui dimensão L.
#
# Sob similaridade geométrica:
#
#       Área ~ Comprimento²
#
# Portanto:
#
#       b0 = 2
#
#
# ------------------------------------------------------------------------------
# POR QUE AAN E BAN NÃO ENTRAM NA CLASSIFICAÇÃO ALOMÉTRICA?
# ------------------------------------------------------------------------------
#
# AAN = Apex angle
# BAN = Base angle
#
# São ângulos e não possuem uma expectativa dimensional simples de isometria.
#
# Portanto, continuam na análise de integração/covariação, mas NÃO recebem
# classificação de:
#
# - Isometric
# - Positive allometry
# - Negative allometry
#
#
# ------------------------------------------------------------------------------
# POR QUE TOL NÃO ENTRA NA ALOMETRIA FORMAL?
# ------------------------------------------------------------------------------
#
# TOL = Total length
#
# É uma variável composta contendo principalmente:
#
#       LBL + PEL
#
# Relações como TOL × LBL ou TOL × PEL possuem, portanto, dependência
# matemática.
#
# TOL permanece na matriz de integração morfológica, mas não é utilizado como
# evidência independente de alometria.
#
#
# ------------------------------------------------------------------------------
# POR QUE NÃO "TODOS CONTRA TODOS"?
# ------------------------------------------------------------------------------
#
# Cada par aparece apenas UMA vez.
#
# SMA é utilizada para estudar a relação de escala entre dois caracteres,
# evitando interpretar X -> Y e Y -> X como duas relações biológicas
# independentes.
#
# ==============================================================================


# ==============================================================================
# 6. HIPÓTESES QUE SERÃO TESTADAS
# ==============================================================================

hipoteses_alometricas <- tibble::tribble(
  
  ~Variavel_X, ~Variavel_Y, ~b_esperado, ~Hipotese,
  
  
  # ---------------------------------------------------------------------------
  # LENGTH × LENGTH
  # ---------------------------------------------------------------------------
  
  "LBL", "PEL", 1, "Length vs length",
  "LBL", "PED", 1, "Length vs length",
  "LBL", "RHD", 1, "Length vs length",
  
  "PEL", "PED", 1, "Length vs length",
  "PEL", "RHD", 1, "Length vs length",
  
  "PED", "RHD", 1, "Length vs length",
  
  
  # ---------------------------------------------------------------------------
  # AREA × LENGTH
  # ---------------------------------------------------------------------------
  
  "LBL", "LFA", 2, "Area vs length",
  "PEL", "LFA", 2, "Area vs length",
  "PED", "LFA", 2, "Area vs length",
  "RHD", "LFA", 2, "Area vs length"
  
)


cat("\n--- Hypotheses to be tested ---\n")

print(
  hipoteses_alometricas
)


# ==============================================================================
# 7. CONFIGURAÇÕES
# ==============================================================================

ALPHA <- 0.05


# Ajuste para múltiplos testes:
#
# "BH" = Benjamini-Hochberg
#
# Se quiser usar Bonferroni:
#
# METODO_AJUSTE <- "bonferroni"

METODO_AJUSTE <- "BH"


# ==============================================================================
# 8. FUNÇÃO PARA EXTRAIR RESULTADOS DO SMA
# ==============================================================================

extrair_resultados_sma <- function(modelo){
  
  resumo_modelo <- modelo$groupsummary
  
  
  # ---------------------------------------------------------------------------
  # Slope SMA
  # ---------------------------------------------------------------------------
  
  b_estimado <- suppressWarnings(
    as.numeric(
      resumo_modelo$Slope[1]
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # Intervalo de confiança de 95%
  # ---------------------------------------------------------------------------
  
  ic_inferior <- suppressWarnings(
    as.numeric(
      resumo_modelo$Slope_lowCI[1]
    )
  )
  
  ic_superior <- suppressWarnings(
    as.numeric(
      resumo_modelo$Slope_highCI[1]
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # p-valor do teste:
  #
  # H0: b = b0
  # ---------------------------------------------------------------------------
  
  p_teste_b <- NA_real_
  
  
  if("Slope_test_p" %in% names(resumo_modelo)){
    
    p_teste_b <- suppressWarnings(
      as.numeric(
        resumo_modelo$Slope_test_p[1]
      )
    )
    
  }
  
  
  # ---------------------------------------------------------------------------
  # Fallback
  #
  # Em algumas versões do smatr o p-valor fica armazenado diretamente dentro
  # de slopetest.
  # ---------------------------------------------------------------------------
  
  if(is.na(p_teste_b)){
    
    if(
      !is.null(modelo$slopetest) &&
      length(modelo$slopetest) >= 1 &&
      !is.null(modelo$slopetest[[1]]) &&
      !is.null(modelo$slopetest[[1]]$p)
    ){
      
      p_teste_b <- suppressWarnings(
        as.numeric(
          modelo$slopetest[[1]]$p
        )
      )
      
    }
    
  }
  
  
  # ---------------------------------------------------------------------------
  # R²
  # ---------------------------------------------------------------------------
  
  r2 <- suppressWarnings(
    as.numeric(
      resumo_modelo$r2[1]
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # p-valor da relação entre X e Y
  # ---------------------------------------------------------------------------
  
  p_relacao <- suppressWarnings(
    as.numeric(
      resumo_modelo$pval[1]
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # Número de observações
  # ---------------------------------------------------------------------------
  
  n_amostras <- suppressWarnings(
    as.numeric(
      resumo_modelo$n[1]
    )
  )
  
  
  # ---------------------------------------------------------------------------
  # Retornar resultados
  # ---------------------------------------------------------------------------
  
  data.frame(
    
    N =
      n_amostras,
    
    R2 =
      r2,
    
    P_Relacao =
      p_relacao,
    
    Coeficiente_b =
      b_estimado,
    
    IC95_inferior =
      ic_inferior,
    
    IC95_superior =
      ic_superior,
    
    P_Teste_b =
      p_teste_b,
    
    stringsAsFactors =
      FALSE
    
  )
  
}


# ==============================================================================
# 9. RODAR OS MODELOS SMA
# ==============================================================================

resultados_sma <- vector(
  "list",
  nrow(hipoteses_alometricas)
)

modelos_sma <- vector(
  "list",
  nrow(hipoteses_alometricas)
)


names(modelos_sma) <- paste(
  hipoteses_alometricas$Variavel_Y,
  hipoteses_alometricas$Variavel_X,
  sep = "_vs_"
)


for(i in seq_len(nrow(hipoteses_alometricas))){
  
  var_x <-
    hipoteses_alometricas$Variavel_X[i]
  
  var_y <-
    hipoteses_alometricas$Variavel_Y[i]
  
  b0 <-
    hipoteses_alometricas$b_esperado[i]
  
  
  cat(
    "\n---------------------------------------------\n"
  )
  
  cat(
    "Analyzing:",
    var_y,
    "~",
    var_x,
    "\n"
  )
  
  cat(
    "Expected slope (b0):",
    b0,
    "\n"
  )
  
  
  # ---------------------------------------------------------------------------
  # Selecionar dados válidos
  # ---------------------------------------------------------------------------
  
  df_par <- df_analise %>%
    
    transmute(
      
      Eixo_X =
        .data[[var_x]],
      
      Eixo_Y =
        .data[[var_y]]
      
    ) %>%
    
    filter(
      
      !is.na(Eixo_X),
      
      !is.na(Eixo_Y),
      
      Eixo_X > 0,
      
      Eixo_Y > 0
      
    )
  
  
  cat(
    "Valid observations:",
    nrow(df_par),
    "\n"
  )
  
  
  # ---------------------------------------------------------------------------
  # Rodar SMA
  # ---------------------------------------------------------------------------
  
  if(nrow(df_par) > 5){
    
    
    modelo <- smatr::sma(
      
      Eixo_Y ~ Eixo_X,
      
      data =
        df_par,
      
      log =
        "xy",
      
      method =
        "SMA",
      
      slope.test =
        b0,
      
      alpha =
        ALPHA
      
    )
    
    
    # -------------------------------------------------------------------------
    # Guardar modelo
    # -------------------------------------------------------------------------
    
    modelos_sma[[i]] <- list(
      
      modelo =
        modelo,
      
      dados =
        df_par,
      
      x =
        var_x,
      
      y =
        var_y,
      
      b0 =
        b0
      
    )
    
    
    # -------------------------------------------------------------------------
    # Extrair resultados
    # -------------------------------------------------------------------------
    
    res_tmp <- extrair_resultados_sma(
      modelo
    )
    
    
    resultados_sma[[i]] <- data.frame(
      
      Variavel_X =
        var_x,
      
      Variavel_Y =
        var_y,
      
      b_esperado =
        b0,
      
      Hipotese =
        hipoteses_alometricas$Hipotese[i],
      
      res_tmp,
      
      stringsAsFactors =
        FALSE
      
    )
    
    
  } else {
    
    
    warning(
      paste(
        "Poucos dados para:",
        var_y,
        "~",
        var_x
      )
    )
    
    
  }
  
}


# ==============================================================================
# 10. JUNTAR RESULTADOS
# ==============================================================================

tabela_sma <- bind_rows(
  resultados_sma
)


# ==============================================================================
# 11. AJUSTAR P-VALORES
# ==============================================================================

tabela_sma <- tabela_sma %>%
  
  mutate(
    
    P_Relacao_Ajustado =
      p.adjust(
        P_Relacao,
        method = METODO_AJUSTE
      ),
    
    P_Teste_b_Ajustado =
      p.adjust(
        P_Teste_b,
        method = METODO_AJUSTE
      )
    
  )


# ==============================================================================
# 12. CLASSIFICAÇÃO ALOMÉTRICA
# ==============================================================================

tabela_sma <- tabela_sma %>%
  
  mutate(
    
    Classificacao = case_when(
      
      is.na(P_Relacao_Ajustado) ~
        "Insufficient data",
      
      P_Relacao_Ajustado >= ALPHA ~
        "Not significant",
      
      is.na(P_Teste_b_Ajustado) ~
        "Slope test unavailable",
      
      P_Teste_b_Ajustado >= ALPHA ~
        "Isometric",
      
      P_Teste_b_Ajustado < ALPHA &
        Coeficiente_b > b_esperado ~
        "Positive allometry",
      
      P_Teste_b_Ajustado < ALPHA &
        Coeficiente_b < b_esperado ~
        "Negative allometry",
      
      TRUE ~
        "Unclassified"
      
    )
    
  )


# ==============================================================================
# 13. FORMATAR INTERVALO DE CONFIANÇA
# ==============================================================================

tabela_sma <- tabela_sma %>%
  
  mutate(
    
    IC95 = paste0(
      
      round(
        IC95_inferior,
        3
      ),
      
      " – ",
      
      round(
        IC95_superior,
        3
      )
      
    ),
    
    Relacao = paste0(
      Variavel_Y,
      " ~ ",
      Variavel_X
    )
    
  )


# ==============================================================================
# 14. TABLE S13
# ==============================================================================

Table_S13 <- tabela_sma %>%
  
  transmute(
    
    `Variable X` =
      Variavel_X,
    
    `Variable Y` =
      Variavel_Y,
    
    N =
      N,
    
    `R²` =
      round(
        R2,
        3
      ),
    
    `SMA slope (b)` =
      round(
        Coeficiente_b,
        3
      ),
    
    `95% CI` =
      IC95,
    
    `Expected slope (b0)` =
      b_esperado,
    
    `p-value relationship` =
      signif(
        P_Relacao,
        4
      ),
    
    `Adjusted p-value relationship` =
      signif(
        P_Relacao_Ajustado,
        4
      ),
    
    `p-value slope test` =
      signif(
        P_Teste_b,
        4
      ),
    
    `Adjusted p-value slope test` =
      signif(
        P_Teste_b_Ajustado,
        4
      ),
    
    Classification =
      Classificacao
    
  )


cat(
  "\n\n--- TABLE S13 ---\n"
)


print(
  dplyr::as_tibble(Table_S13),
  n = Inf
)


# ==============================================================================
# 15. SALVAR TABLE S13
# ==============================================================================

readr::write_csv2(
  
  Table_S13,
  
  "Table_S13_SMA_allometry.csv"
  
)


writexl::write_xlsx(
  
  list(
    
    "Table S13" =
      Table_S13,
    
    "Allometric hypotheses" =
      hipoteses_alometricas
    
  ),
  
  "Table_S13_SMA_allometry.xlsx"
  
)


cat(
  "\nTable S13 salva em CSV e XLSX.\n"
)


# ==============================================================================
# 16. SMA DIAGNOSTICS - PRANCHA ÚNICA EM PNG
# ==============================================================================
#
# MANTÉM A ORGANIZAÇÃO ORIGINAL:
#
# - 4 gráficos por linha;
# - resíduos e Q-Q;
# - sem título geral;
# - títulos individuais maiores;
# - letras dos eixos maiores.
#
# ==============================================================================


modelos_validos <- modelos_sma[
  !sapply(
    modelos_sma,
    is.null
  )
]


n_modelos <- length(
  modelos_validos
)


n_plots <- n_modelos * 2


# ------------------------------------------------------------------------------
# 4 gráficos por linha
# ------------------------------------------------------------------------------

ncol_diag <- 4

nrow_diag <- ceiling(
  n_plots / ncol_diag
)


# ------------------------------------------------------------------------------
# Abrir PNG
# ------------------------------------------------------------------------------

png(
  
  filename =
    "SMA_diagnostics_panel.png",
  
  width =
    4200,
  
  height =
    max(
      3200,
      nrow_diag * 900
    ),
  
  res =
    300
  
)


# ------------------------------------------------------------------------------
# Configuração visual
# ------------------------------------------------------------------------------

par(
  
  mfrow =
    c(
      nrow_diag,
      ncol_diag
    ),
  
  mar =
    c(
      4.2,
      4.2,
      4.5,
      1.2
    ),
  
  cex.main =
    1.35,
  
  cex.axis =
    1.05,
  
  cex.lab =
    1.10
  
)


# ------------------------------------------------------------------------------
# Gerar gráficos
# ------------------------------------------------------------------------------

for(obj in modelos_validos){
  
  
  nome_modelo <- paste0(
    obj$y,
    " ~ ",
    obj$x
  )
  
  
  # ---------------------------------------------------------------------------
  # Resíduos
  # ---------------------------------------------------------------------------
  
  plot(
    
    obj$modelo,
    
    which =
      "residual",
    
    main =
      paste0(
        nome_modelo,
        "\nResiduals"
      )
    
  )
  
  
  # ---------------------------------------------------------------------------
  # Q-Q plot
  # ---------------------------------------------------------------------------
  
  plot(
    
    obj$modelo,
    
    which =
      "qq",
    
    main =
      paste0(
        nome_modelo,
        "\nQ-Q plot"
      )
    
  )
  
}


# ------------------------------------------------------------------------------
# Fechar PNG
# ------------------------------------------------------------------------------

dev.off()


par(
  mfrow =
    c(
      1,
      1
    )
)


cat(
  "\nDiagnósticos salvos em 'SMA_diagnostics_panel.png'.\n"
)


# ==============================================================================
# 17. FIGURA PRINCIPAL
# SMA FOREST PLOT
# ==============================================================================
#
# FIGURA DE ALOMETRIA
#
# - códigos das variáveis na lateral;
# - IC95% horizontal;
# - pequenas barras verticais nas extremidades do IC;
# - ponto = slope SMA;
# - valor de b escrito à direita;
# - linha tracejada = slope esperado;
# - classificação indicada pelas cores.
#
# ==============================================================================


# ------------------------------------------------------------------------------
# Preparar dados
# ------------------------------------------------------------------------------

tabela_fig <- tabela_sma %>%
  
  mutate(
    
    # -------------------------------------------------------------------------
    # Mostrar apenas os códigos
    # -------------------------------------------------------------------------
    
    Relacao =
      paste0(
        Variavel_Y,
        " ~ ",
        Variavel_X
      ),
    
    
    # -------------------------------------------------------------------------
    # Valor do slope
    # -------------------------------------------------------------------------
    
    Label_b =
      paste0(
        "b = ",
        sprintf(
          "%.2f",
          Coeficiente_b
        )
      ),
    
    
    # -------------------------------------------------------------------------
    # Ordem dos painéis
    # -------------------------------------------------------------------------
    
    Hipotese =
      factor(
        
        Hipotese,
        
        levels =
          c(
            "Length vs length",
            "Area vs length"
          )
        
      )
    
  ) %>%
  
  arrange(
    Hipotese,
    Coeficiente_b
  ) %>%
  
  mutate(
    
    Relacao =
      factor(
        Relacao,
        levels = unique(Relacao)
      )
    
  )


# ------------------------------------------------------------------------------
# Valores esperados
# ------------------------------------------------------------------------------

ref_df <- data.frame(
  
  Hipotese =
    factor(
      
      c(
        "Length vs length",
        "Area vs length"
      ),
      
      levels =
        c(
          "Length vs length",
          "Area vs length"
        )
      
    ),
  
  b0 =
    c(
      1,
      2
    )
  
)


# ------------------------------------------------------------------------------
# Paleta das classificações
# ------------------------------------------------------------------------------

paleta_class <- c(
  
  "Isometric" =
    COL_GOLD,
  
  "Positive allometry" =
    COL_RED,
  
  "Negative allometry" =
    COL_NAVY,
  
  "Not significant" =
    COL_PALE,
  
  "Slope test unavailable" =
    COL_PALE,
  
  "Insufficient data" =
    COL_PALE,
  
  "Unclassified" =
    COL_PALE
  
)


# ------------------------------------------------------------------------------
# Posição dos valores de b
# ------------------------------------------------------------------------------

amplitude_x <- diff(
  
  range(
    
    c(
      tabela_fig$IC95_inferior,
      tabela_fig$IC95_superior
    ),
    
    na.rm =
      TRUE
    
  )
  
)


deslocamento_label <- max(
  0.06,
  amplitude_x * 0.025
)


tabela_fig <- tabela_fig %>%
  
  mutate(
    
    Label_x =
      IC95_superior +
      deslocamento_label
    
  )


# ==============================================================================
# GERAR FOREST PLOT
# ==============================================================================

fig_sma_forest <- ggplot(
  
  tabela_fig,
  
  aes(
    x = Coeficiente_b,
    y = Relacao
  )
  
) +
  
  
  # ---------------------------------------------------------------------------
# Linha tracejada = slope esperado
# ---------------------------------------------------------------------------

geom_vline(
  
  data =
    ref_df,
  
  aes(
    xintercept =
      b0
  ),
  
  inherit.aes =
    FALSE,
  
  linetype =
    "dashed",
  
  linewidth =
    0.85,
  
  color =
    COL_NAVY
  
) +
  
  
  # ---------------------------------------------------------------------------
# Intervalo de confiança com barrinhas nas extremidades
# ---------------------------------------------------------------------------

geom_errorbar(
  
  aes(
    
    xmin =
      IC95_inferior,
    
    xmax =
      IC95_superior,
    
    color =
      Classificacao
    
  ),
  
  orientation =
    "y",
  
  width =
    0.22,
  
  linewidth =
    1.0
  
) +
  
  
  # ---------------------------------------------------------------------------
# Ponto central = slope SMA
# ---------------------------------------------------------------------------

geom_point(
  
  aes(
    fill =
      Classificacao
  ),
  
  shape =
    21,
  
  size =
    3.8,
  
  stroke =
    0.75,
  
  color =
    "black"
  
) +
  
  
  # ---------------------------------------------------------------------------
# Valor de b
# ---------------------------------------------------------------------------

geom_text(
  
  aes(
    
    x =
      Label_x,
    
    label =
      Label_b
    
  ),
  
  size =
    3.35,
  
  hjust =
    0,
  
  color =
    "black"
  
) +
  
  
  # ---------------------------------------------------------------------------
# Painéis
# ---------------------------------------------------------------------------

facet_grid(
  
  Hipotese ~ .,
  
  scales =
    "free_y",
  
  space =
    "free_y"
  
) +
  
  
  # ---------------------------------------------------------------------------
# Cores
# ---------------------------------------------------------------------------

scale_fill_manual(
  values =
    paleta_class
) +
  
  scale_color_manual(
    values =
      paleta_class
  ) +
  
  
  # ---------------------------------------------------------------------------
# Espaço extra para valores de b
# ---------------------------------------------------------------------------

scale_x_continuous(
  
  expand =
    expansion(
      mult =
        c(
          0.05,
          0.20
        )
    )
  
) +
  
  
  # ---------------------------------------------------------------------------
# Labels
# ---------------------------------------------------------------------------

labs(
  
  x =
    "SMA slope (b)",
  
  y =
    NULL,
  
  fill =
    NULL,
  
  color =
    NULL
  
) +
  
  
  # ---------------------------------------------------------------------------
# Tema
# ---------------------------------------------------------------------------

theme_minimal(
  base_size =
    12
) +
  
  theme(
    
    panel.grid.minor =
      element_blank(),
    
    panel.grid.major.y =
      element_blank(),
    
    panel.grid.major.x =
      element_line(
        color = "grey88",
        linewidth = 0.4
      ),
    
    strip.text =
      element_text(
        face = "bold",
        size = 11,
        color = COL_NAVY
      ),
    
    strip.background =
      element_rect(
        fill = COL_PALE,
        color = COL_NAVY,
        linewidth = 0.6
      ),
    
    legend.position =
      "bottom",
    
    legend.box =
      "horizontal",
    
    axis.text.y =
      element_text(
        size = 11,
        face = "bold",
        color = "black"
      ),
    
    axis.text.x =
      element_text(
        size = 10,
        color = "black"
      ),
    
    axis.title.x =
      element_text(
        size = 11,
        face = "bold",
        color = "black"
      ),
    
    panel.border =
      element_rect(
        color = "grey78",
        fill = NA,
        linewidth = 0.4
      ),
    
    panel.spacing =
      unit(
        0.7,
        "lines"
      )
    
  )


# ==============================================================================
# EXIBIR FOREST PLOT
# ==============================================================================

print(
  fig_sma_forest
)


# ==============================================================================
# SALVAR FOREST PLOT
# ==============================================================================

ggsave(
  
  filename =
    "Figure_SMA_forestplot.png",
  
  plot =
    fig_sma_forest,
  
  width =
    8.5,
  
  height =
    7.5,
  
  dpi =
    300,
  
  bg =
    "white"
  
)


cat(
  "\nFigura salva como 'Figure_SMA_forestplot.png'.\n"
)


# ==============================================================================
# 18. MATRIZ DE INTEGRAÇÃO MORFOLÓGICA
# ==============================================================================
#
# Esta figura mostra:
#
# - distribuição dos caracteres;
# - dispersão;
# - correlação de Spearman;
# - integração/covariação morfológica.
#
# AAN, BAN e TOL permanecem aqui.
#
# ==============================================================================


df_integracao <- df_analise %>%
  
  select(
    
    AAN,
    BAN,
    LBL,
    TOL,
    PEL,
    PED,
    RHD,
    LFA
    
  )


# ==============================================================================
# GERAR MATRIZ
# ==============================================================================

matriz_integracao <- GGally::ggpairs(
  
  df_integracao,
  
  
  # ---------------------------------------------------------------------------
  # Parte inferior
  # ---------------------------------------------------------------------------
  
  lower = list(
    
    continuous =
      GGally::wrap(
        
        "points",
        
        alpha =
          0.50,
        
        size =
          1.0,
        
        color =
          COL_NAVY
        
      )
    
  ),
  
  
  # ---------------------------------------------------------------------------
  # Diagonal
  # ---------------------------------------------------------------------------
  
  diag = list(
    
    continuous =
      GGally::wrap(
        
        "densityDiag",
        
        fill =
          COL_PALE,
        
        alpha =
          0.85,
        
        color =
          COL_NAVY
        
      )
    
  ),
  
  
  # ---------------------------------------------------------------------------
  # Parte superior
  # ---------------------------------------------------------------------------
  
  upper = list(
    
    continuous =
      GGally::wrap(
        
        "cor",
        
        method =
          "spearman",
        
        size =
          4,
        
        color =
          COL_RED
        
      )
    
  )
  
) +
  
  
  theme_minimal(
    base_size =
      11
  ) +
  
  
  theme(
    
    strip.text =
      element_text(
        face = "bold",
        size = 10,
        color = COL_NAVY
      ),
    
    strip.background =
      element_rect(
        fill = COL_PALE,
        color = COL_NAVY
      ),
    
    panel.grid.minor =
      element_blank(),
    
    panel.grid.major =
      element_line(
        color = "grey92"
      ),
    
    panel.border =
      element_rect(
        color = COL_NAVY,
        fill = NA,
        linewidth = 0.5
      )
    
  )


# ==============================================================================
# EXIBIR MATRIZ
# ==============================================================================

print(
  matriz_integracao
)


# ==============================================================================
# SALVAR MATRIZ
# ==============================================================================

ggsave(
  
  filename =
    "matriz_integracao_morfologica.png",
  
  plot =
    matriz_integracao,
  
  width =
    12,
  
  height =
    12,
  
  dpi =
    300,
  
  bg =
    "white"
  
)


cat(
  "\nFigura suplementar salva como 'matriz_integracao_morfologica.png'.\n"
)


# ==============================================================================
# 19. RESUMO FINAL
# ==============================================================================

cat(
  "\n====================================================\n"
)

cat(
  "ANÁLISE CONCLUÍDA\n"
)

cat(
  "====================================================\n"
)


cat(
  "\nArquivos gerados:\n"
)

cat(
  "1. Table_S13_SMA_allometry.csv\n"
)

cat(
  "2. Table_S13_SMA_allometry.xlsx\n"
)

cat(
  "3. SMA_diagnostics_panel.png\n"
)

cat(
  "4. Figure_SMA_forestplot.png\n"
)

cat(
  "5. matriz_integracao_morfologica.png\n"
)


cat(
  "\nObservações:\n"
)

cat(
  "- AAN e BAN: mantidos apenas na integração/covariação.\n"
)

cat(
  "- TOL: mantido na integração, mas retirado da alometria formal.\n"
)

cat(
  "- Table S13 utiliza a classificação 'Isometric'.\n"
)

cat(
  "- SMA_diagnostics_panel mantém 4 gráficos por linha e letras maiores.\n"
)

cat(
  "- Figure_SMA_forestplot usa apenas os códigos das variáveis.\n"
)

cat(
  "- Os intervalos de confiança possuem barras nas duas extremidades.\n"
)

cat(
  "- A linha tracejada mostra o slope esperado sob isometria.\n"
)

cat(
  "- A matriz de integração pode ser utilizada como figura suplementar.\n"
)


cat(
  "\n--- Fim da análise ---\n"
)

