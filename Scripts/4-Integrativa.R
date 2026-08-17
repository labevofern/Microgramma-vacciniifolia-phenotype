# ==============================================================
# Análises de correlação ambiental
# ==============================================================
# Autor: Juliana Aljahara
# Programa: Doutoranda em Biologia Vegetal
# Instituição: Universidade Federal de Pernambuco (UFPE)
# Orientadora: Dra. Thaís E. Almeida
# ==========================================================
# 1. Pacotes
# ==========================================================

library(vegan)       # Análises de ecologia de comunidades (diversidade, ordenação e estatística multivariada)
library(dplyr)       # Manipulação de dados (filter, group_by, summarise)
library(ggplot2)     # Construção de gráficos
library(readxl)      # Importação de planilhas Excel (.xlsx)
library(missForest)  # Imputação de dados faltantes utilizando Random Forest
library(shapes)      # Análise estatística de formas e métodos de Procrustes
library(ggrepel)     # Melhor posicionamento automático de textos em gráficos (evita sobreposição)
library(patchwork)   # Combinação de múltiplos gráficos em um único layout

# ==========================================================
# 2. DIRETÓRIO E DADOS
# ==========================================================

setwd("C:/Users/nikso/OneDrive/_JULIANA ALJAHARA/1_Análise de macromorfologica/Diretório_Integrativa/2_Tabelas_e_Imagens")

# Matrizes de distância
dist_linear  <- read.csv("dist_linear.csv", row.names = 1)
dist_outline <- read.csv("dist_outline.csv", row.names = 1)

# Variáveis ambientais
env <- read.csv("Especimes_e_Variaveis_Selecionadas_VIF.csv", row.names = 1)

# ==========================================================
# 3. SELEÇÃO DE COLUNAS
# ==========================================================

cols <- c(1, 2, 3, 4, 5, 7, 8, 10:ncol(env))
env <- env[, cols]

# Garantir que rownames = herbário
rownames(env) <- env$Herbario

# ==========================================================
# 4. SINCRONIZAÇÃO ENTRE DATASETS
# ==========================================================

samples <- Reduce(intersect, list(
  rownames(env),
  rownames(dist_linear),
  rownames(dist_outline)
))

env <- env[samples, ]
dist_linear  <- dist_linear[samples, samples]
dist_outline <- dist_outline[samples, samples]

# Converter para dist
dist_linear  <- as.dist(dist_linear)
dist_outline <- as.dist(dist_outline)

# ==========================================================
# 5. DEFINIR VAIAVEIS NUMERICAS
# ==========================================================

# Apenas variáveis numéricas (a partir da coluna 8)
env_num <- env[, 8:ncol(env)]

# ==========================================================
# 6. FUNÇÃO PARA RODAR dbRDA (FILTRO P < 0.1)
# ==========================================================

run_dbrda <- function(dist_mat, env_data){
  
  # Modelo completo
  full <- dbrda(dist_mat ~ ., data = as.data.frame(env_data))
  ani_terms_full <- anova.cca(full, by = "terms", permutations = 9999)
  
  # Captura variáveis com p < 0.1 (inclui ***, **, * e .)
  p_values <- ani_terms_full$`Pr(>F)`
  names(p_values) <- rownames(ani_terms_full)
  sig_vars <- names(p_values[p_values < 0.1 & !is.na(p_values)])
  
  if(length(sig_vars) > 0) {
    form_sig <- as.formula(paste("dist_mat ~", paste(sig_vars, collapse = " + ")))
    model_plot <- dbrda(form_sig, data = as.data.frame(env_data))
  } else {
    model_plot <- dbrda(dist_mat ~ 1, data = as.data.frame(env_data))
  }
  
  ani_model <- anova.cca(model_plot, permutations = 9999)
  ani_axes  <- anova.cca(model_plot, by = "axis", permutations = 9999)
  R2adj     <- RsquareAdj(model_plot)$adj.r.squared
  
  return(list(model = model_plot, 
              ani_terms_full = ani_terms_full, 
              ani_model = ani_model,
              ani_axes = ani_axes,
              R2adj = R2adj))
}

# ==========================================================
# 7. FUNÇÃO PARA GERAR GRÁFICO (CORRIGIDA: LEGENDA DE TEXTO)
# ==========================================================

plot_dbrda <- function(res, titulo, grupo_col){
  
  # 1. Scores dos Sites (Tratamento para Rank 1)
  site <- scores(res$model, display = "sites", choices = c(1,2))
  site <- as.data.frame(site)
  
  if(!"dbRDA2" %in% colnames(site)){
    if("MDS1" %in% colnames(site)) colnames(site)[colnames(site) == "MDS1"] <- "dbRDA2"
    else site$dbRDA2 <- 0
    label_y <- "MDS1 (Resíduos)"
  } else {
    perc_y <- round(res$ani_axes$Variance[2] / sum(res$ani_axes$Variance) * 100, 1)
    label_y <- paste0("dbRDA2 (", perc_y, "%)")
  }
  
  perc_x <- round(res$ani_axes$Variance[1] / sum(res$ani_axes$Variance) * 100, 1)
  label_x <- paste0("dbRDA1 (", perc_x, "%)")
  site$Grupo <- env[rownames(site), grupo_col] 
  
  # 2. Base do gráfico
  p <- ggplot() +
    geom_hline(yintercept = 0, linetype="dashed", color="grey70", linewidth=0.5) +
    geom_vline(xintercept = 0, linetype="dashed", color="grey70", linewidth=0.5) +
    stat_ellipse(data = site, aes(dbRDA1, dbRDA2, color=Grupo), linewidth=0.6, type="t") +
    geom_point(data = site, aes(dbRDA1, dbRDA2, color=Grupo), shape=16, size=2.5, stroke=0.8) +
    scale_color_manual(values = c("#6A0572", "#AB83A1", "#F67280", "#F8B195")) +
    labs(x = label_x, y = label_y, title = titulo,
         subtitle = paste0("R² adj = ", round(res$R2adj*100,2), "% | p-modelo = ", 
                           signif(res$ani_model$`Pr(>F)`[1],3))) +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  # 3. Adição das Setas e Legenda de Significado
  has_vars <- length(labels(terms(res$model))) > 0
  if(has_vars) {
    env_scr <- scores(res$model, display = "bp", choices = c(1,2))
    if(is.null(env_scr) || nrow(env_scr)==0) env_scr <- res$model$CCA$biplot[,1:2]
    env_scr <- as.data.frame(env_scr)
    
    if(!"dbRDA2" %in% colnames(env_scr)) {
      if("MDS1" %in% colnames(env_scr)) colnames(env_scr)[colnames(env_scr) == "MDS1"] <- "dbRDA2"
      else env_scr$dbRDA2 <- 0
    }
    colnames(env_scr)[1:2] <- c("dbRDA1","dbRDA2")
    env_scr$Variable <- rownames(env_scr)
    
    term_table <- as.data.frame(res$ani_terms_full)
    env_scr$p_val <- term_table[rownames(env_scr), "Pr(>F)"]
    
    # Criar rótulos com símbolos para o gráfico
    env_scr$Label <- with(env_scr, ifelse(p_val <= 0.001, paste0(Variable, " ***"),
                                          ifelse(p_val <= 0.01,  paste0(Variable, " **"),
                                                 ifelse(p_val <= 0.05,  paste0(Variable, " *"),
                                                        ifelse(p_val <= 0.1,   paste0(Variable, " ."), Variable)))))
    
    # Mapeamento da Legenda de Texto
    env_scr$Sig_Legenda <- with(env_scr, ifelse(p_val <= 0.001, "*** (p < 0.001)",
                                                ifelse(p_val <= 0.01,  "** (p < 0.01)",
                                                       ifelse(p_val <= 0.05,  "* (p < 0.05)",
                                                              ". (p < 0.1)"))))
    
    mult <- attributes(env_scr)$arrow.mul; if(is.null(mult)) mult <- 1.5
    env_scr$dbRDA1 <- env_scr$dbRDA1 * mult
    env_scr$dbRDA2 <- env_scr$dbRDA2 * mult
    
    p <- p + 
      # Setas: Agora usamos 'colour' mapeado para a legenda, mas fixamos preto no scale
      geom_segment(data = env_scr, 
                   aes(x=0, y=0, xend=dbRDA1, yend=dbRDA2, linetype=Sig_Legenda),
                   arrow=arrow(length=unit(0.2,"cm"), type="closed"),
                   color="black") +
      
      # Nomes das variáveis com os símbolos
      geom_text_repel(data = env_scr, aes(dbRDA1, dbRDA2, label=Label),
                      size=3.5, fontface="italic", segment.color = 'transparent') +
      
      # Força a legenda a mostrar o que você pediu, mas mantém as linhas sólidas no gráfico se preferir
      # Aqui defini 'solid' para todos conforme seu pedido de setas iguais
      scale_linetype_manual(name = "Significância dos Termos",
                            values = setNames(rep("solid", length(unique(env_scr$Sig_Legenda))), 
                                              unique(env_scr$Sig_Legenda))) +
      
      # Ajuste final para a legenda não mostrar uma linha atravessada, apenas o texto
      guides(linetype = guide_legend(override.aes = list(arrow = NULL)))
  }
  
  return(p)
}

# ==========================================================
# 8. RODAR E EXIBIR
# ==========================================================

# Substituir os NAs de cada coluna pela mediana da própria coluna como foi apena 1 valor de alguma colunas fiz isso
env_num_imputado <- env_num
for(i in 1:ncol(env_num_imputado)) {
  if(is.numeric(env_num_imputado[,i])) {
    env_num_imputado[is.na(env_num_imputado[,i]), i] <- median(env_num_imputado[,i], na.rm = TRUE)
  }
}

# Rodar a análise com os dados imputados
res_linear  <- run_dbrda(dist_linear, env_num_imputado)
res_outline <- run_dbrda(dist_outline, env_num_imputado)

grupo <- "Dominio.fitogeografico"

plot_linear  <- plot_dbrda(res_linear,  "dbRDA - Dados Lineares", "Dominio.fitogeografico")
plot_outline <- plot_dbrda(res_outline, "dbRDA - Dados de Forma (Outline)", "Dominio.fitogeografico")

print(plot_linear)
print(plot_outline)

# Para dados lineares 
res_linear$ani_terms_full 
# Para dados outline 
res_outline$ani_terms_full 

# ==========================================================
# 9. EXPORTAÇÃO (ALTA QUALIDADE)
# ==========================================================

ggsave("dbRDA_linear.png", plot_linear,
       width = 18, height = 15, units = "cm", dpi = 600)

ggsave("dbRDA_outline.png", plot_outline,
       width = 18, height = 15, units = "cm", dpi = 600)

# ==========================================================
# 9. COMPARAÇÃO ENTRE MODELOS (PROCRUSTES)
# ==========================================================

cat("\n=========== PROCRUSTES TEST ===========\n")

# 1. Executa o teste (Certifique-se de usar o mesmo nome em todo o bloco)
proc_model <- protest(res_linear$model, res_outline$model, permutations = 999)

# 2. Extração e Preparação dos dados para o ggplot
df_x <- as.data.frame(proc_model$X[, 1:2])
df_y <- as.data.frame(proc_model$Yrot[, 1:2])

# Alterando os nomes das colunas para inglês
colnames(df_x) <- c("Dim1", "Dim2")
colnames(df_y) <- c("Dim1", "Dim2")

df_x$id <- seq_len(nrow(df_x))
df_y$id <- seq_len(nrow(df_y))

seg_data <- merge(df_x, df_y, by="id", suffixes = c(".x", ".y"))

p_proc <- ggplot(data = seg_data) +
  geom_segment(aes(x = Dim1.x, y = Dim2.x, xend = Dim1.y, yend = Dim2.y), 
               color = "grey70", linewidth = 0.6, alpha = 0.6) +
  geom_point(aes(x = Dim1.x, y = Dim2.x, color = "Linear"), size = 2.5) +
  geom_point(aes(x = Dim1.y, y = Dim2.y, color = "Outline"), size = 2.5) +
  scale_color_manual(values = c("Linear" = "#E41A1C", "Outline" = "#377EB8")) +
  labs(title = "Procrustes Comparison: Linear vs. Outline",
       subtitle = paste("Correlation:", round(proc_model$t0, 3), "| p-value:", proc_model$signif),
       x = "Dimension 1", y = "Dimension 2", color = "Model") + # Eixos em inglês
  theme_minimal() +
  theme(legend.position = "bottom")

print(p_proc)

ggsave("PROCRUSTES.png", p_proc,
       width = 18, height = 15, units = "cm", dpi = 600)

# ==========================================================
# 10. PROCESSAMENTO E EXPORTAÇÃO DE TABELAS (ATUALIZADO)
# ==========================================================

# Função para extrair métricas por variável (Termos)
extrair_tabela_resumo <- function(res, nome_modelo) {
  # Extrai a tabela de anova por termos
  tab <- as.data.frame(res$ani_terms_full)
  
  # Limpa linhas que não são variáveis (Residual, Total, etc)
  tab <- tab[!rownames(tab) %in% c("Residual", "Total"), ]
  
  # Cria o dataframe de saída
  df_out <- data.frame(
    Modelo = nome_modelo,
    Variavel = rownames(tab),
    Df = tab$Df,
    F_stat = round(tab$F, 3),
    p_valor = tab$`Pr(>F)`
  )
  
  return(df_out)
}

# Agora você pode rodar o bloco B sem erros:
tabela_variaveis <- rbind(
  extrair_tabela_resumo(res_linear, "Linear"),
  extrair_tabela_resumo(res_outline, "Outline")
)

cat("\n--- DETALHAMENTO POR VARIÁVEL ---\n")
print(tabela_variaveis)

# ==========================================================
# TABELA GLOBAL DOS MODELOS
# ==========================================================

extrair_global <- function(res, nome){
  
  # R2 bruto e ajustado
  r2 <- RsquareAdj(res$model)
  
  # anova global
  an <- as.data.frame(res$ani_model)
  
  # graus de liberdade do modelo e resíduo
  gl_modelo  <- an$Df[1]
  gl_residuo <- an$Df[2]
  
  data.frame(
    Modelo = nome,
    R2_Bruto = round(r2$r.squared,4),
    R2_Ajustado = round(r2$adj.r.squared,4),
    F = round(an$F[1],4),
    df_Modelo = gl_modelo,
    df_Residuo = gl_residuo,
    p_valor = an$`Pr(>F)`[1]
  )
  
}

tabela_modelos_globais <- rbind(
  extrair_global(res_linear,"Linear measurements"),
  extrair_global(res_outline,"Elliptical Fourier Analysis (EFA)")
)

print(tabela_modelos_globais)
cat("\n--- RESUMO GERAL DOS MODELOS (GLOBAL) ---\n")
print(tabela_modelos_globais)

# B. Tabela de Variáveis Ambientais (O que você já tinha)
tabela_variaveis <- rbind(
  extrair_tabela_resumo(res_linear, "Linear"),
  extrair_tabela_resumo(res_outline, "Outline")
)

# C. Tabela Procrustes
proc_res <- data.frame(
  Metrica = c("Correlation_r", "Sum_of_Squares_m2", "p-value"),
  Valor = round(c(proc_model$t0, proc_model$ss, proc_model$signif), 4)
)

# EXPORTAÇÃO FINAL
write.csv(tabela_modelos_globais, "Resumo_Global_Modelos_dbRDA.csv", row.names = FALSE)
write.csv(tabela_variaveis, "Suplemento_Variaveis_dbRDA.csv", row.names = FALSE)
write.csv(proc_res, "Resultado_Procrustes.csv", row.names = FALSE)

cat("\n--- TODAS AS TABELAS FORAM EXPORTADAS COM SUCESSO ---\n")

# ==========================================================
# RELATÓRIO CONSOLIDADO PARA INTERPRETAÇÃO
# ==========================================================

cat("\n************************************************************\n")
cat("RESUMO PARA INTERPRETAÇÃO TÉCNICA\n")
cat("************************************************************\n")

cat("\n--- 1. RESULTADOS: DADOS LINEARES ---\n")
cat("R2 Ajustado:", round(res_linear$R2adj * 100, 2), "%\n")
cat("Significância do Modelo (p-valor):", res_linear$ani_model$`Pr(>F)`[1], "\n")
cat("Variáveis Significativas (p < 0.1):\n")
print(subset(as.data.frame(res_linear$ani_terms_full), `Pr(>F)` < 0.1))

cat("\n--- 2. RESULTADOS: DADOS OUTLINE (FORMA) ---\n")
cat("R2 Ajustado:", round(res_outline$R2adj * 100, 2), "%\n")
cat("Significância do Modelo (p-valor):", res_outline$ani_model$`Pr(>F)`[1], "\n")
cat("Variáveis Significativas (p < 0.1):\n")
print(subset(as.data.frame(res_outline$ani_terms_full), `Pr(>F)` < 0.1))

cat("\n--- 3. TESTE DE PROCRUSTES (CONCORDÂNCIA) ---\n")
cat("Correlação (r):", round(proc_model$t0, 4), "\n")
cat("Soma dos Quadrados (m2):", round(proc_model$ss, 4), "\n")
cat("Significância (p-valor):", proc_model$signif, "\n")

cat("\n************************************************************\n")
cat("FIM DO RELATÓRIO\n")
cat("************************************************************\n")

# ==========================================================
# FAZER UMA UNICA FIGURA
# ==========================================================

# 1. Preparar os gráficos A e B 
# (Damos o mesmo nome de escala para garantir que o patchwork os reconheça como iguais)
p_a <- plot_linear + 
  labs(title = "A. dbRDA: Dados Lineares") +
  theme(legend.position = "bottom")

p_b <- plot_outline + 
  labs(title = "B. dbRDA: Dados de Forma") +
  theme(legend.position = "bottom")

# 2. Preparar o gráfico C
p_c <- p_proc + 
  labs(title = "C. Correlação de Procrustes") +
  theme(legend.position = "bottom")

# 3. Montar a composição final
# (p_a | p_b) coloca os dois primeiros lado a lado
# / p_c coloca o procrustes embaixo
# plot_layout(guides = "collect") unifica as legendas repetidas de A e B
figura_final <- (p_a | p_b) / p_c + 
  plot_layout(guides = "collect", heights = c(1, 1)) & 
  theme(legend.position = "bottom")

# Exibir e Salvar
print(figura_final)

ggsave("Figura_Integrada_Final.png", figura_final,
       width = 22, height = 24, units = "cm", dpi = 600)
