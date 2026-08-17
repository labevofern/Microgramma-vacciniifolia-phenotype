# ==============================================================
# Análise morfométrica de contornos
# ==============================================================
# Autor: Juliana Aljahara
# Programa: Doutoranda em Biologia Vegetal
# Instituição: Universidade Federal de Pernambuco (UFPE)
# Orientadora: Dra. Thaís E. Almeida
# ==========================================================
# 1. Pacotes
# ==========================================================

library(readxl)        # Importação de planilhas Excel (.xlsx)
library(Momocs)        # Análises morfométricas de contorno (EFA, PCA, LDA)
library(tidyverse)     # Manipulação e organização de dados (dplyr, tidyr, ggplot2, etc.)
library(ggplot2)       # Construção de gráficos
library(tidyr)         # Organização e reestruturação de dados
library(dplyr)         # Manipulação de dados (filter, group_by, summarise)
library(ggrepel)       # Melhor posicionamento automático de textos em gráficos
library(ggExtra)       # Adição de gráficos marginais (ex.: boxplots na PCA)
library(vegan)         # Análises ecológicas e estatísticas multivariadas (NMDS, PERMANOVA, dbRDA, Procrustes)
library(sf)            # Manipulação e visualização de dados espaciais (Simple Features)
library(rnaturalearth) # Acesso a mapas-múndi e dados geográficos vetoriais
library(writexl)       # Exportação de tabelas e resultados para arquivos Excel (.xlsx)
library(pheatmap)      # Construção de mapas de calor (heatmaps) e matrizes de confusão
library(RColorBrewer)  # Paletas de cores para gráficos e visualizações

# ==========================================================
# 2. Definir diretório principal
# ==========================================================

setwd("C:/Users/nikso/OneDrive/_JULIANA ALJAHARA/1_Análise de macromorfologica/Diretório_Outline")

# ==========================================================
# 3. Definir as cores
# ==========================================================

my_colors <- c("#6A0572", "#AB83A1", "#F67280", "#F8B195")

# ==========================================================
# 4. Processamento inicial dos dados
# ==========================================================

# Listar arquivos de texto na pasta atual
lf <- list.files("Coordenadas_fourier/",pattern = "\\.txt$", full.names = TRUE)

# Visualizar
print(lf)

# Colocar nome do arquivo excel que to trabalhando
# Ler o arquivo Excel com os nomes das espécies e identificadores dos indivíduos
spp_names <- read_xlsx("Coordenadas_fourier/Contornos.xlsx")

# Exibir os nomes únicos das espécies na coluna 'spp' do arquivo Excel
unique(spp_names$`Dominio fitogeografico`)

# Importar as coordenadas dos arquivos de texto para um objeto de dados
lf_coord <- import_txt(lf)

# Verificar se o número de arquivos e linhas na planilha são iguais
stopifnot(length(lf) == nrow(spp_names))

# Renomear arquivos com base no identificador
names(lf) <- spp_names$Individuo

# Atribuir os nomes das espécies aos objetos
names(lf_coord) <- spp_names$`Dominio fitogeografico`

# Criar data frame com fatores
lf_fac <- data.frame(
  Type = spp_names$`Dominio fitogeografico`,
  Heb = spp_names$Herbario
)

# Criar objeto do Momocs
lf_out <- Out(lf_coord, fac = lf_fac)
lf_out$fac$Type <- as.factor(lf_out$fac$Type)

# ==========================================================
# 5. Forma dos indivíduos (Ajustado)
# ==========================================================

# 1. Forçar a reordenação explicitamente no slot $fac
ordem_desejada <- c("Atlantic Forest", "Caatinga", "Cerrado", "Pampa")

# Certifique-se de que está alterando a coluna 'Type' dentro de $fac
lf_out$fac$Type <- factor(lf_out$fac$Type, levels = ordem_desejada)

# 2. IMPORTANTE: No Momocs, às vezes é necessário reordenar as linhas do objeto 
# para que o gráfico 'panel' siga a ordem visual desejada
lf_out <- slice(lf_out, order(lf_out$fac$Type))

# 3. Iniciar o dispositivo gráfico
png("2_Tabelas_e_Imagens/formas_.png", width = 10, height = 8, units = "in", res = 300)

# 4. Gerar o gráfico
# O panel() agora lerá os dados já ordenados pelo slice
formas_gerais <- panel(lf_out, 
                       fac = "Type", 
                       names = "",
                       col = my_colors)

# 5. Adiciona a legenda (usando os níveis ordenados)
legend("topright", 
       legend = levels(lf_out$fac$Type),
       col = my_colors[1:length(unique(lf_out$fac$Type))],
       pch = 16,
       pt.cex = 2.5,
       cex = 0.8,
       bty = "n",
       inset = c(0, 0),
       xpd = TRUE,
       y.intersp = 2)

dev.off()

# ==========================================================
# 5. Definir o número ideal de harmonicos
# ==========================================================

# Determine optimal number of harmonics
harmonics_info <- calibrate_harmonicpower_efourier(lf_out, 
                                                   nb.h = 20, 
                                                   drop = 1,
                                                   thresh = c(80, 90, 95, 99, 99.9),
                                                   plot = TRUE)

# 1. Extrair a matriz de dados brutos
matriz_dados <- harmonics_info$q

# 2. Calcular a média de potência cumulativa para cada harmônico
# (Transformamos em data frame, calculamos a média das colunas e organizamos)
dados_grafico <- as.data.frame(matriz_dados) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  pivot_longer(cols = everything(), 
               names_to = "harmonics", 
               values_to = "cum") %>%
  mutate(
    # Extrai o número do harmônico (ex: "h1" vira 1, "h9" vira 9)
    h = as.numeric(gsub("h", "", harmonics))
  )

# 3. Construir o gráfico customizado com ggplot2 (Versão em Inglês)
p <- ggplot(dados_grafico, aes(x = h, y = cum)) +
  # Linhas de referência para os limiares (95% e 99%)
  geom_hline(yintercept = c(95, 99), linetype = "dashed", color = "red", alpha = 0.6) +
  # Linha conectando as médias dos harmônicos
  geom_line(color = "#2c3e50", linewidth = 1) +
  # Pontos para destacar cada harmônico individualmente
  geom_point(color = "#16a085", size = 3) +
  
  # === ALTERAÇÃO AQUI: Adicionar rótulos inclinados ===
  geom_text(aes(label = paste0(round(cum, 1), "%")), 
            angle = 45,       # Inclina o texto em 45 graus
            vjust = -0.5,     # Ajusta para cima do ponto
            hjust = -0.1,     # Ajusta levemente para a direita
            size = 3.5, 
            color = "#34495e") +
  
  # Configuração dos eixos
  scale_x_continuous(breaks = 1:max(dados_grafico$h)) +
  scale_y_continuous(limits = c(min(dados_grafico$cum) - 5, 105), breaks = seq(0, 100, 10)) +
  # Títulos e rótulos traduzidos para o inglês
  labs(
    title = "Mean Harmonic Power Calibration",
    subtitle = "Dashed lines indicate 95% and 99% thresholds of cumulative power",
    x = "Number of Harmonics",
    y = "Mean Cumulative Harmonic Power (%)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(face = "italic", hjust = 0.5, color = "gray30"),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold")
  )

# 4. Exibir o gráfico na tela do RStudio
print(p)

# 5. Salvar o gráfico na pasta
ggsave("2_Tabelas_e_Imagens/harmonicos_ggplot2.png", 
       plot = p, 
       width = 10, 
       height = 8, 
       units = "in", 
       dpi = 300)

# ==========================================================
# 5.1 Efourier análise
# ==========================================================

# Realizando a análise de Fourier elíptica (EFA) nas formas armazenadas em lf_out
lf_fou <- efourier(x = lf_out,    #Dados das coordenadas das formas (contornos) 
                   nb.h = 8,     #Quantos harmonicos o comando anterior falou?
                   norm = TRUE) #TRUE = considera apenas a forma/  FALSE = considera tamanho e posição

# ==========================================================
# 6. PCA
# ==========================================================

# Realizar a PCA
lf_pca <- PCA(x = lf_fou,
              fac = lf_fou1$fac)

criar_plot_pca <- function(pca_obj, out_obj) {
  
  pca_data <- data.frame(
    PC1 = pca_obj$x[, 1],
    PC2 = pca_obj$x[, 2],
    especie = out_obj$fac$Type
  )
  
  # % variância explicada
  variance_explained <- round((pca_obj$sdev^2 / sum(pca_obj$sdev^2)) * 100, 1)
  
  p <- ggplot(pca_data,
              aes(x = PC1,
                  y = PC2,
                  color = especie,
                  shape = especie)) +   # 👈 AQUI
    
    geom_point(size = 3.5, alpha = 0.85) +
    
    scale_color_manual(values = my_colors) +
    scale_shape_manual(values = c(17, 15, 16, 7)) +  # 👈 AQUI
    
    labs(
      x = paste0("PC1 (", variance_explained[1], "%)"),
      y = paste0("PC2 (", variance_explained[2], "%)")
    ) +
    
    theme_bw(base_size = 14) +
    
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      legend.direction = "horizontal",
      axis.title = element_text(face = "bold"),
      panel.grid = element_blank()
    )
  
  p_marginal <- ggMarginal(
    p,
    type = "boxplot",
    groupColour = TRUE,
    groupFill = TRUE,
    size = 6
  )
  
  return(p_marginal)
}

# Criar gráfico
p1 <- criar_plot_pca(lf_pca, lf_out)

# Exibir
print(p1)

# Salvar
ggsave("2_Tabelas_e_Imagens/PCA_.png",
       plot = p1,
       width = 10,
       height = 8,
       dpi = 300)

# Criar data frame com scores da PCA
scores_pca <- as.data.frame(lf_pca$x)

# Opcional: arredondar (bom para artigo)
scores_pca <- round(scores_pca, 5)

# Juntar com todas as colunas do spp_names
tabela_scores_completa <- cbind(spp_names, scores_pca)

# Visualizar
head(tabela_scores_completa)

# Salvar em formato Excel (.xlsx)
write_xlsx(tabela_scores_completa, 
           path = "2_Tabelas_e_Imagens/Scores PCA_Outline.xlsx")

# ==========================================================
# TABELA 1 - TODOS OS INDIVÍDUOS (FRONDES) - testar variação intraindíviduo
# ==========================================================

summary(lf_pca)

# Selecionar colunas PC
pcs_cols <- c("PC1", "PC2", "PC3")

# Criar tabela final
tabela_frondes <- tabela_scores_completa[, c("Herbario", "Ordem", pcs_cols)]

# Arredondar PCs (se ainda não estiverem)
tabela_frondes[pcs_cols] <- round(tabela_frondes[pcs_cols], 5)

# Salvar
write.csv(tabela_frondes,
          "2_Tabelas_e_Imagens/Tabela_PCA_todas_as_frondes.csv",
          row.names = FALSE)

# ==========================================================
# TABELA 2 - MÉDIA POR INDIVÍDUO - testar variação entreindividuos
# ==========================================================

# Selecionar apenas PCs
pcs_cols <- c("PC1", "PC2", "PC3")

tabela_media <- tabela_scores_completa %>%
  group_by(Herbario) %>%
  summarise(across(all_of(pcs_cols), mean))

# Arredondar
tabela_media[pcs_cols] <- round(tabela_media[pcs_cols], 5)

# Salvar em formato Excel (.xlsx)
write_xlsx(tabela_media, 
           path = "2_Tabelas_e_Imagens/Scores PCA_Outline_média.xlsx")

# ==========================================================
# 7. LDA
# ==========================================================

# Dados já existentes (supondo que lf_fou já foi criado antes da PCA)
coefs <- lf_fou$coe           # Extrai coeficientes
fac_data <- lf_fou$fac        # Extrai fatores (Type e Heb)

# Calcular MÉDIAS por Heb (indivíduo) -  EU CONSIDERO A MÉDIA PARA EVITAR PSEUDOREPETIÇÃO DENTRO DA LDA
lf_fou_avg <- data.frame(coefs, fac_data) %>%  
  group_by(Heb, Type) %>%  
  summarise(across(where(is.numeric), mean), .groups = "drop")

# Remover as colunas com valores 0 e 1
lf_fou_avg <- lf_fou_avg[, -c(3, 11, 19)]

# Remover colunas constantes (harmônicos que não variam)
# O Momocs gera algumas colunas 0 ou 1 dependendo da normalização
lf_fou_avg_clean <- lf_fou_avg[, -c(1:2)] # Remove Heb e Type temporariamente
lf_fou_avg_clean <- lf_fou_avg_clean[, apply(lf_fou_avg_clean, 2, var) > 0]

# 2. Rodar a LDA (usando a função LDA do Momocs ou MASS::lda)
# Aqui usamos a estrutura do Momocs para manter a compatibilidade
lda_data1 <- LDA(as.matrix(lf_fou_avg_clean), fac = lf_fou_avg$Type)

# 3. Extrair dados para o ggplot2
lda_scores <- as.data.frame(lda_data1$mod.pred$x)

# SOLUÇÃO DEFINITIVA: Forçar a coluna 'Type' vindo direto do dado original
lda_scores$Type <- lf_fou_avg$Type

# 4. Cálculo da Acurácia
acc_geral <- mean(lda_data1$CV.correct)
label_acc <- paste0("Accuracy: ", round(acc_geral * 100, 2), "%")

# 5. Visualizar a Matriz de Confusão
matriz_confusao <- lda_data1$CV.tab
print(matriz_confusao)

# Para uma visualização mais detalhada com as porcentagens de acerto por classe:
# Isso mostra quanto por cento de cada grupo foi classificado corretamente
matriz_prop <- prop.table(matriz_confusao, 1) * 100
print(round(matriz_prop, 2))

# 5. Construção do Gráfico
if ("LD2" %in% colnames(lda_scores)) { 
  
  p_lda <- ggplot(lda_scores, aes(x = LD1, y = LD2, fill = Type, color = Type)) +
    geom_point(shape = 21, size = 5, color = "black", alpha = 0.8) +
    stat_ellipse(show.legend = FALSE, linetype = 2) +
    scale_fill_manual(values = my_colors) +
    scale_color_manual(values = my_colors) +
    theme_classic(base_size = 14) +
    labs(title = "",
         subtitle = label_acc,
         x = "LD1",
         y = "LD2",
         fill = "Domínio",
         color = "Domínio") +
    theme(legend.position = "bottom",
          panel.grid = element_blank(),
          axis.title = element_text(face = "bold"))
  
} else { 
  
  p_lda <- ggplot(lda_scores, aes(x = LD1, fill = Type)) +
    geom_density(alpha = 0.7, color = "black") +
    scale_fill_manual(values = my_colors) +
    theme_classic(base_size = 14) +
    labs(title = "",
         subtitle = label_acc,
         x = "LD1",
         y = "Densidade",
         fill = "Domínio") +
    theme(legend.position = "bottom")
}

# 6. Exibir
print(p_lda)

ggsave("2_Tabelas_e_Imagens/LDA_.png",
       plot = p_lda,
       width = 10,
       height = 8,
       dpi = 300)

# ==============================================================================
# Extrair dados para o ggplot2 e criar DataFrame definitivo de Scores
# ==============================================================================

# Extrair os scores brutos da LDA (LD1, LD2, etc.)
lda_scores <- as.data.frame(lda_data1$mod.pred$x)

# Adicionar as colunas identificadoras vindas do dataframe agrupado (lf_fou_avg)
lda_scores <- lda_scores %>%
  mutate(
    Heb  = lf_fou_avg$Heb,   # Adiciona a coluna do indivíduo/herbário
    Type = lf_fou_avg$Type   # Adiciona a coluna do tipo/domínio
  ) %>%
  # Organiza o DataFrame colocando os identificadores no início
  select(Heb, Type, everything())

# Exibir as primeiras linhas no console para conferir
head(lda_scores)

# Salvar a tabela de scores da LDA do Momocs
write.csv(lda_scores, 
          file = "2_Tabelas_e_Imagens/Scores_LDA_Momocs.csv", 
          row.names = FALSE)

# ==========================================================
# 8. Média das formas
# ==========================================================

# Definir o nome do arquivo e a resolução (300 DPI)
png("2_Tabelas_e_Imagens/formas_médiascv.png", width = 12, height = 8, units = "in", res = 300)

# Definir as cores
my_colors_shape <- c("#6A0572","#F67280")

# Calcular as médias das formas
lf_mean_shape <- MSHAPES(lf_fou, ~Type)

# Plotar usando plot_MSHAPES (sem o argumento 'main')
plot_MSHAPES(lf_mean_shape, palette = pal_manual(my_colors_shape, transp = 0))

# Fechar o dispositivo gráfico e salvar o arquivo
dev.off()

# ==============================================================================
# 9. PERMANOVA (Teste de diferença na forma foliar entre Domínios)
# ==============================================================================

set.seed(123)

# Executar a PERMANOVA
# Usamos scale() nos coeficientes para padronizar as magnitudes dos harmônicos
# e testamos em função da variável 'Type' (Domínios)
permanova_efa <- adonis2(scale(lf_fou_avg[, -c(1:2)]) ~ Type, 
                         data = lf_fou_avg, 
                         method = "euclidean",
                         permutations = 9999)

print("=== Resultados da PERMANOVA (Forma dos Contornos) ===")
print(permanova_efa)

# Opcional: Salvar o resultado da PERMANOVA em um arquivo de texto para usar no artigo
capture.output(permanova_efa, file = "2_Tabelas_e_Imagens/Resultados_PERMANOVA_EFA.txt")

# ==============================================================================
# NOVO: MATRIZ DE DISTÂNCIA (Baseada nos Coeficientes de Fourier Médios)
# ==============================================================================

# 1. Preparar os dados
# Usamos o 'lf_fou_avg' que você já criou (médias por indivíduo/Herbario)
# O scale() aqui é importante se os harmônicos tiverem variações de magnitude muito distintas
efa_scaled <- scale(lf_fou_avg[, -c(1:2)]) # Remove as colunas 'Heb' e 'Type'
rownames(efa_scaled) <- lf_fou_avg$Heb

# 2. Cálculo da Matriz de Distância Euclidiana
# Esta distância representa o quão "longe" uma forma média está da outra no espaço morfométrico
dist_matrix_efa <- dist(efa_scaled, method = "euclidean")

# Converter em matriz para salvar
dist_table_efa <- as.matrix(dist_matrix_efa)

# 3. Salvar Matriz
write.csv(dist_table_efa, "2_Tabelas_e_Imagens/Dados outline_Matriz_Distancia_Euclidiana.csv")

# --- 1. Preparação dos Metadados para o Heatmap ---
# Criamos um dataframe de anotação usando as médias por herbário
annotation_data <- data.frame(Bioma = lf_fou_avg$Type)
rownames(annotation_data) <- lf_fou_avg$Heb

# --- 2. Configuração da Identidade Visual (Cores Exatas) ---
ann_colors = list(
  Bioma = c(
    "Atlantic Forest" = "#6A0572", 
    "Cerrado"         = "#F67280", 
    "Caatinga"        = "#AB83A1", 
    "Pampa"           = "#F8B195"
  )
)

# --- 3. Geração do Heatmap (Mesmo estilo do script anterior) ---
# Usamos a dist_table_efa gerada no passo anterior
pheatmap(dist_table_efa, 
         clustering_method = "ward.D2", 
         annotation_row = annotation_data, 
         annotation_col = annotation_data, 
         annotation_colors = ann_colors,
         
         # Escala Roxo: Branco (Próximo) -> Escuro (Distante)
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
         filename = "2_Tabelas_e_Imagens/Heatmap_Distancia_EFA_Final.png",
         width = 15,   
         height = 10,  
         res = 300
)

# ==========================================================
# 10. MAPA FINAL REFINADO: 6 RÉPLICAS POR BIOMA - DESIGN AVANÇADO
# ==========================================================

library(ggplot2)
library(sf)
library(rnaturalearth)
library(dplyr)
# library(Momocs) # Assumindo que Momocs e lf_out já estão carregados/processados

# 1. Preparação Básica e Cores
brazil <- ne_countries(scale = "medium", returnclass = "sf", country = "brazil")
# Cores vibrantes e distintas
my_colors <- c("#6A0572", "#AB83A1", "#F67280", "#F8B195") 

# Sincronizar nomes (Garantindo que spp_names existe e está correto)
spp_names <- spp_names %>%
  dplyr::rename(
    individuo = 1, 
    dominio_fitogeografico = 2,
    lat = matches("lat", ignore.case = TRUE), # Pega 'Lat' ou 'lat'
    log = matches("log|lon", ignore.case = TRUE) # Pega 'Log', 'log', 'Long'...
  )

# 2. Seleção de 6 exemplares espacialmente distintos por bioma (Total 24)
exemplares_replicas <- spp_names %>%
  group_by(dominio_fitogeografico) %>%
  group_modify(~ {
    # K-means para encontrar 6 centros geográficos por bioma
    set.seed(123) # Para reprodutibilidade
    km <- kmeans(cbind(.x$log, .x$lat), centers = 6)
    .x %>%
      mutate(cluster = km$cluster) %>%
      group_by(cluster) %>%
      slice_sample(n = 1) %>% 
      ungroup()
  }) %>%
  mutate(replica_id = row_number()) %>% 
  ungroup()

indices_replicas <- match(exemplares_replicas$individuo, spp_names$individuo)
n_total <- nrow(exemplares_replicas) # Deve ser 24 (4 biomas x 6 pontos)

# 3. CÁLCULO DE POSIÇÕES RADIAIS (24 Posições no Total)
raio_x <- 40 # Distância horizontal ligeiramente aumentada
raio_y <- 34 # Distância vertical ligeiramente aumentada
centro_log <- -55
centro_lat <- -15

# Ângulos distribuídos para as 24 folhas
angulos <- seq(0, 2*pi, length.out = n_total + 1)[1:n_total]

pos_replicas <- exemplares_replicas %>%
  mutate(
    angulo = angulos,
    dest_log = centro_log + raio_x * cos(angulo),
    dest_lat = centro_lat + raio_y * sin(angulo)
  )

# 4. Criar o Mapa Base Refinado
p_replicas <- ggplot() +
  # Mapa do Brasil com linhas BEM mais marcadas
  geom_sf(data = brazil, fill = "#FDFDFD", color = "gray60", linewidth = 0.8) +
  # Pontos de amostragem total (fundo discreto, mas visível)
  geom_point(data = spp_names, aes(x = log, y = lat, color = dominio_fitogeografico), 
             alpha = 0.08, size = 1.2) +
  scale_color_manual(values = my_colors) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    panel.background = element_rect(fill = "white", color = NA) # Fundo limpo
  )

# 5. Loop de Construção com Novas Estéticas (Sombra e Contorno)
for(i in 1:nrow(exemplares_replicas)) {
  bioma_nome <- exemplares_replicas$dominio_fitogeografico[i]
  origem_log <- exemplares_replicas$log[i]
  origem_lat <- exemplares_replicas$lat[i]
  pos <- pos_replicas[i, ]
  
  # Assumindo Momocs e coo_* funcionais
  folha_coord <- lf_out$coo[[indices_replicas[i]]]
  folha_coord <- coo_center(folha_coord)
  folha_coord <- coo_scale(folha_coord) 
  folha_df <- as.data.frame(folha_coord)
  colnames(folha_df) <- c("x", "y")
  
  # Rotação
  folha_df_rot <- data.frame(x = folha_df$y, y = -folha_df$x)
  
  # AJUSTE DE TAMANHO VISUAL E POSIÇÃO
  fator_visual <- 2.8  
  x_final <- (folha_df_rot$x * fator_visual) + pos$dest_log
  y_final <- (folha_df_rot$y * fator_visual) + pos$dest_lat
  
  folha_final_df <- data.frame(x = x_final, y = y_final)
  
  # Deslocamento para a sombra (drop shadow)
  shadow_offset_x <- 0.5
  shadow_offset_y <- -0.5
  folha_shadow_df <- data.frame(x = x_final + shadow_offset_x, 
                                y = y_final + shadow_offset_y)
  
  cor_atual <- my_colors[match(bioma_nome, sort(unique(spp_names$dominio_fitogeografico)))]
  
  # ORDEM DE DESENHO: Linha -> Sombra -> Polígono
  p_replicas <- p_replicas + 
    # 1. Linha TRACEJADA (linetype = "dashed")
    annotate("segment", 
             x = origem_log, y = origem_lat, 
             xend = pos$dest_log, yend = pos$dest_lat, 
             color = cor_atual, 
             linewidth = 0.45, 
             alpha = 0.5, # Suavizada ligeiramente
             linetype = "dashed") +
    
    # 2. Ponto de origem
    annotate("point", x = origem_log, y = origem_lat, color = cor_atual, size = 1.4) +
    
    # 3. POLÍGONO DA SOMBRA (Efeito de profundidade)
    geom_polygon(data = folha_shadow_df, aes(x = x, y = y), 
                 fill = "black", 
                 alpha = 0.15, # Muito suave
                 color = NA) +
    
    # 4. POLÍGONO DA FOLHA (Alpha = 1, contorno mais marcado)
    geom_polygon(data = folha_final_df, aes(x = x, y = y), 
                 fill = cor_atual, 
                 color = "gray10", # Contorno preto suave
                 linewidth = 0.35, # Ligeiramente mais grosso
                 alpha = 1) 
}

# 6. Finalização e Enquadramento (Ajustado para o novo raio)
p_replicas <- p_replicas + 
  coord_sf(xlim = c(-110, 0), ylim = c(-65, 35))

print(p_replicas)

# 7. Salvar
ggsave("2_Tabelas_e_Imagens/mapa_final_limpo.png", p_replicas, width = 12, height = 10, dpi = 300)






