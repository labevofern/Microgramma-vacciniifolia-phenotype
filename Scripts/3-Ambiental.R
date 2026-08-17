# ==============================================================
# Extração de variáveis ambientais e análises
# ==============================================================
# Autor: Juliana Aljahara
# Programa: Doutoranda em Biologia Vegetal
# Instituição: Universidade Federal de Pernambuco (UFPE)
# Orientadora: Dra. Thaís E. Almeida
# ==========================================================
# 1. Pacotes
# ==========================================================

library(terra)        # Manipulação e análise de dados espaciais, incluindo rasters de clima, solo e outras variáveis ambientais
library(dplyr)        # Manipulação, filtragem e organização de tabelas de dados (data.frames)
library(usdm)         # Diagnóstico de multicolinearidade em dados ecológicos (funções vif e vifstep)
library(readxl)       # Importação de planilhas do Excel (.xls, .xlsx)
library(vegan)        # Análises ecológicas e estatísticas multivariadas (NMDS, PERMANOVA, dbRDA, Procrustes)
library(ggplot2)      # Criação de gráficos estáticos baseados na gramática de gráficos
library(missForest)   # Imputação de dados faltantes (NA) utilizando algoritmos de Random Forest
library(writexl)      # Exportação de tabelas e resultados para arquivos Excel (.xlsx)
library(geobr)        # Download e manipulação de dados espaciais oficiais do Brasil (limites, biomas, estados, municípios etc.)
library(factoextra)   # Visualização e interpretação de análises multivariadas, como PCA e agrupamentos

# ==========================================================
# 2. CONFIGURAÇÕES DE DIRETÓRIO
# ==========================================================
setwd("C:/Users/nikso/OneDrive/_JULIANA ALJAHARA/1_Análise de macromorfologica/Diretório_Ambiental")

# ==========================================================
# 3. LEITURA E NOMEAÇÃO DOS RASTERS
# ==========================================================
arquivos_tif <- list.files(path = ".", pattern = "\\.tif$", recursive = TRUE, full.names = TRUE)
lista_rasters <- lapply(arquivos_tif, rast)

# Atribuir nomes às camadas baseados nos arquivos originais
nomes <- tools::file_path_sans_ext(basename(arquivos_tif))
for(i in 1:length(lista_rasters)) names(lista_rasters[[i]]) <- nomes[i]

# --- Conferência dos Dados ---
cat("Total de matrizes ambientais:", length(lista_rasters), "\n")
cat("Nomes das matrizes:\n")
print(nomes)

# ==========================================================
# 4. PADRONIZAÇÃO ESPACIAL (30 arc-seconds ~ 1 km)
# ==========================================================

# 1. Criar Raster de Referência (Resolução de 30 arc-seconds ~ 0.00833 graus)
# Definimos a extensão do Brasil para garantir que todos os rasters tenham o mesmo "grid"
ref_res <- 0.008333333 
ext_brasil <- ext(-74, -34, -34, 5) # Longitude Oeste/Leste, Latitude Sul/Norte
raster_ref <- rast(extent = ext_brasil, resolution = ref_res, crs = "EPSG:4326")

# 2. Padronizar e Recortar as Matrizes (Resampling + Crop)
# O resample já vai "encaixar" seus dados na extensão e resolução do raster_ref
cat("Iniciando a padronização e recorte para a extensão do Brasil...\n")

lista_final <- lapply(lista_rasters, function(r) {
  # bilinear para variáveis contínuas (clima, altitude, etc)
  # Se tiver dados categóricos (uso do solo), mude para method = "near"
  r_resampled <- resample(r, raster_ref, method = "bilinear")
  
  # Garantindo que o recorte seja exato à extensão definida
  r_crop <- crop(r_resampled, raster_ref)
  
  return(r_crop)
})

cat("Processo concluído. Todas as matrizes estão no mesmo grid e extensão.\n")

# ==========================================================
# 5. CARREGAMENTO DOS PONTOS E EXTRAÇÃO
# ==========================================================

# 1. Carregar os pontos
dados <- read_excel("2_Tabelas_e_Imagens/Contornos.xlsx")

# 2. Converter e Reprojetar
meu_crs <- crs(lista_final[[1]])
pontos_spat <- vect(dados, geom = c("Log", "Lat"), crs = "EPSG:4326")
pontos_spat <- project(pontos_spat, meu_crs)

# 3. CRIAR O BUFFER (Exemplo: 1000 metros = 1km)
# O valor 'width' depende da unidade do seu CRS. 
# Se for Geográfico (graus), use valores como 0.00833 para ~1km.
# Se for Projetado (metros), use 1000 para 1km.
pontos_buffer <- buffer(pontos_spat, width = 0.016666) #2KM

# 3.5. INICIALIZAR O OBJETO
# Criamos uma cópia dos dados originais para receber as novas colunas
dados_completos <- dados 

# 4. Na extração (Loop)
for(i in 1:length(lista_final)) {
  nome_var <- nomes[i]
  
  # Extração com a média do buffer
  extraido <- terra::extract(lista_final[[i]], pontos_buffer, fun = mean, na.rm = TRUE, ID = FALSE)
  
  # Agora o objeto 'dados_completos' existe e pode receber a nova coluna
  dados_completos[[nome_var]] <- extraido[, 1]
}

# -----------------FILTRO DE REPETIÇÕES: Manter apenas um registro por 'Herbario'

# Removendo as linhas duplicadas com base na coluna 'Herbario'
dados_completos <- dados_completos %>%
  distinct(Herbario, .keep_all = TRUE)

# Verificando se as duplicatas sumiram
nrow(dados_completos)

# -----------------ELIMINAR VARIAVEIS AMBIENTAIS COM MUITOS NAS (>15%)

# Calcular a porcentagem de NA por variável (da col 11 em diante)
colunas_ambientais <- dados_completos[, 11:ncol(dados_completos)]
porcentagem_na <- colMeans(is.na(colunas_ambientais)) * 100

# Ver quais variáveis estão acima de um limite (ex: 15%)
vars_ruins <- names(porcentagem_na[porcentagem_na > 15])
cat("Variáveis com excesso de NAs (>15%):", vars_ruins, "\n")

# Manter apenas as colunas boas
dados_completos <- dados_completos[, !(names(dados_completos) %in% vars_ruins)]

# ==========================================================
# 6. CORREÇÃO DE VALORES NA (VIZINHOS MAIS PRÓXIMOS)
# ==========================================================

# Identificar quais colunas receberam os dados ambientais (a partir da coluna 11)
colunas_ambientais <- nomes # Usa a lista de nomes que você já tem

for(i in 1:length(lista_final)) {
  nome_var <- nomes[i]
  raster_atual <- lista_final[[i]]
  
  # 1. Identificar quais linhas ainda estão com NA para a variável atual
  linhas_na <- which(is.na(dados_completos[[nome_var]]))
  
  if(length(linhas_na) > 0) {
    cat(paste0("Corrigindo ", length(linhas_na), " pontos com NA na variável: ", nome_var, "\n"))
    
    # 2. Filtrar os pontos geométricos que deram problema
    pontos_com_na <- pontos_spat[linhas_na]
    
    # 3. Buscar os 2 pixels válidos (não-NA) mais próximos de cada ponto
    # A função terra::distance com 'pairs=FALSE' e fornecendo o raster nos dá os valores correspondentes
    valores_proximos <- terra::extract(
      raster_atual, 
      pontos_com_na, 
      method = "simple", # busca simples baseada na geometria
      xy = FALSE
    )
    
    # Caso o extract simples falhe por estar fora do grid, usamos a função 'near' ou 'distance'
    # Alternativa robusta do pacote terra para vizinhos mais próximos:
    for(j in 1:length(linhas_na)) {
      idx_linha_original <- linhas_na[j]
      ponto_atual <- pontos_com_na[j]
      
      # Transforma o raster em pontos apenas onde há dados (remover NAs do raster de busca)
      # Para não pesar o código, limitamos a busca a uma janela ao redor do ponto
      buffer_busca <- terra::buffer(ponto_atual, width = 0.1) # Janela de busca (~10km)
      raster_cortado <- terra::crop(raster_atual, buffer_busca)
      pontos_raster <- terra::as.points(raster_cortado, values = TRUE, na.rm = TRUE)
      
      if(length(pontos_raster) > 0) {
        # Calcula a distância do nosso ponto para todos os pixels válidos próximos
        distancias <- terra::distance(ponto_atual, pontos_raster)
        
        # Pega o índice dos 2 mais próximos
        indices_2_proximos <- order(distancias)[1:2]
        
        # Extrai os valores desses 2 pixels e calcula a média
        valores_2_proximos <- pontos_raster[[nome_var]][indices_2_proximos, ]
        media_vizinhos <- mean(valores_2_proximos, na.rm = TRUE)
        
        # Atribui a média ao dataframe original
        dados_completos[[nome_var]][idx_linha_original] <- media_vizinhos
      }
    }
  }
}

# ==========================================================
# 5. SALVAR O RESULTADO FINAL ATUALIZADO
# ==========================================================
write.csv(dados_completos, "2_Tabelas_e_Imagens/Dados_Ocorrencia_Ambiental_tudo_semNA.csv", row.names = FALSE)

cat("\nExtração e correção de NAs concluídas com sucesso!\n")

# ==========================================================
#  SELEÇÃO DE VARIÁVEIS AMBIENTAIS COM VIF (th = 7)
# ==========================================================
# Ver a ordem atual das colunas a partir da 11
names(dados_completos)[11:ncol(dados_completos)]

colnames(dados_completos)[11:ncol(dados_completos)] <- c(
  # [1] a [11] Bioclimáticas (WorldClim) - Parte 1
  "Bio1", "Bio10", "Bio11", "Bio12", "Bio13", "Bio14", "Bio15", "Bio16", "Bio17", "Bio18", "Bio19",
  
  # [12] a [19] Bioclimáticas (WorldClim) - Parte 2
  "Bio2", "Bio3", "Bio4", "Bio5", "Bio6", "Bio7", "Bio8", "Bio9",
  
  # [20] a [26] Solo (SoilGrids)
  "Soil_cec",       # soilgrids-isric_myRegion_cec_5-15cm_mean
  "Soil_clay",      # soilgrids-isric_myRegion_clay_5-15cm_mean
  "Soil_nitrogen",  # soilgrids-isric_myRegion_nitrogen_5-15cm_mean
  "Soil_pH",        # soilgrids-isric_myRegion_phh2o_5-15cm_mean
  "Soil_sand",      # soilgrids-isric_myRegion_sand_5-15cm_mean
  "Soil_silt",      # soilgrids-isric_myRegion_silt_5-15cm_mean
  "Soil_soc",       # soilgrids-isric_myRegion_soc_5-15cm_mean
  
  # [27] a [31] Clima, Estrutura e Heterogeneidade - Parte 1
  "Aridity_index",         # ai_v31_yr           
  "Canopy_height",         # altura_copa_Brasil   
  "Texture_variability",   # cv_01_05_5km_uint32   
  
  # [32] a [35] Topografia e Heterogeneidade - Parte 2
  "Altitude",                # elevation_5KMmd_SRTM        
  "Evapotranspiration",      # et0_v31_yr                         
  "Texture_evenness",        # evenness_01_05_1km_uint16
  "Texture_amplitude",       # range_01_05_1km_uint16
  
  # [36] a [39] Topografia Final e Heterogeneidade - Parte 3
  "Terrain_roughness",       # terrain_roughness_index_5KMmd_SRTM   
  "Topo_position",           # topographic_position_index_5KMmd_SRTM  
  "Texture_variance"         # Variance_01_05_1km_uint32
)

# --- Verificação de segurança ---
cat("Total de colunas mapeadas: 39\n")
cat("Total de colunas no DataFrame:", length(11:ncol(dados_completos)), "\n")

# 1. Separar as variáveis ambientais e forçar para data.frame puro
variaveis_ambientais <- as.data.frame(dados_completos[, 11:ncol(dados_completos)])

# 2. Calcular VIF com threshold = 7
vif_resultado <- vifstep(variaveis_ambientais, th = 7)

# 3. Extrair variáveis mantidas
variaveis_mantidas <- vif_resultado@results$Variables

# 4. Criar nova matriz mantendo colunas 1:10 + variáveis selecionadas
nova_matriz <- cbind(
  dados_completos[, 1:10],
  variaveis_ambientais[, variaveis_mantidas]
)

# Salvar matriz de dados limpa
write.csv(nova_matriz, "2_Tabelas_e_Imagens/Especimes_e_Variaveis_Selecionadas_VIF.csv", row.names = FALSE)

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# TABELA COMPLETA PARA PUBLICAÇÃO: VIF INICIAL VS. STATUS FINAL
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# 1. Calcular o VIF inicial de todas (para mostrar a colinearidade original)
v_inicial <- vif(variaveis_ambientais)
colnames(v_inicial) <- c("Variables", "VIF_Inicial")

# 2. Pegar os resultados das que ficaram após o stepwise (VIF Final)
v_final <- vif_resultado@results
colnames(v_final) <- c("Variables", "VIF_Final")

# 3. Unir tudo em uma tabela "Mestra"
tabela_mestra <- v_inicial %>%
  left_join(v_final, by = "Variables") %>%
  mutate(
    `Selection Status` = ifelse(is.na(VIF_Final), "Excluída (VIF > 7)", "Selecionada"),
    # Se excluída, mostra o VIF alto inicial; se ficou, mostra o VIF final baixo
    `VIF Value` = ifelse(is.na(VIF_Final), VIF_Inicial, VIF_Final)
  )

# 4. Limpar, Arredondar e Ordenar
tabela_mestra$`VIF Value` <- round(tabela_mestra$`VIF Value`, 2)

tabela_mestra <- tabela_mestra %>%
  select(Variables, `VIF Value`, `Selection Status`) %>%
  arrange(desc(`Selection Status`), `VIF Value`)

# Renomear colunas para o padrão de Tabelas Suplementares
colnames(tabela_mestra) <- c("Environmental Variable", "VIF", "Selection Status")

# Salvar em formato Excel (.xlsx)
write_xlsx(tabela_mestra, 
           path = "2_Tabelas_e_Imagens/Tabela_Suplementar_VIF_Completa_Final.xlsx")

# Visualizar resultado no console
print(tabela_mestra)

# ==========================================================
# MATRIZ DE DISTÂNCIA PARA NMDS
# ==========================================================

setwd("C:/Users/nikso/OneDrive/_JULIANA ALJAHARA/1_Análise de macromorfologica/Diretório_Ambiental")

# Carregar dados
dados <- read.csv("2_Tabelas_e_Imagens/Especimes_e_Variaveis_Selecionadas_VIF.csv")

# ==========================================================
# 2. Preparação dos Dados
# ==========================================================

# Separar identificadores (colunas 1 a 10)
ids <- dados[, 1:10]
ids$Dominio.fitogeografico <- as.factor(ids$Dominio.fitogeografico)

# Separar variáveis ambientais (colunas 11 em diante)
variaveis_ambientais <- dados[, 11:ncol(dados)]

# Garantir que todas as colunas ambientais sejam numéricas
variaveis_ambientais <- as.data.frame(lapply(variaveis_ambientais, function(x) as.numeric(as.character(x))))

# ==========================================================
# 4. Padronização e Matriz de Distância
# ==========================================================

# Padronizar variáveis (Média 0, Variância 1) 
# Importante para Distância Euclidiana com variáveis em escalas diferentes!
variaveis_pad <- scale(variaveis_ambientais)

rownames(variaveis_pad) <- dados[, 5]

# Gerar Matriz de distância Euclidiana
dist_ambiental <- dist(variaveis_pad, method = "euclidean")

# Converter para dataframe para visualização ou exportação
df_dist <- as.data.frame(as.matrix(dist_ambiental))

# Salvar se necessário
write.csv(df_dist, "dist_ambiental.csv")

# ==========================================================
# 6. NMDS
# ==========================================================

# 1. Rodar NMDS com k = 5 para obter 5 dimensões
nmds_5d <- metaMDS(variaveis_pad, 
                   distance = "euclidean", 
                   k = 5, 
                   trymax = 100)

# ==========================================================
# EXTRAÇÃO DE SCORES E PLOTAGEM
# ==========================================================

# 2. Extrair as coordenadas (NMDS1 até NMDS5)
scores_nmds <- as.data.frame(scores(nmds_5d, display = "sites"))

# 3. Adicionar as colunas de identificação
scores_nmds$Dominio <- ids$Dominio.fitogeografico
scores_nmds$Herbario <- ids$Herbario

# 4. Salvar como Tabela Suplementar
write.csv(scores_nmds, "2_Tabelas_e_Imagens/Tabela_Suplementar_Scores_NMDS.csv", row.names = FALSE)

# Salvar em formato Excel (Ajustado o nome do arquivo)
write_xlsx(scores_nmds, path = "2_Tabelas_e_Imagens/Tabela_Suplementar_Scores_NMDS.xlsx")

# ==========================================================
# GRÁFICO NMDS
# ==========================================================

# Paleta personalizada
my_colors <- c("#6A0572", "#AB83A1", "#F67280", "#F8B195")

# Criar gráfico 
# Corrigido o erro: nmds$stress para nmds_5d$stress
grafico_nmds <- ggplot(scores_nmds, 
                       aes(x = NMDS1, 
                           y = NMDS2, 
                           color = Dominio, 
                           shape = Dominio)) +
  geom_point(size = 3, alpha = 0.8) + # Adicionado alpha para melhor visualização
  scale_color_manual(values = my_colors) +
  scale_shape_manual(values = c(17, 15, 16, 7)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  ) +
  labs(
    x = "NMDS1",
    y = "NMDS2",
    subtitle = paste("Stress =", round(nmds_5d$stress, 3)) # <--- Correção aqui
  )

# Exibir o gráfico
print(grafico_nmds)

ggsave("2_Tabelas_e_Imagens/Grafico_NMDS_Ambiental.png", 
       plot = grafico_nmds, 
       width = 8, height = 6, units = "in")

# ==========================================================
# 7. PERMANOVA (9.999 Permutações)
# ==========================================================
# Testa se os centros dos grupos (centroides) são diferentes

resultado_permanova <- adonis2(dist_ambiental ~ Dominio.fitogeografico,
                               data = ids,
                               permutations = 9999)

print(resultado_permanova)

# ==========================================================
# 8. TESTE DE DISPERSÃO (BETADISPER)
# ==========================================================

# 1. Calcular a dispersão (distância das amostras ao centroide do grupo)
bd <- betadisper(dist_ambiental, ids$Dominio.fitogeografico)

# 2. ANOVA Clássica
# Testa a diferença das médias das distâncias usando a distribuição F
anova_dispersao <- anova(bd)
print("--- ANOVA Clássica de Dispersão ---")
print(anova_dispersao)

# 3. Permutest (O que você pediu especificamente)
# Mais robusto para dados ecológicos/ambientais
perm_dispersao <- permutest(bd, permutations = 9999)
print("--- Teste de Permutação de Dispersão (9.999 perm) ---")
print(perm_dispersao)

# 4. Visualização Rápida (Opcional, mas ajuda a ver o espalhamento)
plot(bd, main = "Dispersão Multivariada por Domínio", ellipse = TRUE, hull = FALSE)

# ==========================================================
# SALVAMENTO DAS TABELAS SUPLEMENTARES
# ==========================================================

# 1. Salvar PERMANOVA
# Transformamos em dataframe e garantimos que o CSV inclua os nomes das colunas/linhas
tabela_permanova <- as.data.frame(resultado_permanova)
write.csv(tabela_permanova, 
          "2_Tabelas_e_Imagens/Tabela_Suplementar_PERMANOVA.csv", 
          row.names = TRUE)

# 2. Salvar BETADISPER (ANOVA)
tabela_betadisper_anova <- as.data.frame(anova_dispersao)
write.csv(tabela_betadisper_anova, 
          "2_Tabelas_e_Imagens/Tabela_Suplementar_Betadisper_ANOVA.csv", 
          row.names = TRUE)

# 3. Salvar BETADISPER (Permutest)
# Usamos o objeto 'perm_dispersao' que já criamos com 9.999 permutações
tabela_betadisper_perm <- as.data.frame(perm_dispersao$tab)
write.csv(tabela_betadisper_perm, 
          "2_Tabelas_e_Imagens/Tabela_Suplementar_Betadisper_permutest.csv", 
          row.names = TRUE)

print("Tabelas salvas com sucesso na pasta 2_Tabelas_e_Imagens!")

# ==============================================================================
# ANÁLISE MORFOMÉTRICA GEOGRÁFICA: LINEAR VS OUTLINE (Com Validação de Silhueta)
# ==============================================================================

# 1. CARREGAMENTO DE PACOTES
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse, cluster, factoextra, ggthemes, 
  rnaturalearth, sf, patchwork, vegan, viridis
)

# 2. DETERMINAÇÃO DO NÚMERO IDEAL DE CLUSTERS (Silhueta)
# ------------------------------------------------------------------------------
# Dica: O pico do gráfico indica o melhor K (mais coesão e separação)
# ------------------------------------------------------------------------------
dados <- read_excel("2_Tabelas_e_Imagens/Contornos.xlsx")

# Removendo as linhas duplicadas com base na coluna 'Herbario'
dados <- dados %>%
  distinct(Herbario, .keep_all = TRUE)

# Verificando se as duplicatas sumiram
nrow(dados)

dist_linear_df <- read.csv("2_Tabelas_e_Imagens/dist_linear.csv", row.names = 1, check.names = FALSE)
dist_outline_df <- read.csv("2_Tabelas_e_Imagens/dist_outline.csv", row.names = 1, check.names = FALSE)
dist_ambiental_df <- read.csv("2_Tabelas_e_Imagens/dist_ambiental.csv", row.names = 1, check.names = FALSE)

dist_linear_mat  <- as.matrix(dist_linear_df)
dist_outline_mat <- as.matrix(dist_outline_df)
dist_ambiental_mat <- as.matrix(dist_ambiental_df)

dist_linear  <- as.dist(dist_linear_mat)
dist_outline <- as.dist(dist_outline_mat)
dist_ambiental <- as.dist(dist_ambiental_mat)

p_sil_lin <- fviz_nbclust(as.matrix(dist_linear), pam, method = "silhouette", k.max = 8) +
  labs(title = "Linear", x = "Clusters", y = "Silhouette Width") +
  theme_minimal()

p_sil_out <- fviz_nbclust(as.matrix(dist_outline), pam, method = "silhouette", k.max = 8) +
  labs(title = "Outline", x = "Clusters", y = "Silhouette Width") +
  theme_minimal()

# 1. Armazenar a composição em um objeto
diagnostico_k <- (p_sil_lin | p_sil_out)

# 2. Exibir na tela
print(diagnostico_k)

# 3. Salvar o arquivo
# Formato PNG para apresentações
ggsave("Grafico_Diagnostico_K_Ideal.png", 
       plot = diagnostico_k, 
       width = 20, 
       height = 10, 
       units = "cm", 
       dpi = 300)

# 3. DEFINIÇÃO DE PARÂMETROS FINAIS (Ajuste conforme os gráficos acima)
K_LIN    <- 2
K_OUT    <- 2  
SET_SEED <- 123
PALETA   <- "Set1"

# 4. EXECUÇÃO DA CLUSTERIZAÇÃO (PAM)
set.seed(SET_SEED)
res_pam_linear  <- pam(dist_linear, k = K_LIN)
res_pam_outline <- pam(dist_outline, k = K_OUT)

# Define colors based on the number of clusters (K_LIN)
# 1. Definindo as cores manuais
cores_sil <- c("#6A0572", "#84A59D")

# 5. GRÁFICOS DE PERFIL DE SILHUETA (Opcional: Validação individual)
# --- Perfil Linear ---
sil_plot_lin <- fviz_silhouette(res_pam_linear, label = FALSE, print.summary = FALSE) +
  scale_fill_manual(values = cores_sil) +
  scale_color_manual(values = cores_sil) +
  labs(
    title = "Silhouette Profile: Linear Traits",
    subtitle = paste("Average width:", round(res_pam_linear$silinfo$avg.width, 3)),
    x = "Specimens",
    y = "Silhouette Width"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(), 
    legend.position = "none",              
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, face = "italic"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1) # Nomes em 90 graus
  )

# --- Perfil Outline ---
sil_plot_out <- fviz_silhouette(res_pam_outline, label = FALSE, print.summary = FALSE) +
  scale_fill_manual(values = cores_sil) +
  scale_color_manual(values = cores_sil) +
  labs(
    title = "Silhouette Profile: Geometric Shape",
    subtitle = paste("Average width:", round(res_pam_outline$silinfo$avg.width, 3)),
    x = "Specimens",
    y = "Silhouette Width"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, face = "italic"),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1) # Nomes em 90 graus
  )

print(sil_plot_lin)
print(sil_plot_out)

# ==========================================================
# SALVAMENTO INDIVIDUAL EM PNG (ALTA RESOLUÇÃO)
# ==========================================================

# 1. Salvar o gráfico do Perfil Linear
ggsave(
  filename = "2_Tabelas_e_Imagens/Perfil_Silhueta_Linear.png", 
  plot = sil_plot_lin, 
  width = 25,           # Largura suficiente para os nomes dos espécimes
  height = 16,          # Altura proporcional
  units = "cm", 
  dpi = 600,            # Qualidade de impressão (600 DPI é excelente)
  bg = "white"          # Garante o fundo branco
)

# 2. Salvar o gráfico do Perfil Outline (Morfometria Geométrica)
ggsave(
  filename = "2_Tabelas_e_Imagens/Perfil_Silhueta_Outline.png", 
  plot = sil_plot_out, 
  width = 25, 
  height = 16, 
  units = "cm", 
  dpi = 600, 
  bg = "white"
)

cat("As imagens PNG individuais foram salvas na pasta '2_Tabelas_e_Imagens'!")
# 6. PREPARAÇÃO DA BASE DE DADOS ESPACIAL
map_data <- dados %>%
  mutate(
    Cluster_Linear  = as.factor(res_pam_linear$clustering),
    Cluster_Outline = as.factor(res_pam_outline$clustering)
  ) %>%
  st_as_sf(coords = c("Log", "Lat"), crs = 4326)

library(patchwork)

# --------------------------------
# PERFIS DE SILHUETA LIMPOS
# --------------------------------

sil_plot_lin <- fviz_silhouette(
  res_pam_linear,
  label = FALSE,
  print.summary = FALSE
) +
  scale_fill_manual(values = cores_sil) +
  scale_color_manual(values = cores_sil) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "Specimens",
    y = "Silhouette Width"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    
    # remove títulos
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    
    # ajusta nomes dos espécimes
    axis.text.x = element_text(
      angle = 90,
      size = 5,       # diminui fonte
      vjust = 0.5,
      hjust = 1
    )
  )


sil_plot_out <- fviz_silhouette(
  res_pam_outline,
  label = FALSE,
  print.summary = FALSE
) +
  scale_fill_manual(values = cores_sil) +
  scale_color_manual(values = cores_sil) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "Specimens",
    y = "Silhouette Width"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = "none",
    
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    
    axis.text.x = element_text(
      angle = 90,
      size = 5,
      vjust = 0.5,
      hjust = 1
    )
  )


# --------------------------------
# FIGURA COMPOSTA
# --------------------------------

fig_silhouette <- (
  p_sil_lin | p_sil_out
) /
  (
    sil_plot_lin | sil_plot_out
  ) +
  plot_annotation(
    tag_levels = "A"
  )

print(fig_silhouette)

ggsave(
  "Figura_Composta_Clustering.png",
  plot = fig_silhouette,
  width = 36,   # aumentei largura
  height = 22,
  units = "cm",
  dpi = 600
)

# 0. CORES E MAPEAMENTO (Nomes em Inglês para a Legenda)
# ------------------------------------------------------------------------------
my_colors <- c("#6A0572", "#AB83A1", "#F67280", "#F8B195", "#84A59D")

cores_biomas <- c(
  "Amazon Forest"   = my_colors[5], # Verde Musgo/Sage
  "Caatinga"        = my_colors[2], # Lilás
  "Cerrado"         = my_colors[3], # Coral
  "Atlantic Forest" = my_colors[1], # Violeta
  "Pampa"           = my_colors[4], # Pêssego
  "Pantanal"        = "#D3D3D3"     # Cinza
)

# ==============================================================================
# 7. PREPARAÇÃO DOS DADOS ESPACIAIS (Link Corrigido IBGE - Sem geobr)
# ==============================================================================
library(sf)
library(dplyr)
library(ggplot2)
library(patchwork)
library(vegan) # Para o teste de Mantel

# URL direta e oficial do IBGE para os Biomas
url_biomas <- "https://geoftp.ibge.gov.br/informacoes_ambientais/estudos_ambientais/biomas/vetores/Biomas_250mil.zip"
temp_zip   <- tempfile(fileext = ".zip")
temp_dir   <- tempdir()

message("Baixando dados de biomas direto do IBGE...")
download.file(url_biomas, temp_zip, mode = "wb", quiet = TRUE)
unzip(temp_zip, exdir = temp_dir)

# Carrega o arquivo shapefile descompactado
shp_path   <- list.files(temp_dir, pattern = "\\.shp$", full.names = TRUE)
biomas_raw <- st_read(shp_path, quiet = TRUE)

# Cria e padroniza o mapa de fundo usando a coluna exata: 'Bioma'
biomas_mapa <- biomas_raw %>%
  st_transform(4326) %>%
  filter(Bioma != "Sistema Costeiro") %>%
  mutate(name_biome = case_when(
    Bioma == "Amazônia"        ~ "Amazon Forest",
    Bioma == "Mata Atlântica"  ~ "Atlantic Forest",
    Bioma == "Cerrado"         ~ "Cerrado",
    Bioma == "Caatinga"        ~ "Caatinga",
    Bioma == "Pampa"           ~ "Pampa",
    Bioma == "Pantanal"        ~ "Pantanal",
    TRUE ~ as.character(Bioma)
  ))

# Limpa e traduz a base dos seus espécimes (pontos amostrais)
map_data_clean <- map_data %>%
  filter(`Dominio fitogeografico` != "Sistema Costeiro") %>%
  mutate(`Dominio fitogeografico` = case_when(
    `Dominio fitogeografico` == "Amazônia"        ~ "Amazon Forest",
    `Dominio fitogeografico` == "Mata Atlântica"  ~ "Atlantic Forest",
    `Dominio fitogeografico` == "Cerrado"         ~ "Cerrado",
    `Dominio fitogeografico` == "Caatinga"        ~ "Caatinga",
    `Dominio fitogeografico` == "Pampa"           ~ "Pampa",
    `Dominio fitogeografico` == "Pantanal"        ~ "Pantanal",
    TRUE ~ `Dominio fitogeografico`
  ))

print("--- DADOS ESPACIAIS PREPARADOS COM SUCESSO ---")

# ==============================================================================
# 8. CONFIGURAÇÕES ESTÉTICAS (Padrão de Publicação Internacional)
# ==============================================================================

# Cores suaves para os Biomas de fundo (não ofuscam os pontos dos clusters)
my_colors <- c("#6A0572", "#AB83A1", "#F67280", "#F8B195", "#84A59D")
cores_biomas <- c(
  "Amazon Forest"   = my_colors[5], # Verde Sage
  "Caatinga"        = my_colors[2], # Lilás
  "Cerrado"         = my_colors[3], # Coral
  "Atlantic Forest" = my_colors[1], # Violeta
  "Pampa"           = my_colors[4], # Pêssego
  "Pantanal"        = "#D3D3D3"     # Cinza Claro
)

# Cores marcantes para os Clusters (Tons Premium)
cores_clusters <- c("#DAA520", "#800000", "#4682B4", "#D2691E")

# ==============================================================================
# 9. CONSTRUÇÃO DOS MAPAS (ggplot2)
# ==============================================================================

# Mapa A: Cluster Linear
mapa_linear <- ggplot() +
  geom_sf(data = biomas_mapa, aes(fill = name_biome), color = "white", size = 0.1, alpha = 0.6) +
  geom_sf(data = map_data_clean, aes(color = Cluster_Linear), size = 2.5, alpha = 0.9) +
  scale_fill_manual(values = cores_biomas, name = "Biomes") +
  scale_color_manual(values = cores_clusters, name = "Linear Clusters") +
  theme_minimal(base_size = 11) +
  labs(title = "Linear Morphometrics Clustering") +
  theme(panel.grid = element_blank(), axis.text = element_blank())

# Mapa B: Cluster Outline
mapa_outline <- ggplot() +
  geom_sf(data = biomas_mapa, aes(fill = name_biome), color = "white", size = 0.1, alpha = 0.6) +
  geom_sf(data = map_data_clean, aes(color = Cluster_Outline), size = 2.5, alpha = 0.9) +
  scale_fill_manual(values = cores_biomas, name = "Biomes") +
  scale_color_manual(values = cores_clusters, name = "Outline Clusters") +
  theme_minimal(base_size = 11) +
  labs(title = "Outline Geometric Clustering") +
  theme(panel.grid = element_blank(), axis.text = element_blank())

# ==============================================================================
# 10. ANÁLISE DE MANTEL (Geográfica vs Ambiental)
# ==============================================================================
message("Calculando testes de Mantel...")
mantel_linear  <- mantel(dist_linear, dist_ambiental, method = "pearson", permutations = 9999)
mantel_outline <- mantel(dist_outline, dist_ambiental, method = "pearson", permutations = 9999)

# Armazena os resultados formatados para conferência rápida no console
legenda_mantel_lin <- paste0("Linear - Mantel r: ", round(mantel_linear$statistic, 3), " | p-val: ", mantel_linear$signif)
legenda_mantel_out <- paste0("Outline - Mantel r: ", round(mantel_outline$statistic, 3), " | p-val: ", mantel_outline$signif)

cat("\n--- RESULTADOS DO TESTE DE MANTEL ---\n")
print(legenda_mantel_lin)
print(legenda_mantel_out)
cat("-------------------------------------\n")

# ==============================================================================
# 11. MONTAGEM E SALVAMENTO DA FIGURA COMPOSTA FINAL (Patchwork)
# ==============================================================================

# Monta o painel combinando os Mapas (Top) com os Gráficos de Silhueta da Seção 6 (Bottom)
fig_final <- (mapa_linear | mapa_outline) / (sil_plot_lin | sil_plot_out) + 
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")

plot(fig_final)

# Salva o arquivo final em alta resolução (600 DPI) para publicação
ggsave(
  "Figura_Final_Completa_Biomas.png",
  plot   = fig_final,
  width  = 38,
  height = 25,
  units  = "cm",
  dpi    = 600
)

message("Concluído! A imagem 'Figura_Final_Completa_Biomas.png' foi gerada com sucesso.")
























# 8. FUNÇÃO NATURE STYLE (Correção do Erro de Jitter)
# ------------------------------------------------------------------------------
gerar_mapa_nature <- function(data_pontos, data_biomas, coluna_cluster, titulo, res_mantel, res_pam) {
  
  label_stats <- paste0(
    "Mantel r = ", round(res_mantel$statistic, 3), " (P = ", res_mantel$signif, ")\n",
    "Si = ", round(res_pam$silinfo$avg.width, 3)
  )
  
  # APLICANDO JITTER ESPACIAL NOS DADOS (Evita o erro do geom_sf)
  # Isso desloca os pontos em aproximadamente 0.1 graus de forma aleatória
  data_jittered <- st_jitter(data_pontos, factor = 0.008) 
  
  ggplot() +
    # 1. BIOMAS (Fundo)
    geom_sf(data = data_biomas, 
            aes(fill = name_biome), 
            color = "white", size = 0.05, alpha = 0.25) + 
    scale_fill_manual(values = cores_biomas, name = "Biomes") +
    
    ggnewscale::new_scale_fill() +
    
    # 2. PONTOS (Usando os dados com Jitter aplicado e Transparência)
    geom_sf(data = data_jittered, 
            aes(fill = !!sym(coluna_cluster)), 
            shape = 21,      
            color = "white", 
            size = 4.5,      
            stroke = 0.4,    # Borda um pouco mais fina para transparência
            alpha = 0.6) +   # Transparência para ver sobreposição (Overplotting)
    
    scale_fill_manual(values = cores_clusters, name = "Morph. Clusters") +
    
    # 3. ESTATÍSTICAS
    annotate("text", x = -35, y = -33, label = label_stats, 
             size = 3.2, fontface = "bold.italic", color = "grey20", hjust = 1) +
    
    theme_void() +
    labs(title = titulo) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b=15)),
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      legend.position = "right",
      plot.margin = margin(10, 10, 10, 10)
    )
}

# 9. EXECUÇÃO
p1 <- gerar_mapa_nature(map_data_clean, biomas_mapa, "Cluster_Linear", 
                        "Linear Traits Analysis", mantel_linear, res_pam_linear)

p2 <- gerar_mapa_nature(map_data_clean, biomas_mapa, "Cluster_Outline", 
                        "Geometric Shape Analysis", mantel_outline, res_pam_outline)

mapa_final <- (p1 | p2) + plot_layout(guides = "collect")
print(mapa_final)

# Salvando o Mapa Final
ggsave("Figure_1_Morphological_Variation.png", 
       plot = mapa_final, width = 35, height = 18, units = "cm", dpi = 600)

# Criando a tabela de dados suplementares
supplementary_data <- map_data_clean %>%
  # EXTRAIR coordenadas da geometria antes de descartá-la
  mutate(
    Longitude = st_coordinates(.)[,1],
    Latitude  = st_coordinates(.)[,2]
  ) %>%
  # Agora removemos a coluna de geometria para salvar como tabela comum
  st_drop_geometry() %>%
  # Seleciona e renomeia
  select(
    `Phytogeographical Domain` =  `Dominio fitogeografico`,
    Latitude, 
    Longitude,
    `Cluster Linear` = Cluster_Linear,
    `Cluster Outline` = Cluster_Outline
  )

# Salvando em CSV
write_csv(supplementary_data, "Supplementary_Data_Clustering.csv")

# Criando a Tabela Geral de Resultados Estatísticos
estatisticas_resumo <- data.frame(
  Metodo = c("Linear Morphometry", "Geometric Morphometry (Outline)"),
  `Number of Clusters (K)` = c(K_LIN, K_OUT),
  `Average Silhouette Width` = c(
    round(res_pam_linear$silinfo$avg.width, 3),
    round(res_pam_outline$silinfo$avg.width, 3)
  ),
  `Mantel Statistic (r)` = c(
    round(mantel_linear$statistic, 3),
    round(mantel_outline$statistic, 3)
  ),
  `Mantel p-value` = c(
    mantel_linear$signif,
    mantel_outline$signif
  )
)

# Salvar a tabela principal de estatísticas
write_csv(estatisticas_resumo, "Table_1_Statistical_mapa.csv")

# Gerar Tabela de Contingência (Frequência de Clusters por Bioma)
# Isso mostra como a morfologia se distribui nos domínios
contingencia_linear <- table(map_data_clean$`Dominio fitogeografico`, map_data_clean$Cluster_Linear)
contingencia_outline <- table(map_data_clean$`Dominio fitogeografico`, map_data_clean$Cluster_Outline)

# Salvar tabelas de frequência
write.csv(contingencia_linear, "Table_S2_Frequency_Linear.csv")
write.csv(contingencia_outline, "Table_S3_Frequency_Outline.csv")

cat("\n--- Tabelas estatísticas exportadas com sucesso! ---\n")
