<div align="center">

# 🌿 Overview of the Analyses

This document provides a general overview of the analyses performed by each script in the project.

**Juliana Aljahara · Niksoney A. Mendonça · Thaís E. Almeida**

</div>

---

## 📏 `Linear_.R` — Linear morphometrics

This script analyzes measurements describing the size and proportions of leaf structures.

The data are initially organized, missing values are estimated, and repeated measurements are summarized by specimen. Extreme values are then identified to prevent highly discrepant observations from disproportionately affecting the comparisons.

The script compares leaf measurements among phytogeographic domains, identifies differences between groups, and summarizes the main patterns of morphological variation. It also assesses how accurately linear measurements distinguish the domains and calculates morphological distances among specimens.

Each leaf measurement is subsequently related to environmental variables. Finally, allometric analyses examine how different structures grow relative to one another and determine whether these relationships are proportional or allometric.

**Main outputs:** comparisons among domains, morphological ordination, specimen classification, a distance matrix, morphology–environment relationships, and allometric patterns.

---

## 🍂 `Outline_.R` — Leaf outline morphometrics

This script analyzes leaf shape from the coordinates of leaf outlines.

The outlines are imported, organized by specimen and phytogeographic domain, and converted into numerical shape descriptors. The number of harmonics required to represent leaf shape adequately is also evaluated.

The main trends in leaf-shape variation are summarized and compared among domains. The script also evaluates whether leaf outlines can correctly classify specimens and calculates representative mean shapes for each group.

Global differences among domains are tested, a shape-based distance matrix is produced, and spatial representations of the analyzed specimens are generated.

**Main outputs:** outline representations, harmonic diagnostics, shape ordination, specimen classification, mean shapes, differences among domains, and a distance matrix.

---

## 🌎 `Ambiental_.R` — Environmental and spatial characterization

This script characterizes the environmental conditions of specimen occurrence sites.

Climatic, edaphic, topographic, and vegetation-structure information is extracted from environmental layers for each collection coordinate. Duplicate records are removed, variables with excessive missing data are excluded, and the remaining missing values are estimated using information from nearby areas.

The environmental dataset is then reduced to avoid retaining strongly redundant variables. Environmental differences among phytogeographic domains are summarized, visualized, and statistically tested.

The script also identifies morphological clusters based on linear measurements and leaf outlines, evaluates the consistency of these clusters, and displays their geographic distribution. Finally, it tests whether morphological differences among specimens are associated with environmental differences.

**Main outputs:** an environmental dataset for the specimens, selected environmental variables, environmental ordination, differences among domains, morphological clusters, maps, and morphology–environment relationships.

---

## 🔗 `Integrativa_.R` — Integration of morphology and environment

This script integrates the results obtained from the linear, outline, and environmental analyses.

The datasets are first synchronized so that only specimens shared by all analyses are retained. The script then evaluates how much variation in linear measurements and leaf outlines can be explained by environmental conditions.

The environmental variables contributing most strongly to morphological variation are identified and represented graphically. Patterns derived from linear measurements and leaf outlines are also compared to determine the level of agreement between the two morphometric approaches.

**Main outputs:** environmental effects on linear morphology, environmental effects on leaf shape, identification of the most relevant environmental variables, and agreement between the two morphometric approaches.

---

> 🎯 **Summary:** together, the four scripts help explain how leaf size and shape vary across phytogeographic domains, which environmental conditions are associated with this variation, and whether different morphometric approaches reveal concordant patterns.

