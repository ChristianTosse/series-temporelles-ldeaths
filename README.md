# Décès mensuels par maladies pulmonaires au Royaume-Uni (1974-1979)

Projet personnel de **séries temporelles** sous R, sur la série `ldeaths` (fournie avec R).

**Question** : comment distinguer une hausse anormale des décès de la saisonnalité attendue ?

## Démarche
- Exploration : chronogramme, saisonnalité, stationnarité, ACF/PACF, décomposition additive (choix justifié).
- Référence de prévision : naïf saisonnier (RMSE 266).
- Option 1, rejetée : AR sur la série brute (RMSE 317).
- Option 2, testée : SARIMA (meilleure RMSE 162 ; le modèle de meilleur AIC prévoit mal).
- **Option 3, retenue : régression harmonique (tendance + 2 harmoniques) à erreurs AR(2)**, avec une RMSE de 161 sur 1979, soit environ 40 % d'erreur en moins que la référence saisonnière ; résidus en bruit blanc (Ljung-Box, p = 0,70).
- **Modèle de Serfling** : niveau attendu et seuil épidémique à 95 %. Trois alertes, dont l'hiver 1975-1976.

## Fichiers
| Fichier | Contenu |
|---|---|
| `projet_ldeaths_complet.html` / `.pdf` | rapport complet : notions, code commenté, sorties, interprétations |
| `projet_ldeaths_complet.Rmd` | source R Markdown (à compiler avec Knit dans RStudio) |
| `projet_ldeaths.R` | script R seul, commenté ligne par ligne |
| `figures/` | figures produites par le script |
| `header.tex`, `boxes.lua`, `style.css` | fichiers de mise en forme nécessaires à la compilation |

## Reproduire
R ≥ 4.0, packages `rmarkdown` et `knitr`. Le modèle n'utilise que R de base (`stats`, `graphics`).

Auteur : Christian TOSSE — M2 Mathématiques des Données, Université Sorbonne Paris Nord.
