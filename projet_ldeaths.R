# =============================================================================
#  PROJET SÉRIES TEMPORELLES : décès mensuels par maladies pulmonaires (UK)
#  Données : ldeaths (R, package datasets), janvier 1974 – décembre 1979
#  Auteur  : Christian TOSSE
#
#  Plan :
#   0. Les données et leur structure
#   1. Exploration : graphiques, tendance, saisonnalité, ACF / PACF
#   2. Pourquoi un AR "brut" ne suffit pas (expérience instructive)
#   3. Le modèle retenu : régression harmonique + erreurs AR(p)
#   4. Validation : diagnostic des résidus
#   5. Prévision de 1979 : RMSE, intervalles de prévision
#   6. Bonus : modèle de Serfling et seuil épidémique
#
#  Aucun package externe : uniquement R de base (stats, graphics).
# =============================================================================


# -----------------------------------------------------------------------------
# 0. LES DONNÉES ET LEUR STRUCTURE
# -----------------------------------------------------------------------------

dir.create("figures", showWarnings = FALSE)   # dossier où l'on sauvegarde les graphiques

y <- ldeaths                                  # la série est déjà dans R : on la copie dans y
print(y)                                      # tableau années x mois : 6 lignes, 12 colonnes
class(y)                                      # "ts" : objet série temporelle (valeurs + calendrier)
start(y); end(y)                              # début 1974 mois 1, fin 1979 mois 12
frequency(y)                                  # 12 observations par an -> données mensuelles
length(y)                                     # 72 observations au total (6 ans x 12 mois)
summary(y)                                    # min 1300, médiane 1870, moyenne 2057, max 3891
# Interprétation : ~2000 décès par mois en moyenne, avec un écart énorme
# entre le mois le plus bas (1300) et le plus haut (3891, février 1976).

# Découpage apprentissage / test : on cache la dernière année au modèle
y_train <- window(y, end   = c(1978, 12))     # 1974-1978 : 60 mois pour estimer
y_test  <- window(y, start = c(1979, 1))      # 1979 : 12 mois pour juger les prévisions
n_train <- length(y_train)                    # 60
h       <- length(y_test)                     # horizon de prévision : 12 mois


# -----------------------------------------------------------------------------
# 1. EXPLORATION
# -----------------------------------------------------------------------------

# 1.a Chronogramme : la série dans le temps
png("figures/01_chronogramme.png", width = 900, height = 450)  # ouvre un fichier image
plot(y, main = "Décès mensuels par maladies pulmonaires, UK",
     ylab = "Nombre de décès", xlab = "Année", lwd = 2)          # trace la série
abline(v = 1979, lty = 2, col = "grey40")                        # sépare apprentissage et test
dev.off()                                                        # ferme et écrit le fichier
# Lecture : un pic chaque hiver, un creux chaque été, le même motif chaque
# année -> SAISONNALITÉ de période 12. Les pics varient (1976 très haut) :
# ce sont ces écarts à la saison "normale" qui intéressent la surveillance.

# 1.b Graphique par mois : la saisonnalité vue directement
png("figures/02_monthplot.png", width = 900, height = 450)
monthplot(y, main = "Profil par mois (trait horizontal = moyenne du mois)",
          ylab = "Nombre de décès")                              # une sous-série par mois
dev.off()
# Lecture : janvier-février ~2900, juillet-septembre ~1400. Le mois explique
# à lui seul l'essentiel des variations.

# 1.c Décomposition : tendance + saison + reste
dec <- decompose(y_train)                     # moyenne mobile 2x12 pour la tendance, puis
                                              # moyenne par mois des écarts pour la saison
png("figures/03_decomposition.png", width = 900, height = 700)
plot(dec)                                     # 4 panneaux : observé, tendance, saison, aléatoire
dev.off()
round(dec$figure, 0)                          # les 12 coefficients saisonniers (écart au niveau moyen)
# Lecture : tendance légèrement décroissante ; saison d'amplitude ~ +/-700 ;
# le reste est irrégulier (gros résidu en février 1976).

# 1.d Autocorrélations de la série brute
png("figures/04_acf_pacf_brute.png", width = 900, height = 400)
par(mfrow = c(1, 2))                          # deux graphiques côte à côte
acf(y_train,  lag.max = 36, main = "ACF série brute")    # corrélation entre y_t et y_{t-k}
pacf(y_train, lag.max = 36, main = "PACF série brute")   # corrélation "nette" à k
dev.off()
round(acf(y_train, plot = FALSE, lag.max = 12)$acf[-1], 2)
# ATTENTION à l'axe "Lag" des graphiques : il est en ANNÉES pour un objet ts
# mensuel (1.0 = 12 mois, 0.5 = 6 mois). Les pointillés bleus = +/- 2/sqrt(n).
# Lecture : l'ACF a la forme d'une SINUSOÏDE de période 12 (+0.74 à 1 mois,
# -0.69 à 6 mois, +0.68 à 12 mois) qui décroît très lentement.
# Conclusion : la moyenne de la série dépend du mois -> la série N'EST PAS
# STATIONNAIRE. Or un modèle AR suppose la stationnarité. Il faut d'abord
# traiter la saisonnalité.


# -----------------------------------------------------------------------------
# 2. EXPÉRIENCE : un AR directement sur la série brute
# -----------------------------------------------------------------------------

# On fait "comme si" la série était stationnaire, pour voir ce qui se passe.
ar_brut <- ar(y_train, aic = TRUE, order.max = 15)   # choisit l'ordre p par AIC (Yule-Walker)
ar_brut$order                                        # p = 5
prev_ar_brut <- predict(ar_brut, n.ahead = h)$pred   # prévisions 1979

# Prévision de référence très simple : "naïf saisonnier"
# (janvier 1979 = janvier 1978, etc.). Tout modèle sérieux doit faire mieux.
prev_naif <- ts(as.numeric(window(y_train, start = c(1978, 1))),
                start = c(1979, 1), frequency = 12)

rmse <- function(obs, prev) sqrt(mean((obs - prev)^2))  # erreur quadratique moyenne, en décès
rmse(y_test, prev_ar_brut)                   # ~317
rmse(y_test, prev_naif)                      # ~266
# Leçon : l'AR brut fait PIRE que la méthode naïve. Un AR(5) ne "voit" que
# 5 mois en arrière : il ne peut pas reproduire un cycle de 12 mois, et ses
# prévisions reviennent vite vers la moyenne globale (~2000).


# -----------------------------------------------------------------------------
# 3. LE MODÈLE RETENU : régression harmonique + erreurs AR(p)
# -----------------------------------------------------------------------------
# Idée : y_t = [tendance + saison déterministe] + u_t,   u_t ~ AR(p)
#   - tendance  : a + b*t
#   - saison    : somme de K paires (cosinus, sinus) de période 12, 12/2, ...
#   - u_t       : écart à la saison normale, stationnaire, modélisé par un AR

# Fonction qui fabrique les régresseurs harmoniques pour des temps tt
harmoniques <- function(tt, K) {
  M <- NULL                                            # matrice vide au départ
  for (k in 1:K) {                                     # pour chaque harmonique k
    M <- cbind(M, cos(2 * pi * k * tt / 12),           # cosinus de période 12/k
                  sin(2 * pi * k * tt / 12))           # sinus de période 12/k
  }
  colnames(M) <- paste0(rep(c("cos", "sin"), K), rep(1:K, each = 2))  # noms lisibles
  M                                                    # renvoie la matrice
}

t_train <- 1:n_train                         # temps 1..60 pour l'apprentissage
t_test  <- (n_train + 1):(n_train + h)       # temps 61..72 pour 1979

# 3.a Choix de K (nombre d'harmoniques) par régression ordinaire
for (K in 1:3) {
  fit_lm <- lm(y_train ~ t_train + harmoniques(t_train, K))     # MCO
  cat("K =", K, " AIC =", round(AIC(fit_lm), 1),
      " R2 ajusté =", round(summary(fit_lm)$adj.r.squared, 3), "\n")
}
K <- 2                                        # K = 2 : AIC minimal ; K = 3 n'apporte rien
reg_saison <- lm(y_train ~ t_train + harmoniques(t_train, K))
summary(reg_saison)
# Lecture : tendance -4.9 décès/mois (significative, p = 0.008) ;
# cos1, sin1, sin2 très significatifs ; R2 ajusté = 0.85 :
# tendance + saison expliquent 85 % de la variance.

# 3.b Les résidus de cette régression : sont-ils indépendants ?
u <- resid(reg_saison)                       # u_t = écart à la saison normale
png("figures/05_acf_pacf_residus_reg.png", width = 900, height = 400)
par(mfrow = c(1, 2))
acf(u,  lag.max = 24, main = "ACF de u (après retrait saison)")
pacf(u, lag.max = 24, main = "PACF de u")
dev.off()
round(pacf(u, plot = FALSE, lag.max = 6)$acf, 2)
2 / sqrt(n_train)                            # bande de significativité approx. : 0.26
# Lecture : la sinusoïde a disparu (saison bien retirée). Il reste une
# PACF à -0.26 au retard 2, à la limite de la bande. Règle : pour un AR(p),
# la PACF s'annule après le retard p -> candidat p = 2.

# 3.c Estimation conjointe (maximum de vraisemblance) pour p = 0..3
X_train <- cbind(t = t_train, harmoniques(t_train, K))    # régresseurs apprentissage
X_test  <- cbind(t = t_test,  harmoniques(t_test,  K))    # régresseurs 1979 (connus d'avance)
for (p in 0:3) {
  f <- arima(y_train, order = c(p, 0, 0), xreg = X_train, method = "ML")
  cat("p =", p, " AIC =", round(AIC(f), 1), " BIC =", round(BIC(f), 1), "\n")
}
# AIC minimal pour p = 2 (832.0) ; le BIC, plus sévère, hésite avec p = 0.
# On retient p = 2 : AIC + PACF + résidus plus propres (voir section 4).

modele <- arima(y_train, order = c(2, 0, 0), xreg = X_train, method = "ML")
modele                                       # coefficients et écarts-types
# ar1 = 0.14 (e.t. 0.12, non significatif) ; ar2 = -0.31 (e.t. 0.13, |z| = 2.3,
# significatif à 5 %) : un mois anormalement haut tend à être suivi, deux
# mois plus tard, d'un mois un peu plus bas (après un pic, une partie des
# personnes fragiles est déjà décédée : effet de "moisson").
# Tendance : -5.1 décès/mois ; sigma^2 = 45455 -> écart-type d'innovation ~213.

# Stationnarité de la partie AR : racines de 1 - phi1 z - phi2 z^2 hors du disque unité
phi <- coef(modele)[c("ar1", "ar2")]
Mod(polyroot(c(1, -phi)))                    # modules = 1.79 > 1 -> AR stationnaire, OK


# -----------------------------------------------------------------------------
# 4. VALIDATION : diagnostic des résidus
# -----------------------------------------------------------------------------
e <- resid(modele)                           # résidus = innovations estimées

png("figures/06_diagnostic.png", width = 900, height = 700)
par(mfrow = c(2, 2))
plot(e, main = "Résidus", ylab = ""); abline(h = 0, lty = 2)     # centrés ? motif ?
acf(e, lag.max = 24, main = "ACF des résidus")                    # tout dans la bande ?
hist(e, breaks = 15, main = "Histogramme", xlab = "")             # forme
qqnorm(e); qqline(e)                                              # normalité
dev.off()

Box.test(e, lag = 24, type = "Ljung-Box", fitdf = 2)
# H0 : pas d'autocorrélation jusqu'au retard 24. p ~ 0.70 -> on ne rejette
# pas : les résidus sont un bruit blanc, le modèle a capté la dépendance.
shapiro.test(e)
# H0 : normalité. p ~ 1e-05 -> normalité REJETÉE. Cause : quelques mois
# extrêmes (février 1976 surtout) alourdissent la queue droite.
# Conséquence : les estimations restent valables, mais les intervalles de
# prévision (qui supposent la normalité) sont à lire avec prudence.
# C'est une limite à mentionner honnêtement.


# -----------------------------------------------------------------------------
# 5. PRÉVISION DE 1979 ET ÉVALUATION
# -----------------------------------------------------------------------------
pr <- predict(modele, n.ahead = h, newxreg = X_test)  # prévisions + écarts-types
prev  <- pr$pred                                      # prévision ponctuelle
borne_inf <- prev - 1.96 * pr$se                      # intervalle de prévision à 95 %
borne_sup <- prev + 1.96 * pr$se

rmse(y_test, prev)                                    # ~161 décès
# Comparaison : naïf saisonnier 266, AR brut 317 -> erreur réduite de ~40 %.
mean(y_test >= borne_inf & y_test <= borne_sup)       # = 1 : les 12 mois sont dans l'IP 95 %
# (décembre 1979, 1915, frôle la borne basse 1906 : hiver doux).

resultat <- data.frame(mois = month.abb,              # tableau récapitulatif
                       observe = as.numeric(y_test),
                       prevu = round(as.numeric(prev)),
                       inf95 = round(as.numeric(borne_inf)),
                       sup95 = round(as.numeric(borne_sup)))
print(resultat)

png("figures/07_prevision_1979.png", width = 900, height = 450)
plot(window(y, start = c(1977, 1)), ylim = c(800, 3600), lwd = 2,
     main = "Prévision 1979 : harmoniques + AR(2)", ylab = "Décès", xlab = "")
polygon(c(time(prev), rev(time(prev))),
        c(borne_inf, rev(borne_sup)), col = rgb(0, 0, 1, 0.15), border = NA)  # bande 95 %
lines(prev, col = "blue", lwd = 2, lty = 2)           # prévision
lines(prev_naif, col = "darkorange", lty = 3, lwd = 2) # référence naïve
legend("topright", c("Observé", "Prévision AR", "Naïf saisonnier", "IP 95 %"),
       col = c("black", "blue", "darkorange", rgb(0, 0, 1, 0.3)),
       lty = c(1, 2, 3, 1), lwd = c(2, 2, 2, 8), bty = "n")
dev.off()


# -----------------------------------------------------------------------------
# 6. BONUS : MODÈLE DE SERFLING ET SEUIL ÉPIDÉMIQUE
# -----------------------------------------------------------------------------
# Serfling (1963) : niveau de base = tendance + une harmonique annuelle,
# estimé par MCO SANS les périodes épidémiques ; seuil = borne supérieure
# de l'intervalle de prévision à 95 %. Un mois au-dessus = excès de décès.

t_all <- 1:length(y)                                  # temps 1..72 (toute la série)
base  <- data.frame(y = as.numeric(y), t = t_all,
                    cos1 = cos(2 * pi * t_all / 12),
                    sin1 = sin(2 * pi * t_all / 12))

# Étape 1 : ajustement sur tout, pour repérer les mois clairement épidémiques
s1 <- lm(y ~ t + cos1 + sin1, data = base)
seuil1 <- predict(s1, newdata = base, interval = "prediction", level = 0.95)[, "upr"]
garder <- base$y <= seuil1                            # on exclut les mois au-dessus du seuil

sum(!garder)                                          # 1 mois exclu : février 1976
# Étape 2 : on ré-estime le niveau de base sans eux (base "non épidémique")
serfling <- lm(y ~ t + cos1 + sin1, data = base[garder, ])
summary(serfling)
ip <- predict(serfling, newdata = base, interval = "prediction", level = 0.95)
base$attendu <- ip[, "fit"]                           # décès attendus hors épidémie
base$seuil   <- ip[, "upr"]                           # seuil épidémique
base$alerte  <- base$y > base$seuil                   # TRUE = excès
base$exces   <- ifelse(base$alerte, base$y - base$attendu, 0)  # décès en excès

base$date <- format(as.Date(paste(floor(time(y)), cycle(y), 1, sep = "-")), "%Y-%m")
print(base[base$alerte, c("date", "y", "attendu", "seuil", "exces")])
# Résultat : 3 alertes. Février 1976 (+1074 décès en excès), mars 1976 (+500),
# janvier 1979 (+504). Février-mars 1976 : hiver 1975-76, probablement une
# épidémie de grippe (hypothèse à vérifier dans la littérature avant de l écrire).
sum(base$exces)                                       # total des décès en excès (~2078)

png("figures/08_serfling.png", width = 900, height = 450)
plot(y, lwd = 2, ylim = c(1000, 4000), main = "Modèle de Serfling : attendu et seuil",
     ylab = "Décès", xlab = "")
lines(ts(base$attendu, start = c(1974, 1), frequency = 12), col = "darkgreen", lwd = 2)
lines(ts(base$seuil,   start = c(1974, 1), frequency = 12), col = "red", lty = 2, lwd = 2)
points(time(y)[base$alerte], base$y[base$alerte], pch = 19, col = "red", cex = 1.4)
legend("topright", c("Observé", "Attendu (base)", "Seuil 95 %", "Alerte"),
       col = c("black", "darkgreen", "red", "red"), lty = c(1, 1, 2, NA),
       pch = c(NA, NA, NA, 19), lwd = 2, bty = "n")
dev.off()
# =============================================================================
# FIN
# =============================================================================
