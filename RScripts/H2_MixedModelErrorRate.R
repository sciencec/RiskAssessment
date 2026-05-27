# ============================================================
# Hypothese 2: Anwendungsgenauigkeit (Error Rate)
# Generalisiertes Mixed-Effects-Modell (GLMM): Error Rate
# Zufallseffekte: Versuchsperson (SUBJECT) und
#                 Stimuluskombination 
#                 (SINGLE_OCC x HIGH_SPEED x INTELL_MATCH)
# Linkfunktion: logistisch (da ERROR binär: 0/1)
# ============================================================

library(lme4)
library(lmerTest)
library(ggplot2)

# ------------------------------------------------------------
# 1. Daten einlesen
# ------------------------------------------------------------
ErrorDaten <- learning_timings_errors_neu

# ------------------------------------------------------------
# 2. Variablen vorbereiten
# ------------------------------------------------------------

# Leere Zeilen entfernen (letzte Zeile ist komplett NA)
ErrorDaten <- ErrorDaten[!is.na(ErrorDaten$SUBJECT), ]

# CONDITION als Faktor mit sprechenden Labels (FFT = Referenzkategorie)
ErrorDaten$CONDITION <- factor(
  ErrorDaten$CONDITION,
  levels = c(1, 2, 3),
  labels = c("FFT", "Complete Tree", "Risk Table")
)

# Stimuluskombination als Item-Faktor
ErrorDaten$ITEM <- factor(
  paste(ErrorDaten$SINGLE_OCC,
        ErrorDaten$HIGH_SPEED,
        ErrorDaten$INTELL_MATCH,
        sep = "_")
)

# SUBJECT als Faktor
ErrorDaten$SUBJECT <- factor(ErrorDaten$SUBJECT)

# ERROR sicherstellen als Integer (0/1)
ErrorDaten$ERROR <- as.integer(ErrorDaten$ERROR)

cat("Datenbeschreibung:\n")
cat("  Gesamtbeobachtungen:          ", nrow(ErrorDaten), "\n")
cat("  Versuchspersonen:             ", nlevels(ErrorDaten$SUBJECT), "\n")
cat("  Items (Stimuluskombinationen):", nlevels(ErrorDaten$ITEM), "\n")
cat("  Fehlerrate gesamt:            ", round(mean(ErrorDaten$ERROR) * 100, 1), "%\n")
cat("\n  Fehlerrate je Bedingung:\n")
print(tapply(ErrorDaten$ERROR, ErrorDaten$CONDITION, function(x) paste0(round(mean(x) * 100, 1), "%")))

# ------------------------------------------------------------
# 3. Modell 1: GLMM mit allen Beobachtungen
#    Random intercepts für SUBJECT und ITEM
#    Logistische Linkfunktion (family = binomial)
# ------------------------------------------------------------
cat("\n--- Modell 1: Alle Beobachtungen ---\n")

m1 <- glmer(
  ERROR ~ CONDITION +
    (1 | SUBJECT) +
    (1 | ITEM),
  data = ErrorDaten,
  family = binomial(link = "logit"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

print(summary(m1))

# Odds Ratios mit 95%-KI
cat("\nOdds Ratios (exp(Koeffizient)) mit 95%-KI:\n")
or_ci <- exp(cbind(OR = fixef(m1), confint(m1, method = "Wald", parm = "beta_")))
print(round(or_ci, 3))

# ------------------------------------------------------------
# 4. Varianzdekomposition
# ------------------------------------------------------------
cat("\n--- Varianzdekomposition (Modell 1) ---\n")
vc <- as.data.frame(VarCorr(m1))
vc$pct <- round(vc$vcov / sum(vc$vcov) * 100, 1)
print(vc[, c("grp", "vcov", "pct")])

# ------------------------------------------------------------
# 5. Likelihood-Ratio-Test: Effekt von CONDITION
# ------------------------------------------------------------
cat("\n--- Likelihood-Ratio-Test: Effekt von CONDITION ---\n")

m1_null <- glmer(
  ERROR ~ 1 +
    (1 | SUBJECT) +
    (1 | ITEM),
  data = ErrorDaten,
  family = binomial(link = "logit"),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)

print(anova(m1_null, m1))

# ------------------------------------------------------------
# 6. Visualisierung: Fehlerrate je Bedingung
# ------------------------------------------------------------
fehlerrate_summary <- aggregate(ERROR ~ CONDITION, data = ErrorDaten,
                                 FUN = function(x) c(mean = mean(x), se = sd(x)/sqrt(length(x))))
fehlerrate_df <- data.frame(
  CONDITION = fehlerrate_summary$CONDITION,
  mean = fehlerrate_summary$ERROR[, "mean"],
  se   = fehlerrate_summary$ERROR[, "se"]
)

p <- ggplot(fehlerrate_df, aes(x = CONDITION, y = mean * 100, fill = CONDITION)) +
  geom_bar(stat = "identity", alpha = 0.8, width = 0.5) +
  geom_errorbar(aes(ymin = (mean - se) * 100, ymax = (mean + se) * 100),
                width = 0.15, linewidth = 0.8) +
  scale_fill_manual(values = c("FFT" = "#2196F3",
                               "Complete Tree" = "#F44336",
                               "Risk Table" = "#4CAF50")) +
  labs(
    title = "Error rate by Condition",
    subtitle = "Mean ± Standard error (trial level)",
    x = "Condition",
    y = "Error rate (%)",
    fill = "Condition"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

ggsave("Plot_H2_MixedModelErrorRate_EN.png", p, width = 7, height = 5, dpi = 150)
cat("\nAbbildung gespeichert: Plot_H2_MixedModelErrorRate_EN.png\n")

# ------------------------------------------------------------
# 7. Hinweise zur Interpretation
# ------------------------------------------------------------
cat("
HINWEISE ZUR INTERPRETATION:
- Da ERROR binär ist (0/1), wird ein GLMM mit logistischer
  Linkfunktion verwendet (kein lineares Modell).
- Die Koeffizienten (b) liegen auf der Log-Odds-Skala.
  exp(b) = Odds Ratio: Werte > 1 bedeuten höhere Fehlerwahrschein-
  lichkeit im Vergleich zur Referenzkategorie FFT.
- Ein Odds Ratio von z.B. 1.5 bedeutet, dass die Chance einen
  Fehler zu machen in dieser Bedingung 1.5-mal so hoch ist
  wie in der FFT-Bedingung.
- Der Likelihood-Ratio-Test vergleicht das Modell mit CONDITION
  gegen ein Nullmodell ohne Bedingungseffekt.
- Die Varianzdekomposition zeigt, wie viel Varianz auf
  Personen- vs. Item-Unterschiede entfällt.
")
