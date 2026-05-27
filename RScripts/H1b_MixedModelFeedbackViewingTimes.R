# ============================================================
# Hypothese 1b: Feedback Viewing Times
# Mixed-Effects-Modell: Feedback Viewing Time
# Zufallseffekte: Versuchsperson (SUBJECT) und
#                 Stimuluskombination 
#                 (SINGLE_OCC x HIGH_SPEED x INTELL_MATCH)
# ============================================================

library(lme4)
library(lmerTest)   # liefert p-Werte für lmer via Satterthwaite
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# 1. Daten einlesen
# ------------------------------------------------------------

FeedbackDaten<-learning_timings_feedback_neu

# ------------------------------------------------------------
# 2. Variablen vorbereiten
# ------------------------------------------------------------

# CONDITION als Faktor mit sprechenden Labels
FeedbackDaten$CONDITION <- factor(
  FeedbackDaten$CONDITION,
  levels = c(1, 2, 3),
  labels = c("FFT", "Complete Tree", "Risk Table")
)

# Stimuluskombination als Item-Faktor
FeedbackDaten$ITEM <- factor(
  paste(FeedbackDaten$SINGLE_OCC,
        FeedbackDaten$HIGH_SPEED,
        FeedbackDaten$INTELL_MATCH,
        sep = "_")
)

# SUBJECT als Faktor
FeedbackDaten$SUBJECT <- factor(FeedbackDaten$SUBJECT)

# FEEDBACK_TIME liegt in Millisekunden vor und ist stark rechtsschief.
# Wir verwenden log-transformierte Werte als abhängige Variable.
# Achtung: Werte von 15000ms sind zensiert (Teilnehmer haben die
# gesamte Zeit gewartet). Diese sollten im Idealfall mit einem
# Zensierungsmodell (Tobit) behandelt werden. Hier werden sie
# zunächst eingeschlossen und die Sensitivität in einem
# zweiten Modell (ohne zensierte Werte) geprüft.
FeedbackDaten$LOG_FEEDBACK_TIME <- log(FeedbackDaten$FEEDBACK_TIME)

# Zensierungs-Indikator (>= 14900ms als konservative Grenze)
FeedbackDaten$CENSORED <- FeedbackDaten$FEEDBACK_TIME >= 14900

cat("Anteil zensierter Beobachtungen:", 
    round(mean(FeedbackDaten$CENSORED) * 100, 1), "%\n")

# ------------------------------------------------------------
# 3. Modell 1: Vollständiges Modell (alle Beobachtungen)
#    Random intercepts für SUBJECT und ITEM
# ------------------------------------------------------------
cat("\n--- Modell 1: Alle Beobachtungen ---\n")

m1 <- lmer(
  LOG_FEEDBACK_TIME ~ CONDITION +
    (1 | SUBJECT) +
    (1 | ITEM),
  data = FeedbackDaten,
  REML = TRUE
)

print(summary(m1))

# Post-hoc-Kontraste: FFT vs. Complete Tree und FFT vs. Risk Table
# (FFT ist die Referenzkategorie durch die Faktor-Kodierung oben)
cat("\nKoeffizienten (Referenz: FFT):\n")
print(coef(summary(m1)))

# ------------------------------------------------------------
# 4. Modell 2: Sensitivitätsanalyse ohne zensierte Beobachtungen
# ------------------------------------------------------------
cat("\n--- Modell 2: Ohne zensierte Beobachtungen (>= 14900ms) ---\n")

FeedbackDaten_unzensiert <- FeedbackDaten[FeedbackDaten$CENSORED == FALSE, ]

cat("Verbleibende Beobachtungen:", nrow(FeedbackDaten_unzensiert), "\n")

m2 <- lmer(
  LOG_FEEDBACK_TIME ~ CONDITION +
    (1 | SUBJECT) +
    (1 | ITEM),
  data = FeedbackDaten_unzensiert,
  REML = TRUE
)

print(summary(m2))

# ------------------------------------------------------------
# 5. Varianzdekomposition: Wie viel Varianz entfällt auf
#    Personen vs. Items vs. Residual?
# ------------------------------------------------------------
cat("\n--- Varianzdekomposition (Modell 1) ---\n")
vc <- as.data.frame(VarCorr(m1))
vc$pct <- round(vc$vcov / sum(vc$vcov) * 100, 1)
print(vc[, c("grp", "vcov", "pct")])

cat("\n--- Varianzdekomposition (Modell 2) ---\n")
vc2 <- as.data.frame(VarCorr(m2))
vc2$pct <- round(vc2$vcov / sum(vc2$vcov) * 100, 1)
print(vc2[, c("grp", "vcov", "pct")])

# ------------------------------------------------------------
# 6. Modellvergleich: Mit vs. ohne CONDITION (Likelihood-Ratio-Test)
# ------------------------------------------------------------
cat("\n--- Likelihood-Ratio-Test: Effekt von CONDITION ---\n")

m1_null <- lmer(
  LOG_FEEDBACK_TIME ~ 1 +
    (1 | SUBJECT) +
    (1 | ITEM),
  data = FeedbackDaten,
  REML = FALSE   # ML für Modellvergleich
)

m1_ml <- lmer(
  LOG_FEEDBACK_TIME ~ CONDITION +
    (1 | SUBJECT) +
    (1 | ITEM),
  data = FeedbackDaten,
  REML = FALSE
)

print(anova(m1_null, m1_ml))

# ------------------------------------------------------------
# 7. Visualisierung: Verteilung der FEEDBACK_TIME je Bedingung
# ------------------------------------------------------------
p <- ggplot(FeedbackDaten, aes(x = CONDITION, y = FEEDBACK_TIME / 1000,
                                fill = CONDITION)) +
  geom_violin(alpha = 0.4, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7) +
  scale_fill_manual(values = c("FFT" = "#2196F3",
                               "Complete Tree" = "#F44336",
                               "Risk Table" = "#4CAF50")) +
  labs(
    title = "Feedback Viewing Time by Condition",
    subtitle = "Raw data (trial level); dashed line = 15-second upper limit",
    x = "Condition",
    y = "Feedback Viewing Time (seconds)",
    fill = "Bedingung"
  ) +
  geom_hline(yintercept = 15, linetype = "dashed", color = "grey40") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

ggsave("Plot_H1b_FBTbyCondition_EN.png", p, width = 8, height = 5, dpi = 150)
cat("\nAbbildung gespeichert: Plot_H1b_FBTbyCondition_EN.png\n")

# ------------------------------------------------------------
# 8. Hinweis zur Interpretation
# ------------------------------------------------------------
cat("
HINWEISE ZUR INTERPRETATION:
- Die Koeffizienten beziehen sich auf log(FEEDBACK_TIME).
  Zur Rücktransformation: exp(Koeffizient) gibt den
  multiplikativen Effekt auf die ursprüngliche Skala an.
- Modell 1 enthält zensierte Beobachtungen (15s-Deckeneffekt).
  Modell 2 ist die Sensitivitätsanalyse ohne diese Werte.
  Wenn beide Modelle zu denselben Schlussfolgerungen führen,
  stärkt das die Robustheit der Ergebnisse.
- Die Varianzdekomposition zeigt, wie viel Varianz auf
  Personen- vs. Item-Unterschiede vs. Residualvarianz entfällt.
  Ein hoher Item-Anteil würde die Notwendigkeit des
  Mixed-Effects-Ansatzes gegenüber Reviewer 3 untermauern.
")
