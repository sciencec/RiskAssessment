# ============================================================
# Hypothese 1b: Feedback Viewing Times (kumuliert)
# Analyse kumulativer Feedback Viewing Times
# Abhängige Variable: Summe der Feedback-Betrachtungszeit
#                     pro Versuchsperson (in Sekunden)
# Gruppenvergleich: FFT vs. Complete Tree vs. Risk Table
# Methode: Kruskal-Wallis-Test + post-hoc Mann-Whitney-U-Tests
# ============================================================

library(ggplot2)

# ------------------------------------------------------------
# 1. Daten einlesen
# ------------------------------------------------------------
FeedbackDaten <- learning_timings_feedback_neu

# ------------------------------------------------------------
# 2. Variablen vorbereiten
# ------------------------------------------------------------

# CONDITION als Faktor mit sprechenden Labels (FFT = Referenzkategorie)
FeedbackDaten$CONDITION <- factor(
  FeedbackDaten$CONDITION,
  levels = c(1, 2, 3),
  labels = c("FFT", "Complete Tree", "Risk Table")
)

# SUBJECT als Faktor
FeedbackDaten$SUBJECT <- factor(FeedbackDaten$SUBJECT)

# FEEDBACK_TIME in Sekunden umrechnen (liegt in ms vor)
# Umrechnung ms -> s erfolgt direkt bei der Aggregation (siehe unten)

# ------------------------------------------------------------
# 3. Kumulative Betrachtungszeit pro Versuchsperson berechnen
# ------------------------------------------------------------
# Summe aller Feedback-Viewing-Times über alle Trials
# bis zum Erreichen des Lernkriteriums
KumulDaten <- aggregate(
  FEEDBACK_TIME ~ SUBJECT + CONDITION,
  data = FeedbackDaten,
  FUN = sum
)
colnames(KumulDaten)[3] <- "CUM_FEEDBACK_TIME_S"
KumulDaten$CUM_FEEDBACK_TIME_S <- KumulDaten$CUM_FEEDBACK_TIME_S / 1000

cat("Datenbeschreibung (kumulative Ebene):\n")
cat("  Versuchspersonen gesamt:", nrow(KumulDaten), "\n")
cat("\n  N je Bedingung:\n")
print(table(KumulDaten$CONDITION))

# ------------------------------------------------------------
# 4. Deskriptive Statistik
# ------------------------------------------------------------
cat("\n--- Deskriptive Statistik (kumulative Feedback Viewing Time in Sekunden) ---\n")

deskriptiv <- aggregate(
  CUM_FEEDBACK_TIME_S ~ CONDITION,
  data = KumulDaten,
  FUN = function(x) c(
    N      = length(x),
    Mean   = round(mean(x), 1),
    Median = round(median(x), 1),
    SD     = round(sd(x), 1),
    Min    = round(min(x), 1),
    Max    = round(max(x), 1)
  )
)
print(do.call(data.frame, deskriptiv))

# Gesamtzeiten (Minuten) je Bedingung – praktisch interpretierbar
cat("\nGesamte Betrachtungszeit über alle Personen je Bedingung (Minuten):\n")
gesamt <- aggregate(CUM_FEEDBACK_TIME_S ~ CONDITION, data = KumulDaten, FUN = sum)
gesamt$Minuten <- round(gesamt$CUM_FEEDBACK_TIME_S / 60, 1)
print(gesamt[, c("CONDITION", "Minuten")])

# ------------------------------------------------------------
# 5. Voraussetzungsprüfung: Normalverteilung
# ------------------------------------------------------------
cat("\n--- Shapiro-Wilk-Test auf Normalverteilung je Bedingung ---\n")
for (bed in levels(KumulDaten$CONDITION)) {
  x <- KumulDaten$CUM_FEEDBACK_TIME_S[KumulDaten$CONDITION == bed]
  sw <- shapiro.test(x)
  cat(sprintf("  %s: W = %.4f, p = %.4f\n", bed, sw$statistic, sw$p.value))
}
cat("  (Bei p < .05: Normalverteilung abgelehnt -> nicht-parametrische Tests)\n")

# ------------------------------------------------------------
# 6. Kruskal-Wallis-Test: Globaler Gruppeneffekt
# ------------------------------------------------------------
cat("\n--- Kruskal-Wallis-Test: Globaler Effekt von CONDITION ---\n")

kw <- kruskal.test(CUM_FEEDBACK_TIME_S ~ CONDITION, data = KumulDaten)
print(kw)

# Effektgröße: Eta-squared für Kruskal-Wallis
# eta^2 = (H - k + 1) / (n - k), wobei k = Anzahl Gruppen, n = Gesamt-N
k <- 3
n <- nrow(KumulDaten)
eta_sq <- (kw$statistic - k + 1) / (n - k)
cat(sprintf("  Effektgröße eta² = %.3f\n", eta_sq))
cat("  (Interpretation: klein >= .01, mittel >= .06, groß >= .14)\n")

# ------------------------------------------------------------
# 7. Post-hoc Mann-Whitney-U-Tests (paarweise Vergleiche)
# ------------------------------------------------------------
cat("\n--- Post-hoc Mann-Whitney-U-Tests (paarweise, zweiseitig) ---\n")

paare <- list(
  c("FFT", "Complete Tree"),
  c("FFT", "Risk Table"),
  c("Complete Tree", "Risk Table")
)

ergebnisse <- data.frame()

for (paar in paare) {
  g1 <- KumulDaten$CUM_FEEDBACK_TIME_S[KumulDaten$CONDITION == paar[1]]
  g2 <- KumulDaten$CUM_FEEDBACK_TIME_S[KumulDaten$CONDITION == paar[2]]
  mw <- wilcox.test(g1, g2, exact = FALSE)

  # Effektgröße r = Z / sqrt(N)
  z_val <- qnorm(mw$p.value / 2)
  r_eff <- abs(z_val) / sqrt(length(g1) + length(g2))

  ergebnisse <- rbind(ergebnisse, data.frame(
    Vergleich  = paste(paar[1], "vs.", paar[2]),
    W          = round(mw$statistic, 1),
    p          = round(mw$p.value, 4),
    r          = round(r_eff, 3),
    sig        = ifelse(mw$p.value < .001, "***",
                 ifelse(mw$p.value < .01,  "**",
                 ifelse(mw$p.value < .05,  "*",
                 ifelse(mw$p.value < .1,   ".", "n.s."))))
  ))
}

print(ergebnisse)
cat("\nSignifikanzcodes: *** p<.001  ** p<.01  * p<.05  . p<.1\n")
cat("Effektgröße r: klein >= .1, mittel >= .3, groß >= .5\n")

# ------------------------------------------------------------
# 8. Visualisierung
# ------------------------------------------------------------

# Violin + Boxplot
p <- ggplot(KumulDaten, aes(x = CONDITION, y = CUM_FEEDBACK_TIME_S / 60,
                             fill = CONDITION)) +
  geom_violin(alpha = 0.4, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.07, alpha = 0.3, size = 1.2, color = "grey30") +
  scale_fill_manual(values = c("FFT"           = "#2196F3",
                               "Complete Tree" = "#F44336",
                               "Risk Table"    = "#4CAF50")) +
  labs(
    title    = "Cumulative Feedback Viewing Time by Condition",
    subtitle = "Total per subject; points = individual values",
    x        = "Condition",
    y        = "Cumulative Feedback Viewing Time (minutes)",
    fill     = "Condition"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

ggsave("Plot_H1b_CumulativeFBTbyCondition_en.png", p, width = 7, height = 5, dpi = 150)
cat("\nAbbildung gespeichert: Plot_H1b_CumulativeFBTbyCondition_en.png\n")

# ------------------------------------------------------------
# 9. Hinweise zur Interpretation
# ------------------------------------------------------------
cat("
HINWEISE ZUR INTERPRETATION:
- Die abhängige Variable ist die Summe der Feedback-Viewing-Times
  pro Versuchsperson über alle Trials (in Sekunden).
- Da die kumulativen Werte rechtsschief verteilt sind (Schiefe ~1.1-1.5),
  werden nicht-parametrische Tests verwendet.
- Der Kruskal-Wallis-Test prüft den globalen Gruppeneffekt.
- Die paarweisen Mann-Whitney-U-Tests entsprechen dem Vorgehen
  der ursprünglichen Revision und sind direkt vergleichbar.
- Effektgröße r wird als |Z| / sqrt(N) berechnet.
- Diese Analyse ist auf Personenebene und ergänzt das
  Mixed-Effects-Modell auf Trial-Ebene.
")
