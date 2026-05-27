# ============================================================
# Hypothese 2: Anwendungsgenauigkeit (Error Rate)
# Abhängige Variable: Fehlerrate pro Versuchsperson
#                     (aggregiert über alle Trials)
# Methode: Kruskal-Wallis-Test + post-hoc Mann-Whitney-U-Tests
# Entspricht dem Vorgehen im revidierten Manuskript
# ============================================================

library(ggplot2)

# ------------------------------------------------------------
# 1. Daten einlesen
# ------------------------------------------------------------
ErrorDaten <- learning_timings_errors_neu

# ------------------------------------------------------------
# 2. Variablen vorbereiten
# ------------------------------------------------------------

# Leere Zeile entfernen
ErrorDaten <- ErrorDaten[!is.na(ErrorDaten$SUBJECT), ]

ErrorDaten$CONDITION <- factor(
  ErrorDaten$CONDITION,
  levels = c(1, 2, 3),
  labels = c("FFT", "Complete Tree", "Risk Table")
)
ErrorDaten$SUBJECT <- factor(ErrorDaten$SUBJECT)
ErrorDaten$ERROR   <- as.integer(ErrorDaten$ERROR)

# Fehlerrate pro Person: Anteil fehlerhafter Trials
PersonDaten <- aggregate(ERROR ~ SUBJECT + CONDITION, data = ErrorDaten, FUN = mean)
colnames(PersonDaten)[3] <- "ERROR_RATE"

# Anzahl Trials pro Person als Zusatzinfo
n_trials <- aggregate(ERROR ~ SUBJECT + CONDITION, data = ErrorDaten, FUN = length)
colnames(n_trials)[3] <- "N_TRIALS"
PersonDaten <- merge(PersonDaten, n_trials, by = c("SUBJECT","CONDITION"))

# ANT_E pro Person hinzufügen
ant_per_subj <- ErrorDaten[!duplicated(ErrorDaten$SUBJECT),
                            c("SUBJECT","ANT_E","ANT_E_BUCKET")]
PersonDaten <- merge(PersonDaten, ant_per_subj, by = "SUBJECT")

cat("Datenbeschreibung:\n")
cat("  Versuchspersonen:", nrow(PersonDaten), "\n")
cat("\n  N je Bedingung:\n")
print(table(PersonDaten$CONDITION))

# ------------------------------------------------------------
# 3. Deskriptive Statistik
# ------------------------------------------------------------
cat("\n--- Deskriptive Statistik (Fehlerrate pro Person, in %) ---\n")
deskriptiv <- aggregate(ERROR_RATE ~ CONDITION, data = PersonDaten,
  FUN = function(x) c(
    N      = length(x),
    Mean   = round(mean(x) * 100, 1),
    Median = round(median(x) * 100, 1),
    SD     = round(sd(x) * 100, 1),
    Min    = round(min(x) * 100, 1),
    Max    = round(max(x) * 100, 1)
  )
)
print(do.call(data.frame, deskriptiv))

# ------------------------------------------------------------
# 4. Voraussetzungsprüfung: Normalverteilung
# ------------------------------------------------------------
cat("\n--- Shapiro-Wilk-Test auf Normalverteilung je Bedingung ---\n")
for (bed in levels(PersonDaten$CONDITION)) {
  x  <- PersonDaten$ERROR_RATE[PersonDaten$CONDITION == bed]
  sw <- shapiro.test(x)
  cat(sprintf("  %s: W = %.4f, p = %.4f\n", bed, sw$statistic, sw$p.value))
}
cat("  (Bei p < .05: Normalverteilung abgelehnt -> nicht-parametrische Tests)\n")

# ------------------------------------------------------------
# 5. Kruskal-Wallis-Test: Globaler Gruppeneffekt
# ------------------------------------------------------------
cat("\n--- Kruskal-Wallis-Test: Globaler Effekt von CONDITION ---\n")
kw <- kruskal.test(ERROR_RATE ~ CONDITION, data = PersonDaten)
print(kw)

# Effektgröße eta²
k    <- 3
n    <- nrow(PersonDaten)
eta_sq <- (kw$statistic - k + 1) / (n - k)
cat(sprintf("  Effektgröße eta² = %.3f\n", eta_sq))
cat("  (Interpretation: klein >= .01, mittel >= .06, groß >= .14)\n")

# ------------------------------------------------------------
# 6. Post-hoc Mann-Whitney-U-Tests
# ------------------------------------------------------------
cat("\n--- Post-hoc Mann-Whitney-U-Tests (paarweise, zweiseitig) ---\n")

paare <- list(
  c("FFT", "Complete Tree"),
  c("FFT", "Risk Table"),
  c("Complete Tree", "Risk Table")
)

ergebnisse <- data.frame()
for (paar in paare) {
  g1  <- PersonDaten$ERROR_RATE[PersonDaten$CONDITION == paar[1]]
  g2  <- PersonDaten$ERROR_RATE[PersonDaten$CONDITION == paar[2]]
  mw  <- wilcox.test(g1, g2, exact = FALSE)
  z   <- qnorm(mw$p.value / 2)
  r   <- abs(z) / sqrt(length(g1) + length(g2))
  ergebnisse <- rbind(ergebnisse, data.frame(
    Vergleich = paste(paar[1], "vs.", paar[2]),
    W  = round(mw$statistic, 1),
    p  = round(mw$p.value, 4),
    r  = round(r, 3),
    sig = ifelse(mw$p.value < .001, "***",
          ifelse(mw$p.value < .01,  "**",
          ifelse(mw$p.value < .05,  "*",
          ifelse(mw$p.value < .1,   ".", "n.s."))))
  ))
}
print(ergebnisse)
cat("\nSignifikanzcodes: *** p<.001  ** p<.01  * p<.05  . p<.1\n")
cat("Effektgröße r: klein >= .1, mittel >= .3, groß >= .5\n")

# ------------------------------------------------------------
# 7. Visualisierung
# ------------------------------------------------------------
p <- ggplot(PersonDaten, aes(x = CONDITION, y = ERROR_RATE * 100, fill = CONDITION)) +
  geom_violin(alpha = 0.4, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.07, alpha = 0.3, size = 1.2, color = "grey30") +
  scale_fill_manual(values = c("FFT"           = "#2196F3",
                               "Complete Tree" = "#F44336",
                               "Risk Table"    = "#4CAF50")) +
  labs(
    title    = "Error rate by Condition",
    subtitle = "Proportion of incorrect trials per subject",
    x = "Condition", y = "Error rate (%)", fill = "Condition"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

ggsave("Plot_H2_ErrorRateByCondition_EN.png", p, width = 7, height = 5, dpi = 150)
cat("\nAbbildung gespeichert: Plot_H2_ErrorRateByCondition_EN.png\n")

# ------------------------------------------------------------
# 8. Hinweise zur Interpretation
# ------------------------------------------------------------
cat("
HINWEISE ZUR INTERPRETATION:
- Die Fehlerrate ist der Anteil fehlerhafter Trials pro Person
  (0 = keine Fehler, 1 = alle Trials fehlerhaft).
- Kruskal-Wallis prüft den globalen Gruppeneffekt.
- eta² als Effektgröße: klein >= .01, mittel >= .06, groß >= .14.
- Effektgröße r = |Z| / sqrt(N) für Mann-Whitney-U:
  klein >= .1, mittel >= .3, groß >= .5.
- Diese Analyse entspricht dem Vorgehen auf Personenebene
  im revidierten Manuskript (Mann-Whitney-U-Tests auf
  aggregierten Fehlerraten).
- Das GLMM auf Trial-Ebene ist eine ergänzende, robustere
  Analyse, die Item- und Personenvarianz berücksichtigt.
")
