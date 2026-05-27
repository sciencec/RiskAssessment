# ============================================================
# Hypothese 3: Gedächtnisleistung (Memory Test Score)
# Abhängige Variable: SMT_Score pro Versuchsperson
#                     (kontinuierlich, 0.0 - 1.0)
#                     = Anteil der Stimuluskombinationen,
#                     die der reproduzierte Entscheidungsbaum
#                     korrekt klassifiziert
# Methode: Kruskal-Wallis-Test + post-hoc Mann-Whitney-U-Tests
#           + Interaktion mit ANT_E (Berlin Numeracy Test)
# Entspricht dem Vorgehen im revidierten Manuskript
# ============================================================

library(ggplot2)

# ------------------------------------------------------------
# 1. Daten einlesen
# ------------------------------------------------------------
# Daten direkt aus der R-Umgebung verwenden
RecallDaten <- surprise_recall

# ------------------------------------------------------------
# 2. Variablen vorbereiten
# ------------------------------------------------------------

# Spaltennamen bereinigen (BOM-Zeichen in Participant #)
colnames(RecallDaten) <- c("SUBJECT","CONDITION","ANT_E","ANT_E_BUCKET","SMT_Score")

# CONDITION als Faktor mit sprechenden Labels (FFT = Referenzkategorie)
RecallDaten$CONDITION <- factor(
  RecallDaten$CONDITION,
  levels = c(1, 2, 3),
  labels = c("FFT", "Complete Tree", "Risk Table")
)

# ANT_E_BUCKET als Faktor
RecallDaten$ANT_E_BUCKET <- factor(
  RecallDaten$ANT_E_BUCKET,
  levels = c(1, 2, 3),
  labels = c("low (0-1)",
             "medium (2)",
             "high (3-4)")
)

# SMT_Score in Prozent (0-100) für bessere Lesbarkeit
RecallDaten$SMT_PCT <- RecallDaten$SMT_Score * 100

cat("Datenbeschreibung:\n")
cat("  Versuchspersonen gesamt:", nrow(RecallDaten), "\n")
cat("\n  N je Bedingung:\n")
print(table(RecallDaten$CONDITION))

# ------------------------------------------------------------
# 3. Deskriptive Statistik
# ------------------------------------------------------------
cat("\n--- Deskriptive Statistik (SMT_Score in %) ---\n")
deskriptiv <- aggregate(SMT_PCT ~ CONDITION, data = RecallDaten,
  FUN = function(x) c(
    N        = length(x),
    Mean     = round(mean(x), 1),
    Median   = round(median(x), 1),
    SD       = round(sd(x), 1),
    Min      = round(min(x), 1),
    Max      = round(max(x), 1),
    Perfekt  = round(mean(x == 100) * 100, 1)
  )
)
print(do.call(data.frame, deskriptiv))

cat("\nAnteil perfekter Scores (SMT = 100%) je Bedingung:\n")
print(tapply(RecallDaten$SMT_Score, RecallDaten$CONDITION,
             function(x) paste0(round(mean(x == 1.0) * 100, 1), "%")))

# ------------------------------------------------------------
# 4. Voraussetzungsprüfung: Normalverteilung
# ------------------------------------------------------------
cat("\n--- Shapiro-Wilk-Test auf Normalverteilung je Bedingung ---\n")
for (bed in levels(RecallDaten$CONDITION)) {
  x  <- RecallDaten$SMT_PCT[RecallDaten$CONDITION == bed]
  sw <- shapiro.test(x)
  cat(sprintf("  %s: W = %.4f, p = %.4f\n", bed, sw$statistic, sw$p.value))
}
cat("  (Starke negative Schiefe durch Deckeneffekt -> nicht-parametrische Tests)\n")

# ------------------------------------------------------------
# 5. Kruskal-Wallis-Test: Globaler Gruppeneffekt
# ------------------------------------------------------------
cat("\n--- Kruskal-Wallis-Test: Globaler Effekt von CONDITION ---\n")
kw <- kruskal.test(SMT_PCT ~ CONDITION, data = RecallDaten)
print(kw)

# Effektgröße eta²
k      <- 3
n      <- nrow(RecallDaten)
eta_sq <- (kw$statistic - k + 1) / (n - k)
cat(sprintf("  Effektgröße eta² = %.3f\n", eta_sq))
cat("  (Interpretation: klein >= .01, mittel >= .06, groß >= .14)\n")

# ------------------------------------------------------------
# 6. Post-hoc Mann-Whitney-U-Tests (paarweise)
# ------------------------------------------------------------
cat("\n--- Post-hoc Mann-Whitney-U-Tests (paarweise, zweiseitig) ---\n")

paare <- list(
  c("FFT", "Complete Tree"),
  c("FFT", "Risk Table"),
  c("Complete Tree", "Risk Table")
)

ergebnisse <- data.frame()
for (paar in paare) {
  g1  <- RecallDaten$SMT_PCT[RecallDaten$CONDITION == paar[1]]
  g2  <- RecallDaten$SMT_PCT[RecallDaten$CONDITION == paar[2]]
  mw  <- wilcox.test(g1, g2, exact = FALSE)
  z   <- qnorm(mw$p.value / 2)
  r   <- abs(z) / sqrt(length(g1) + length(g2))
  ergebnisse <- rbind(ergebnisse, data.frame(
    Vergleich = paste(paar[1], "vs.", paar[2]),
    W   = round(mw$statistic, 1),
    p   = round(mw$p.value, 6),
    r   = round(r, 3),
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
# 7. Interaktion: ANT_E x CONDITION
# ------------------------------------------------------------
cat("\n--- Interaktion: ANT_E_BUCKET x CONDITION (SMT_Score in %) ---\n")

# Deskriptive Übersicht
recall_ant <- aggregate(SMT_PCT ~ CONDITION + ANT_E_BUCKET,
                        data = RecallDaten, FUN = mean)
recall_ant$SMT_PCT <- round(recall_ant$SMT_PCT, 1)
cat("\nMittlerer SMT-Score (%) nach Bedingung und ANT_E-Bucket:\n")
print(reshape(recall_ant,
              idvar = "CONDITION", timevar = "ANT_E_BUCKET",
              direction = "wide"))

# Kruskal-Wallis je ANT_E_BUCKET-Gruppe
cat("\nKruskal-Wallis je ANT_E_BUCKET-Gruppe:\n")
for (bucket in levels(RecallDaten$ANT_E_BUCKET)) {
  sub    <- RecallDaten[RecallDaten$ANT_E_BUCKET == bucket, ]
  kw_sub <- kruskal.test(SMT_PCT ~ CONDITION, data = sub)
  cat(sprintf("  %s: H(%d) = %.3f, p = %.4f\n",
              bucket, kw_sub$parameter, kw_sub$statistic, kw_sub$p.value))
}

# Lineare Regression mit Interaktionsterm (ANT_E metrisch zentriert)
cat("\n--- Lineare Regression: SMT_Score ~ CONDITION * ANT_E ---\n")
RecallDaten$ANT_E_C <- scale(RecallDaten$ANT_E, center = TRUE, scale = FALSE)
lm_inter <- lm(SMT_Score ~ CONDITION * ANT_E_C, data = RecallDaten)
print(summary(lm_inter))

# Modellvergleich: Mit vs. ohne Interaktion
lm_main <- lm(SMT_Score ~ CONDITION + ANT_E_C, data = RecallDaten)
cat("\n--- F-Test: Interaktionsterm ---\n")
print(anova(lm_main, lm_inter))

# ------------------------------------------------------------
# 8. Visualisierung
# ------------------------------------------------------------

# Violin + Boxplot je Bedingung
p1 <- ggplot(RecallDaten, aes(x = CONDITION, y = SMT_PCT, fill = CONDITION)) +
  geom_violin(alpha = 0.4, trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.07, alpha = 0.3, size = 1.2, color = "grey30") +
  scale_fill_manual(values = c("FFT"           = "#2196F3",
                               "Complete Tree" = "#F44336",
                               "Risk Table"    = "#4CAF50")) +
  labs(
    title    = "Memory test score by Condition",
    subtitle = "SMT score in %",
    x = "Condition", y = "SMT score (%)", fill = "Condition"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")
ggsave("Plot_H3_SMTbyCondition_EN.png", p1, width = 7, height = 5, dpi = 150)

# Interaktion: ANT_E_BUCKET x CONDITION
p2 <- ggplot(recall_ant,
             aes(x = ANT_E_BUCKET, y = SMT_PCT,
                 color = CONDITION, group = CONDITION)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_color_manual(values = c("FFT"           = "#2196F3",
                                "Complete Tree" = "#F44336",
                                "Risk Table"    = "#4CAF50")) +
  labs(
    title    = "Interaction: Memory test score x Berlin Numeracy Test",
    subtitle = "Mean SMT score (%) by Condition and ANT_E",
    x = "Berlin Numeracy Test (ANT_E)",
    y = "SMT score (%)",
    color = "Condition"
  ) +
  theme_minimal(base_size = 13) +
  theme(axis.text.x = element_text(angle = 15, hjust = 1))
ggsave("Plot_H3_Interaction_SMT_ANT-E_EN.png", p2, width = 8, height = 5, dpi = 150)

cat("\nAbbildungen gespeichert: Plot_H3_SMTbyCondition_EN.png, Plot_H3_SMTbyCondition_EN.png\n")

# ------------------------------------------------------------
# 9. Hinweise zur Interpretation
# ------------------------------------------------------------
cat("
HINWEISE ZUR INTERPRETATION:
- SMT_Score ist der Anteil der Stimuluskombinationen, die der
  reproduzierte Entscheidungsbaum korrekt klassifiziert
  (0 = alle falsch, 1.0 = alle richtig; entspricht 0-100%).
- Die starke negative Schiefe (Deckeneffekt: 49-94% perfekte
  Scores je Bedingung) macht nicht-parametrische Tests
  zwingend erforderlich.
- Kruskal-Wallis prüft den globalen Gruppeneffekt;
  eta² als Effektgröße: klein >= .01, mittel >= .06, groß >= .14.
- Mann-Whitney-U für paarweise Vergleiche;
  r als Effektgröße: klein >= .1, mittel >= .3, groß >= .5.
- Die lineare Regression für den Interaktionsterm ist eine
  Annäherung; eine Beta-Regression wäre für proportionale
  Daten [0,1] streng genommen angemessener.
- Hypothese 3: FFT-Nutzer erzielen höhere SMT-Scores als
  Complete-Tree- und Risk-Table-Nutzer.
- Interaktionshypothese: Der FFT-Vorteil ist robust gegenüber
  Unterschieden im ANT_E – die FFT-Linie verläuft in der
  Interaktionsgrafik flacher als die anderen Bedingungen.
- ACHTUNG: N = 186 in dieser Datei vs. N = 153 in den anderen
  Dateien. Diese Datei enthält alle ursprünglich rekrutierten
  Versuchspersonen vor dem Outlier-Ausschluss für die anderen
  abhängigen Variablen. Dies sollte im Rebuttal explizit
  erwähnt und die unterschiedlichen N begründet werden.
")
