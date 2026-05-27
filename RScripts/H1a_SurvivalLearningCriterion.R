# ============================================================
# Hypothese 1a: Lerngeschwindigkeit
# Abhängige Variable: Anzahl Trials bis zum Lernkriterium
#                     (MAX_TRIAL pro Versuchsperson)
# Methode: Survival-Analyse (Cox-Regression) +
#           post-hoc Mann-Whitney-U-Tests
# Entspricht dem Vorgehen im revidierten Manuskript
# ============================================================

library(survival)
library(survminer)
library(ggplot2)

# ------------------------------------------------------------
# 1. Daten einlesen
# ------------------------------------------------------------
FeedbackDaten <- learning_timings_feedback_neu

# ------------------------------------------------------------
# 2. Variablen vorbereiten
# ------------------------------------------------------------

FeedbackDaten$CONDITION <- factor(
  FeedbackDaten$CONDITION,
  levels = c(1, 2, 3),
  labels = c("FFT", "Complete Tree", "Risk Table")
)
FeedbackDaten$SUBJECT <- factor(FeedbackDaten$SUBJECT)

# Anzahl Trials bis Lernkriterium pro Person
TrialDaten <- aggregate(TRIAL ~ SUBJECT + CONDITION, data = FeedbackDaten, FUN = max)
colnames(TrialDaten)[3] <- "MAX_TRIAL"

# ANT_E pro Person hinzufügen
ant_per_subj <- FeedbackDaten[!duplicated(FeedbackDaten$SUBJECT),
                               c("SUBJECT","ANT_E","ANT_E_BUCKET")]
TrialDaten <- merge(TrialDaten, ant_per_subj, by = "SUBJECT")

cat("Datenbeschreibung:\n")
cat("  Versuchspersonen:", nrow(TrialDaten), "\n")
cat("\n  N je Bedingung:\n")
print(table(TrialDaten$CONDITION))

# ------------------------------------------------------------
# 3. Deskriptive Statistik
# ------------------------------------------------------------
cat("\n--- Deskriptive Statistik (Trials bis Lernkriterium) ---\n")
deskriptiv <- aggregate(MAX_TRIAL ~ CONDITION, data = TrialDaten,
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

# ------------------------------------------------------------
# 4. Survival-Analyse
# ------------------------------------------------------------
# Das Lernkriterium kann nach minimal 40 Trials erreicht werden.
# "Ereignis" = Lernkriterium erreicht (alle Versuchspersonen
# haben es erreicht, daher keine Zensierung).
# Wir modellieren die Zeit (Trials) bis zum Ereignis.

# Ereignis-Indikator: alle haben das Kriterium erreicht
TrialDaten$EVENT <- 1

# Survival-Objekt
surv_obj <- Surv(time = TrialDaten$MAX_TRIAL, event = TrialDaten$EVENT)

# -- Kaplan-Meier-Kurven je Bedingung --
km_fit <- survfit(surv_obj ~ CONDITION, data = TrialDaten)

cat("\n--- Kaplan-Meier-Überblick ---\n")
print(summary(km_fit)$table)

# -- Cox-Regression: Zeitfenster 40-60 --
cat("\n--- Cox-Regression: Trials 40-60 ---\n")
dat_60 <- TrialDaten
dat_60$TIME  <- pmin(dat_60$MAX_TRIAL, 60)
dat_60$EVENT <- as.integer(dat_60$MAX_TRIAL <= 60)

cox_60 <- coxph(Surv(TIME, EVENT) ~ CONDITION, data = dat_60)
print(summary(cox_60))

cat("\nHazard Ratios (HR > 1 = schneller, d.h. früher Lernkriterium erreicht):\n")
print(round(exp(coef(cox_60)), 3))

# -- Cox-Regression: Zeitfenster 40-55 --
cat("\n--- Cox-Regression: Trials 40-55 ---\n")
dat_55 <- TrialDaten
dat_55$TIME  <- pmin(dat_55$MAX_TRIAL, 55)
dat_55$EVENT <- as.integer(dat_55$MAX_TRIAL <= 55)

cox_55 <- coxph(Surv(TIME, EVENT) ~ CONDITION, data = dat_55)
print(summary(cox_55))

cat("\nHazard Ratios:\n")
print(round(exp(coef(cox_55)), 3))

# -- Log-Rank-Test (globaler Gruppenvergleich) --
cat("\n--- Log-Rank-Test: Globaler Gruppenvergleich ---\n")
logrank <- survdiff(surv_obj ~ CONDITION, data = TrialDaten)
print(logrank)

# Effektgröße: Pseudo-R² (Nagelkerke) für Cox-Modell
r2_cox <- function(model) {
  n <- model$n
  r2 <- 1 - exp(-model$score / n)
  r2_max <- 1 - exp(2 * model$loglik[1] / n)
  round(r2 / r2_max, 3)
}
cat(sprintf("\nNagelkerke R² (Cox 40-60): %.3f\n", r2_cox(cox_60)))
cat(sprintf("Nagelkerke R² (Cox 40-55): %.3f\n", r2_cox(cox_55)))

# ------------------------------------------------------------
# 5. Post-hoc Mann-Whitney-U-Tests
# ------------------------------------------------------------
cat("\n--- Post-hoc Mann-Whitney-U-Tests (paarweise) ---\n")

paare <- list(
  c("FFT", "Complete Tree"),
  c("FFT", "Risk Table"),
  c("Complete Tree", "Risk Table")
)

ergebnisse <- data.frame()
for (paar in paare) {
  g1 <- TrialDaten$MAX_TRIAL[TrialDaten$CONDITION == paar[1]]
  g2 <- TrialDaten$MAX_TRIAL[TrialDaten$CONDITION == paar[2]]
  mw <- wilcox.test(g1, g2, exact = FALSE)
  z_val <- qnorm(mw$p.value / 2)
  r_eff <- abs(z_val) / sqrt(length(g1) + length(g2))
  ergebnisse <- rbind(ergebnisse, data.frame(
    Vergleich = paste(paar[1], "vs.", paar[2]),
    W = round(mw$statistic, 1),
    p = round(mw$p.value, 4),
    r = round(r_eff, 3),
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
# 6. Visualisierung: Kaplan-Meier-Kurven
# Beide Abbildungen beginnen erst ab Trial 40 (frühestmöglicher
# Zeitpunkt, zu dem das Lernkriterium erreicht werden kann).
# ------------------------------------------------------------

# Hilfsfunktion: KM-Plot für ein bestimmtes Zeitfenster
# Hinweis: ggsurvplot() gibt kein regulaeres ggplot-Objekt zurueck,
# daher kann ggsave() es nicht direkt verarbeiten. Stattdessen
# wird das PNG-Device manuell geoeffnet und geschlossen.
plot_km <- function(km_obj, data, xlim_max, title_text, dateiname) {

  p <- ggsurvplot(
    km_obj,
    data          = data,
    fun           = "event",         # Kumulativer Anteil mit erreichtem Kriterium
    conf.int      = TRUE,
    pval          = TRUE,
    xlim          = c(40, xlim_max), # X-Achse beginnt bei Trial 40
    break.time.by = 5,               # Tick-Abstaende alle 5 Trials
    legend.labs   = c("FFT", "Complete Tree", "Risk Table"),
    palette       = c("#2196F3", "#F44336", "#4CAF50"),
    xlab          = "Trial",
    ylab          = "Cumulative percentage of participants\nwho met the learning criterion",
    title         = title_text,
    ggtheme       = theme_minimal(base_size = 13),
    legend.title  = "Condition"
  )

  png(dateiname, width = 8, height = 5, units = "in", res = 150)
  print(p)
  dev.off()
  cat(sprintf("Abbildung gespeichert: %s\n", dateiname))
}

# Abbildung 1: Zeitfenster Trials 40-60
plot_km(
  km_obj    = km_fit,
  data      = TrialDaten,
  xlim_max  = 60,
  title_text = "Learning curves by condition: Time window trials 40-60 (Kaplan-Meier)",
  dateiname  = "Plot_H1a_Survival_40_60_EN.png"
)

# Abbildung 2: Zeitfenster Trials 40-55
plot_km(
  km_obj    = km_fit,
  data      = TrialDaten,
  xlim_max  = 55,
  title_text = "Learning curves by condition: Time window trials 40–55 (Kaplan-Meier)",
  dateiname  = "Plot_H1a_Survival_40_55_EN.png"
)


# Abbildung 3: Alle Trials ab 40 (vollständiger Verlauf)
plot_km(
  km_obj    = km_fit,
  data      = TrialDaten,
  xlim_max  = max(TrialDaten$MAX_TRIAL),
  title_text = "Learning curves by condition: All trials with n ≥ 40 (Kaplan-Meier)",
  dateiname  = "Plot_H1a_Survival_40_all_EN.png"
)
# ------------------------------------------------------------
# 7. Hinweise zur Interpretation
# ------------------------------------------------------------
cat("
HINWEISE ZUR INTERPRETATION:
- Hazard Ratio (HR) > 1: Bedingung erreicht das Lernkriterium
  SCHNELLER als die Referenz (FFT).
  HR < 1: Bedingung erreicht es LANGSAMER.
- Achtung: In der Survival-Analyse für Lernkurven bedeutet
  ein hohes HR, dass in dieser Bedingung mehr Personen früher
  das Kriterium erreichen (= günstig).
- Das Zeitfenster 40-55 ist restriktiver und erfasst den
  Bereich, in dem die Gruppenunterschiede am größten sind.
- Der Log-Rank-Test prüft den globalen Gruppenunterschied
  über alle Trials.
- Nagelkerke R² ist ein Maß für die Modellgüte des
  Cox-Modells (0 = kein Effekt, 1 = perfekte Vorhersage).
- Die paarweisen Mann-Whitney-U-Tests ergänzen die
  Survival-Analyse auf Personenebene.
")
