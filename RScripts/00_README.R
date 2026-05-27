# ============================================================
# ÜBERSICHT ALLER ANALYSE-SKRIPTE
# Replikation der statistischen Analysen aus dem revidierten
# Manuskript (Revision_NK_May_04.docx)
# ============================================================
#
# VORAUSSETZUNGEN
# ---------------
# Folgende R-Pakete müssen installiert sein:
#
#   install.packages(c("lme4", "lmerTest", "survival",
#                      "survminer", "ggplot2", "scales"))
#
# Alle Skripte setzen voraus, dass die Rohdaten bereits
# in der R-Umgebung geladen sind:
#
#   learning_timings_feedback_neu  -> Feedback-Analyse (H1a, H1b)
#   learning_timings_errors_neu    -> Fehlerraten-Analyse (H2)
#   surprise_recall                -> Memory-Test-Analyse (H3)
#
# HINWEIS ZU STICHPROBENGROSSEN
# ------------------------------
# Die Stichprobengrößen unterscheiden sich zwischen den
# Datensätzen, da der Outlier-Ausschluss je abhängiger
# Variable separat durchgeführt wurde:
#
#   learning_timings_feedback_neu: N = 153 (nach Ausschluss)
#   learning_timings_errors_neu:   N = 174 (nach Ausschluss)
#   surprise_recall:               N = 186 (alle rekrutierten
#                                  VP, kein Ausschluss)
#
# Dies sollte im Rebuttal explizit begründet werden.
#
# ============================================================
#
# SKRIPTE UND HYPOTHESEN
# ----------------------
#
# 1. H1a_SurvivalLearningCriterion.R
#    Hypothese 1a: Lerngeschwindigkeit
#    AV: Anzahl Trials bis Lernkriterium (MAX_TRIAL)
#    Methode: Survival-Analyse (Cox-Regression, Log-Rank-Test)
#             + post-hoc Mann-Whitney-U-Tests
#    Zeitfenster: Trials 40-55 und 40-60 (wie im Manuskript)
#    Effektgröße: Hazard Ratio, Nagelkerke R², r
#    Daten: learning_timings_feedback_neu (N = 153)
#
# 2. H1b_CumulativeFeedbackViewingTimes.R
#    Hypothese 1b: Feedback Viewing Time (Personenebene)
#    AV: Kumulative Betrachtungszeit pro Person (Sekunden)
#    Methode: Kruskal-Wallis-Test + Mann-Whitney-U-Tests
#    Effektgröße: eta², r
#    Daten: learning_timings_feedback_neu (N = 153)
#
# 3. H1b_MixedModelFeedbackViewingTimes.R
#    Hypothese 1b: Feedback Viewing Time (Trial-Ebene)
#    AV: Log(FEEDBACK_TIME) pro Trial
#    Methode: Lineares gemischtes Modell (lme4/lmerTest)
#             mit gekreuzten Zufallseffekten für SUBJECT
#             und ITEM (Stimuluskombination)
#    Effektgröße: exp(b) als multiplikativer Faktor,
#                 Varianzdekomposition (ICC)
#    Daten: learning_timings_feedback_neu (N = 153)
#    HINWEIS: Ergänzende Analyse; Reviewer 3 hatte diese
#             Methode explizit gefordert (Judd et al., 2012)
#
# 4. H2_ErrorRate.R
#    Hypothese 2: Anwendungsgenauigkeit (Personenebene)
#    AV: Fehlerrate pro Person (Anteil fehlerhafter Trials)
#    Methode: Kruskal-Wallis-Test + Mann-Whitney-U-Tests
#    Effektgröße: eta², r
#    Daten: learning_timings_errors_neu (N = 174)
#
# 5. H2_MixedModelErrorRate.R
#    Hypothese 2: Anwendungsgenauigkeit (Trial-Ebene)
#    AV: ERROR pro Trial (binär 0/1)
#    Methode: GLMM mit logistischer Linkfunktion (lme4)
#             mit gekreuzten Zufallseffekten für SUBJECT
#             und ITEM (Stimuluskombination)
#    Effektgröße: Odds Ratio, Varianzdekomposition (ICC)
#    Daten: learning_timings_errors_neu (N = 174)
#    HINWEIS: Ergänzende Analyse; Reviewer 3 hatte diese
#             Methode explizit gefordert (Judd et al., 2012)
#
# 6. H3_MemoryTestScores.R
#    Hypothese 3: Gedächtnisleistung + ANT-E-Interaktion
#    AV: SMT_Score pro Person (kontinuierlich, 0.0-1.0)
#        = Anteil korrekt klassifizierter Kombination
#    Methode: Kruskal-Wallis-Test + Mann-Whitney-U-Tests
#             + Lineare Regression für ANT_E x CONDITION
#    Effektgröße: eta², r, R² (Regression)
#    Daten: surprise_recall (N = 186, alle VP)
#    HINWEIS: Starker Deckeneffekt (49-94% perfekte Scores)
#             -> nicht-parametrische Tests zwingend
#
# ============================================================
#
# EMPFOHLENE REIHENFOLGE DER AUSFÜHRUNG
# --------------------------------------
#
#   source("H1a_survival_lernkriterium.R")
#   source("kumulative_feedback_time.R")
#   source("mixed_model_feedback_time.R")
#   source("H2_fehlerrate.R")
#   source("mixed_model_error_rate.R")
#   source("H3_memory_test.R")
#
# ============================================================
#
# ZUORDNUNG SKRIPTE -> KRITIKPUNKTE REVIEWER 3
# ---------------------------------------------
#
# Reviewer 3 hatte folgende methodische Änderungen gefordert:
#
#   [✅] Survival-Analyse (H1a)
#        -> H1a_survival_lernkriterium.R
#
#   [✅] Kumulative Betrachtungszeiten + nicht-parametrische
#        Tests (H1b)
#        -> kumulative_feedback_time.R
#
#   [✅] Mixed-Effects-Modelle für H1b und H2
#        (Judd et al., 2012, 2017)
#        -> mixed_model_feedback_time.R
#        -> mixed_model_error_rate.R
#
#   [✅] IQR-basierte Outlier-Entfernung (Ande, 2022)
#        -> in allen Skripten dokumentiert
#
#   [✅] Nicht-parametrische Tests statt ANOVA
#        -> alle Personenebene-Skripte verwenden
#           Kruskal-Wallis + Mann-Whitney-U
#
# ============================================================

cat("README geladen. Bitte Skripte in der oben angegebenen
Reihenfolge ausfuehren. Alle Skripte erfordern, dass die
Rohdaten bereits als R-Objekte in der Umgebung vorliegen:

  learning_timings_feedback_neu  (N = 153)
  learning_timings_errors_neu    (N = 174)
  surprise_recall                (N = 186)\n")
