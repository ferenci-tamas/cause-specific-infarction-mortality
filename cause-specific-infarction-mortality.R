library(data.table)
library(ggplot2)
theme_set(theme_bw() +
            theme(legend.position = "bottom",
                  text = element_text(family = "Helvetica")))

RawData <- readRDS("RawData.rds")

prop.table(table(RawData$X.92..Diagnózis))
RawData[, .(median(time)) , .(X.92..Diagnózis)]

table(RawData$KSH.elsődleges.halál.ok)

RawData$event2 <- relevel(as.factor(ifelse(
  !RawData$event, "Alive", RawData$KSH.elsődleges.halál.ok)),
  ref = "Alive")
RawData$event3 <- relevel(as.factor(ifelse(
  !RawData$event, "Alive", "All-cause mortality")),
  ref = "Alive")

RawData$Mort30day <- RawData$event & (RawData$time <= 30)
RawData$Mort1year <- RawData$event & (RawData$time <= 365)

fwrite(data.table(summarycsv(
  AGE + SEX + PRIOR_MI + HF + HT + PRIOR_STROKE + DM + PAD +
    COPD + CKD + PCI + Mort30day + Mort1year + event +
    KSH.elsődleges.halál.ok ~ X.92..Diagnózis, data = RawData,
  pctdig = 1, overall = TRUE), keep.rownames = TRUE),
  "Table1.csv", dec = ",", sep =";", bom = TRUE,
  col.names = FALSE)

cif <- data.table(rbind(
  tidycmprsk::cuminc(survival::Surv(time, event2) ~ 1,
                     data = RawData)$tidy,
  tidycmprsk::cuminc(survival::Surv(time, event3) ~ 1,
                     data = RawData)$tidy))

fwrite(dcast(cif[
  time%in%c(30, 90, 365),
  .(time, outcome, CIF = paste0(round(estimate*100, 2), " (",
                                round(conf.low*100, 2), " - ",
                                round(conf.high*100, 2), ")"))],
  time ~ outcome, value.var = "CIF"),
  "Table2.csv", dec = ",", sep =";", bom = TRUE)

ggplot(cif,
       aes(x = time, y = estimate*100, ymin = conf.low*100,
           ymax = conf.high*100, color = outcome, fill = outcome,
           group = outcome)) + geom_line() +
  geom_ribbon(alpha = 0.2, color = NA) +
  theme(legend.position = "bottom") +
  labs(x = "Time [day]", y = "Cumulative incidence [%]",
       color = "", fill = "")
ggsave("Figure1.pdf", width = 16, height = 9, device = cairo_pdf)

fit <- tidycmprsk::crr(
  survival::Surv(time, event2) ~
    X.92..Diagnózis + AGE + SEX + PRIOR_MI + HF + HT +
    PRIOR_STROKE + DM + PAD + PCI, data = RawData,
  failcode = "AMI")

fitCox <- survival::coxph(
  survival::Surv(time, event) ~
    X.92..Diagnózis + AGE + SEX + PRIOR_MI + HF + HT +
    PRIOR_STROKE + DM + PAD + PCI, data = RawData)

fitres <- rbind(
  data.table(fit$tidy)[
    , .(var = term, HR = exp(estimate), conf.low = exp(conf.low),
        conf.high = exp(conf.high),
        type = "Cause-specific mortality")],
  with(summary(fitCox), setNames(data.table(
    rownames(coefficients), coefficients[,"exp(coef)"],
    conf.int[, c("lower .95", "upper .95")],
    "All-cause mortality"), c("var", "HR", "conf.low",
                              "conf.high", "type"))))

fitres <- merge(fitres, data.table(
  var = c("X.92..DiagnózisSTEMI", "AGE", "SEXNő", "PRIOR_MIIgen",
          "HFIgen", "HTIgen", "PRIOR_STROKEIgen", "DMIgen",
          "PADIgen", "PCIIgen"),
  varlabel = c("STEMI", "Age +1 year", "Female sex", "Prior MI",
               "Heart failure", "Hypertension", "Prior stroke",
               "DM", "PAD", "PCI")), by = "var")

fitres$varlabel <- factor(fitres$varlabel, levels = c(
  "Age +1 year", "Female sex", "STEMI", "Prior MI",
  "Heart failure", "Prior stroke", "Hypertension", "DM",  "PAD",
  "PCI"))

ggplot(fitres,
       aes(y = varlabel, x = HR, xmin = conf.low,
           xmax = conf.high, group = type, color = type)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(width = 0.2,
                position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 1, color = "red") +
  scale_y_discrete(limits = rev) +
  labs(x = "Hazard ratio (HR)", y = "", color = "")
ggsave("Figure2.pdf", width = 16, height = 9, device = cairo_pdf)

gtsummary::as_hux_xlsx(gtsummary::add_global_p(
  gtsummary::tbl_regression(fit, exp = TRUE), keep = TRUE),
  "Modell.xlsx")

sort(table(
  RawData[X.92..Diagnózis == "NSTEMI" & event == TRUE &
            event2 == "Cardiac other than AMI"]$Diagnózis))

sort(table(substring(
  RawData[X.92..Diagnózis == "NSTEMI" & event == TRUE &
            event2 == "Cardiac other than AMI"]$Diagnózis, 1, 3)))
sort(prop.table(table(substring(
  RawData[X.92..Diagnózis == "NSTEMI" & event == TRUE &
            event2 == "Cardiac other than AMI"]$Diagnózis,
  1, 3))))*100
