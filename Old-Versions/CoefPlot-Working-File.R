##################################
library("broom")
library("ggplot2")

coefs_all <- tidy(tc_pwp, conf.int = TRUE, exponentiate = F)
coefs_all$HazardRatio <- exp(coefs_all$estimate)
coefs_all$Variable <- c(
  "Peace w/ TC", "Fighting w/ TC", "Favorable Outcome for State", "Favorable Outcome for TC", "Absorbed",
  "Liberal Democracy Index", "State Area (logged)", "Mountains", "ELF", "TC Tally"
)

coefs1 <- tidy(tc1_pwp, conf.int = TRUE, exponentiate = F)
coefs1$HazardRatio <- exp(coefs1$estimate)
coefs1$Variable <- c("Liberal Democracy Index", "State Area (logged)", "Mountains", "ELF")

coefs23 <- tidy(tc23_pwp_lin20, conf.int = TRUE, exponentiate = F)
coefs23$HazardRatio <- exp(coefs23$estimate)
coefs23$Variable <- c(
  "Peace w/ TC", "Fighting w/ TC", "Favorable Outcome for State", "Favorable Outcome for TC", "Absorbed",
  "Liberal Democracy Index", "State Area (logged)", "Mountains", "ELF"
)

coefs45 <- tidy(tc45_pwp_lin20, conf.int = TRUE, exponentiate = F)
coefs45$HazardRatio <- exp(coefs45$estimate)
coefs45$Variable <- c(
  "Peace w/ TC", "Fighting w/ TC", "Favorable Outcome for State", "Favorable Outcome for TC", "Absorbed",
  "Liberal Democracy Index", "State Area (logged)", "Mountains", "ELF"
)

coefs6plus <- tidy(tc6_more_pwp_lin20, conf.int = TRUE, exponentiate = F)
coefs6plus$HazardRatio <- exp(coefs6plus$estimate)
coefs6plus$Variable <- c(
  "Peace w/ TC", "Fighting w/ TC", "Favorable Outcome for State", "Favorable Outcome for TC", "Absorbed",
  "Liberal Democracy Index", "State Area (logged)", "Mountains", "ELF"
)


combined <- data.frame(
  Variable = coefs_all$Variable,
  HazardRatio = coefs_all$HazardRatio,
  Beta = coefs_all$estimate,
  SE = summary(tc_pwp)$coef[, 4],
  TC = "All TCs"
)
combined$Variable <- factor(combined$Variable, levels = unique(as.character(combined$Variable)))

first <- data.frame(
  Variable = coefs1$Variable,
  HazardRatio = coefs1$HazardRatio,
  Beta = coefs1$estimate,
  SE = summary(tc1_pwp)$coef[, 4],
  TC = "First TC"
)

mod23 <- data.frame(
  Variable = coefs23$Variable,
  HazardRatio = coefs23$HazardRatio,
  Beta = coefs23$estimate,
  SE = summary(tc23_pwp_lin20)$coef[, 4],
  TC = "2-3"
)

mod45 <- data.frame(
  Variable = coefs45$Variable,
  HazardRatio = coefs45$HazardRatio,
  Beta = coefs45$estimate,
  SE = summary(tc45_pwp_lin20)$coef[, 4],
  TC = "4-5"
)

mod6plus <- data.frame(
  Variable = coefs6plus$Variable,
  HazardRatio = coefs6plus$HazardRatio,
  Beta = coefs6plus$estimate,
  SE = summary(tc6_more_pwp_lin20)$coef[, 4],
  TC = "6+"
)


allmod <- data.frame(rbind(combined, first, mod23, mod45, mod6plus))
# allmod <- subset(allmod, Variable!="TC Tally") #TC Tally Removed!!!!!!
allmod$Variable <- factor(allmod$Variable, levels = unique(as.character(allmod$Variable)))


interval1 <- -qnorm((1 - 0.9) / 2) # 90% multiplier
interval2 <- -qnorm((1 - 0.95) / 2) # 95% multiplier


allmod$TC <- factor(allmod$TC, levels = c("All TCs", "First TC", "2-3", "4-5", "6+"))
allmod$TC <- factor(allmod$TC, levels = c("6+", "4-5", "2-3", "First TC", "All TCs"))

pwp_fig <- ggplot(allmod, aes(colour = TC))
pwp_fig <- pwp_fig + scale_color_viridis_d(option = "C", breaks = c("All TCs", "First TC", "2-3", "4-5", "6+"), name = "TC Number")
pwp_fig <- pwp_fig + geom_hline(yintercept = 0, colour = gray(1 / 2), lty = 2)
pwp_fig <- pwp_fig + geom_linerange(aes(
  x = Variable, ymin = Beta - SE * interval1,
  ymax = Beta + SE * interval1
),
lwd = 1.5, position = position_dodge(width = 1 / 2)
)
pwp_fig <- pwp_fig + geom_pointrange(aes(
  x = Variable, y = Beta, ymin = Beta - SE * interval2,
  ymax = Beta + SE * interval2, size = 2
),
lwd = 1, position = position_dodge(width = 1 / 2),
shape = 18
)
pwp_fig <- pwp_fig + xlim(rev(levels(allmod$Variable))) + coord_flip() + theme_minimal()
pwp_fig <- pwp_fig + labs(y = "Coefficient Estimate")
pwp_fig <- pwp_fig + theme(
  axis.title.y = element_blank(),
  axis.title.x = element_text(family = "serif", size = 14),
  axis.text.y = element_text(family = "serif", size = 14),
  legend.position = c(.9, .3),
  legend.title = element_text(family = "serif", size = 16),
  legend.text = element_text(family = "serif", size = 14)
)

# pwp_fig <- pwp_fig + ggtitle("PWP Model Results")
print(pwp_fig)

## Version without All TCs

allmod <- subset(allmod, Variable != "TC Tally")
allmod$Variable <- factor(allmod$Variable, levels = unique(as.character(allmod$Variable)))
allmod$TC <- factor(allmod$TC, levels = c("First TC", "2-3", "4-5", "6+"))
allmod$TC <- factor(allmod$TC, levels = c("6+", "4-5", "2-3", "First TC"))

pwp_fig2 <- ggplot(allmod, aes(colour = TC))
pwp_fig2 <- pwp_fig2 + scale_color_viridis_d(option = "C", breaks = c("First TC", "2-3", "4-5", "6+"), name = "TC Number")
pwp_fig2 <- pwp_fig2 + geom_hline(yintercept = 0, colour = gray(1 / 2), lty = 2)
pwp_fig2 <- pwp_fig2 + geom_linerange(aes(
  x = Variable, ymin = Beta - SE * interval1,
  ymax = Beta + SE * interval1
),
lwd = 1.5, position = position_dodge(width = 1 / 2)
)
pwp_fig2 <- pwp_fig2 + geom_pointrange(aes(
  x = Variable, y = Beta, ymin = Beta - SE * interval2,
  ymax = Beta + SE * interval2, size = 2
),
lwd = 1, position = position_dodge(width = 1 / 2),
shape = 18
)
pwp_fig2 <- pwp_fig2 + xlim(rev(levels(allmod$Variable))) + coord_flip() + theme_minimal()
pwp_fig2 <- pwp_fig2 + labs(y = "Coefficient Estimate")
pwp_fig2 <- pwp_fig2 + theme(
  axis.title.y = element_blank(),
  axis.title.x = element_text(family = "serif", size = 14),
  axis.text.y = element_text(family = "serif", size = 14),
  legend.position = c(.9, .3),
  legend.title = element_text(family = "serif", size = 16),
  legend.text = element_text(family = "serif", size = 14)
)

print(pwp_fig2)



#### Elapsed time version

coefs_all <- tidy(tc_pwp_elapse, conf.int = TRUE, exponentiate = F)
coefs_all$HazardRatio <- exp(coefs_all$estimate)
coefs_all$Variable <- c(
  "Peace w/ TC", "Fighting w/ TC", "Favorable Outcome for State", "Favorable Outcome for TC", "Absorbed",
  "Liberal Democracy Index", "State Area (logged)", "Mountains", "ELF", "TC Tally"
)

coefs1 <- tidy(tc1_pwp_elapse, conf.int = TRUE, exponentiate = F)
coefs1$HazardRatio <- exp(coefs1$estimate)
coefs1$Variable <- c("Liberal Democracy Index", "State Area (logged)", "Mountains", "ELF")

coefs23 <- tidy(tc23_pwp_elapse, conf.int = TRUE, exponentiate = F)
coefs23$HazardRatio <- exp(coefs23$estimate)
coefs23$Variable <- c(
  "Peace w/ TC", "Fighting w/ TC", "Favorable Outcome for State", "Favorable Outcome for TC", "Absorbed",
  "Liberal Democracy Index", "State Area (logged)", "Mountains", "ELF"
)

coefs4plus <- tidy(tc4_more_pwp_elapse, conf.int = TRUE, exponentiate = F)
coefs4plus$HazardRatio <- exp(coefs4plus$estimate)
coefs4plus$Variable <- c(
  "Peace w/ TC", "Fighting w/ TC", "Favorable Outcome for State", "Favorable Outcome for TC", "Absorbed",
  "Liberal Democracy Index", "State Area (logged)", "Mountains", "ELF"
)



combined <- data.frame(
  Variable = coefs_all$Variable,
  Beta = coefs_all$estimate,
  HazardRatio = coefs_all$HazardRatio,
  SE = summary(tc_pwp_elapse)$coef[, 4],
  TC = "All TCs"
)

first <- data.frame(
  Variable = coefs1$Variable,
  Beta = coefs1$estimate,
  HazardRatio = coefs1$HazardRatio,
  SE = summary(tc1_pwp_elapse)$coef[, 4],
  TC = "First TC"
)

mod23 <- data.frame(
  Variable = coefs23$Variable,
  HazardRatio = coefs23$HazardRatio,
  Beta = coefs23$estimate,
  SE = summary(tc23_pwp_elapse)$coef[, 4],
  TC = "2-3"
)

mod4plus <- data.frame(
  Variable = coefs4plus$Variable,
  HazardRatio = coefs4plus$HazardRatio,
  Beta = coefs4plus$estimate,
  SE = summary(tc4_more_pwp_elapse)$coef[, 4],
  TC = "4+"
)



allmod <- data.frame(rbind(combined, first, mod23, mod4plus))
# allmod <- subset(allmod, Variable!="TC Tally") #TC Tally Removed!!!!!!
allmod$Variable <- factor(allmod$Variable, levels = unique(as.character(allmod$Variable)))

allmod$TC <- factor(allmod$TC, levels = c("All TCs", "First TC", "2-3", "4+"))
allmod$TC <- factor(allmod$TC, levels = c("4+", "2-3", "First TC", "All TCs"))

elapse_fig <- ggplot(allmod, aes(colour = TC))
elapse_fig <- elapse_fig + scale_color_viridis_d(option = "C", breaks = c("All TCs", "First TC", "2-3", "4+"), name = "TC Number")
elapse_fig <- elapse_fig + geom_hline(yintercept = 0, colour = gray(1 / 2), lty = 2)
elapse_fig <- elapse_fig + geom_linerange(aes(
  x = Variable, ymin = Beta - SE * interval1,
  ymax = Beta + SE * interval1
),
lwd = 1.5, position = position_dodge(width = 1 / 2)
)
elapse_fig <- elapse_fig + geom_pointrange(aes(
  x = Variable, y = Beta, ymin = Beta - SE * interval2,
  ymax = Beta + SE * interval2, size = 2
),
lwd = 1, position = position_dodge(width = 1 / 2),
shape = 18
)
elapse_fig <- elapse_fig + xlim(rev(levels(allmod$Variable))) + coord_flip() + theme_minimal()
elapse_fig <- elapse_fig + labs(y = "Coefficient Estimate")
elapse_fig <- elapse_fig + theme(
  axis.title.y = element_blank(),
  axis.title.x = element_text(family = "serif", size = 14),
  axis.text.y = element_text(family = "serif", size = 14),
  legend.position = c(.8, .3),
  legend.title = element_text(family = "serif", size = 16),
  legend.text = element_text(family = "serif", size = 14)
)

print(elapse_fig)

