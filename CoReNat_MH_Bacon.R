# CoReNat Age-depth modelling

# 22.11.2022

# load library

library(rbacon) # for c14 cal + age-depth model
library(rplum) # for Pb210 cal + age-depth model
library(rintcal) # getting c14 from F14C dates
library(tidyverse) # data wrangling
# using Martin Blauuw tutorials: https://maarten14c.github.io/

# First get post-bomb c14 ages

F14C.age(1.28394, sdev = c(0.005909), decimals = 0) # MH20

F14C.age(10.075, sdev = c(0.0038), decimals = 0) # MH30

## calibrate 14C with rintcal

MH600 <- caldist(age= 7053, error = 35, "intcal20")

stat600 <- point.estimates(MH600) %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  rename("value" = ".") %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate("Depth" = 600)

MH600 <- MH600 %>%
  as.data.frame() %>%
  rename(cal_BP = "cal BP") %>%
  summarise("Core" = "MH", "Depth" = 600, min = min(cal_BP), max = max(cal_BP))
  
MH540 <- caldist(age= 6950, error = 30, "intcal20")

stat540 <- point.estimates(MH540) %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  rename("value" = ".") %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate("Depth" = 540)

MH540 <- MH540 %>%
  as.data.frame() %>%
  rename(cal_BP = "cal BP") %>%
  summarise("Core" = "MH", "Depth" = 540, min = min(cal_BP), max = max(cal_BP))



MH500 <- caldist(age= 5938, error = 35, "intcal20")

stat500 <- point.estimates(MH500) %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  rename("value" = ".") %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate("Depth" = 500)

MH500 <- MH500 %>%
  as.data.frame() %>%
  rename(cal_BP = "cal BP") %>%
  summarise("Core" = "MH", "Depth" = 500, min = min(cal_BP), max = max(cal_BP))

MH410 <- caldist(age= 6528, error = 39, "intcal20")

stat410 <- point.estimates(MH410) %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  rename("value" = ".") %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate("Depth" = 410)

MH410 <- MH410 %>%
  as.data.frame() %>%
  rename(cal_BP = "cal BP") %>%
  summarise("Core" = "MH", "Depth" = 410, min = min(cal_BP), max = max(cal_BP))

MH360 <- caldist(age= 5968, error = 23, "intcal20")

stat360 <- point.estimates(MH360) %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  rename("value" = ".") %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate("Depth" = 360)

MH360 <- MH360 %>%
  as.data.frame() %>%
  rename(cal_BP = "cal BP") %>%
  summarise("Core" = "MH", "Depth" = 360, min = min(cal_BP), max = max(cal_BP))

MH300 <- caldist(age= 6092, error = 37, "intcal20")

stat300 <- point.estimates(MH300) %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  rename("value" = ".") %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate("Depth" = 300)

MH300 <- MH300 %>%
  as.data.frame() %>%
  rename(cal_BP = "cal BP") %>%
  summarise("Core" = "MH", "Depth" = 300, min = min(cal_BP), max = max(cal_BP))


MH220 <- caldist(age= 5951, error = 36, "intcal20")

stat220 <- point.estimates(MH220) %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  rename("value" = ".") %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate("Depth" = 220)

MH220 <- MH220 %>%
  as.data.frame() %>%
  rename(cal_BP = "cal BP") %>%
  summarise("Core" = "MH", "Depth" = 220, min = min(cal_BP), max = max(cal_BP))


MH190 <- caldist(age= 4820, error = 30, "intcal20")

stat190 <- point.estimates(MH190) %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  rename("value" = ".") %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate("Depth" = 190)

MH190 <- MH190 %>%
  as.data.frame() %>%
  rename(cal_BP = "cal BP") %>%
  summarise("Core" = "MH", "Depth" = 190, min = min(cal_BP), max = max(cal_BP))

  
MH150 <- caldist(age= 3450, error = 38, "intcal20")

stat150 <- point.estimates(MH150) %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  rename("value" = ".") %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate("Depth" = 150)

MH150 <- MH150 %>%
  as.data.frame() %>%
  rename(cal_BP = "cal BP") %>%
  summarise("Core" = "MH", "Depth" = 150, min = min(cal_BP), max = max(cal_BP))

  
MH80 <- caldist(age= 2595, error = 37, "intcal20")

stat80 <- point.estimates(MH80) %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  rename("value" = ".") %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate("Depth" = 80)

MH80 <- MH80 %>%
  as.data.frame() %>%
  rename(cal_BP = "cal BP") %>%
  summarise("Core" = "MH", "Depth" = 80, min = min(cal_BP), max = max(cal_BP))

MH66 <- caldist(age= 1610, error = 30, "intcal20")

stat66 <- point.estimates(MH66) %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  rename("value" = ".") %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate("Depth" = 66)

MH66 <- MH66 %>%
  as.data.frame() %>%
  rename(cal_BP = "cal BP") %>%
  summarise("Core" = "MH", "Depth" = 66, min = min(cal_BP), max = max(cal_BP))

MH60 <- caldist(age= 521, error = 35, "intcal20")

stat60 <- point.estimates(MH60) %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  rename("value" = ".") %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate("Depth" = 60)

MH60 <- MH60 %>%
  as.data.frame() %>%
  rename(cal_BP = "cal BP") %>%
  summarise("Core" = "MH", "Depth" = 60, min = min(cal_BP), max = max(cal_BP))

MH56 <- caldist(age= 420, error = 30, "intcal20")

stat56 <- point.estimates(MH56) %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  rename("value" = ".") %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate("Depth" = 56)

MH56 <- MH56 %>%
  as.data.frame() %>%
  rename(cal_BP = "cal BP") %>%
  summarise("Core" = "MH", "Depth" = 56, min = min(cal_BP), max = max(cal_BP))

MH30 <- caldist(age= -2000, error = 3, postbomb = TRUE,
                "nh3", BCAD=FALSE)

stat30 <- point.estimates(MH30) %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  rename("value" = ".") %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate("Depth" = 30)

MH30 <- MH30 %>%
  as.data.frame() %>%
  rename(cal_BP = "cal BP") %>%
  summarise("Core" = "MH", "Depth" = 30, min = min(cal_BP), max = max(cal_BP))


MH20 <- caldist(age= -2000, error = 40, postbomb = TRUE,
                "nh3", BCAD=FALSE)

stat20 <- point.estimates(MH20) %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  rename("value" = ".") %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate("Depth" = 20)

MH20 <- MH20 %>%
  as.data.frame() %>%
  rename(cal_BP = "cal BP") %>%
  #names the other way around due to - numbers for min and max
  summarise("Core" = "MH", "Depth" = 20, max = min(cal_BP), min = max(cal_BP))


MH <- MH600 %>%
  full_join(MH540) %>%
  full_join(MH500) %>%
  full_join(MH410) %>%
  full_join(MH360) %>%
  full_join(MH300) %>%
  full_join(MH220) %>%
  full_join(MH190) %>%
  full_join(MH150) %>%
  full_join(MH80) %>%
  full_join(MH66) %>%
  full_join(MH60) %>%
  full_join(MH56) %>%
  full_join(MH30) %>%
  full_join(MH20) 

MHstat <- stat600 %>%
  full_join(stat540) %>%
  full_join(stat500) %>%
  full_join(stat410) %>%
  full_join(stat360) %>%
  full_join(stat300) %>%
  full_join(stat220) %>%
  full_join(stat190) %>%
  full_join(stat150) %>%
  full_join(stat80) %>%
  full_join(stat66) %>%
  full_join(stat60) %>%
  full_join(stat56) %>%
  full_join(stat30) %>%
  full_join(stat20) 

MH <- MH %>%
  left_join(MHstat, by = "Depth")

####### Modelling MH with rbacon ##########

# MH from a mangrove lagoon on Mantehage Island, 600 cm core.
# Based on pollen, 600 - 500 mangroves expanding, 
# accretion due to transgression should be increasing between start and end 
# 500 - 100 covers the mid-Holocene high-stand - likely sedimentation higher 
# during this point in time - also more bioturbation
# 60 - 64 significant lithological, palynological, diatom change due to :
        #1 chosen volcanic eruption - slump
      
# 60 - 0 mangrove recovering following event/period, 

# post 30 modern sediments

# Use this information to set priors

Bacon("MH", postbomb = 3, slump = c(60,63), 
      boundary = c(150, 220, 540), acc.mean=c(10, 20, 10, 1, 5))

agedepth()


d0 <- mean(accrate.depth(0))
d3 <- mean(accrate.depth(3))
d10 <- mean(accrate.depth(10))
d12 <- mean(accrate.depth(12))
d14 <- mean(accrate.depth(14))
d20 <- mean(accrate.depth(20))
d25 <- mean(accrate.depth(25))
d30 <- mean(accrate.depth(30))
d43 <- mean(accrate.depth(43))
d47 <- mean(accrate.depth(47))
d50 <- mean(accrate.depth(50))
d57 <- mean(accrate.depth(57))
d60 <- mean(accrate.depth(60))
d66 <- mean(accrate.depth(66))
d68 <- mean(accrate.depth(68))
d70 <- mean(accrate.depth(70))
d75 <- mean(accrate.depth(75))
d80 <- mean(accrate.depth(80))
d90 <- mean(accrate.depth(90))
d100<- mean(accrate.depth(100))
d130<- mean(accrate.depth(130))
d140<- mean(accrate.depth(140))
d150<- mean(accrate.depth(150))
d170<- mean(accrate.depth(170))
d190<- mean(accrate.depth(190))
d199<- mean(accrate.depth(199))
d210<- mean(accrate.depth(210))
d220<- mean(accrate.depth(220))
d230<- mean(accrate.depth(230))
d240<- mean(accrate.depth(240))
d260<- mean(accrate.depth(260))
d270<- mean(accrate.depth(270))
d280<- mean(accrate.depth(280))
d300<- mean(accrate.depth(300))
d320<- mean(accrate.depth(320))
d330<- mean(accrate.depth(330))
d340<- mean(accrate.depth(340))
d350<- mean(accrate.depth(350))
d360<- mean(accrate.depth(360))
d390<- mean(accrate.depth(390))
d410<- mean(accrate.depth(410))
d430<- mean(accrate.depth(430))
d449<- mean(accrate.depth(449))
d470<- mean(accrate.depth(470))
d480<- mean(accrate.depth(480))
d500<- mean(accrate.depth(500))
d520<- mean(accrate.depth(520))
d540<- mean(accrate.depth(530))
d570<- mean(accrate.depth(570))
d580<- mean(accrate.depth(580))
d590<- mean(accrate.depth(590))
d600<- mean(accrate.depth(600))



mean_acc <- vctrs::vec_c(d0,d3,d10,d12,d14,d20,d25,d30,d43,d47,d50,d57,d60,
                         d66,d68,d70,d75,d80,d90,d100,d130,d140,d150,d170,d190,
                         d199,d210,d220,d230,d240,d260,d270,d280,d300,d320,
                         d330,d340,d350,d360,d390,d410,d430,d449,d470,d480,d500,
                         d520,d540,d570,d580,d590,d600)


depth <- read.csv("Pol_depth.csv")

chron <- read.csv("MH_ages_slump_selected.csv")%>%
  mutate(Depth = round(depth, digits = 0))%>%
  mutate(Age_BP = mean) %>%
  select(Depth, Age_BP)

comb <- depth %>%
  left_join(chron)%>%
  group_by(Depth)%>%
  slice_head()%>%
  ungroup()

MH_avg_sed_acc <- comb %>%
  arrange(Depth)%>%
  mutate(Mean_sed = mean_acc) 

#write.csv(MH_avg_sed_acc, file = "MH_avg_sed_acc_slump.csv")

plot(MH_avg_sed_acc$Mean_sed)

accrate.depth.ghost()

accrate.age.ghost()

agedepth(rotate.axes = T, rev.age = T) # plots agedepth models using last call

ages_bacon <- as_tibble(info$ranges)
ages_bacon %>% head() # gets the ages in table form for use in plotting diags

calib_bacon <- map2_dfr(
  info$calib$d, 
  info$calib$probs, 
  function(depth, prob_matrix) {
    colnames(prob_matrix) <- c("cal_age_bp", "density")
    df <- as_tibble(prob_matrix)
    df$depth <- depth
    df
  }) %>%
  select(depth, cal_age_bp, density)

calib_bacon %>% head() # all calibrated ages for each depth (large)

#write.csv(calib_bacon, "calib_bacon_slump.csv")
