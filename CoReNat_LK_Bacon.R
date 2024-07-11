# CoReNat Age-depth modelling

# 29.06.2023

# load library

library(rbacon) # for c14 cal + age-depth model
library(rplum)
library(rintcal) # getting c14 from F14C dates
library(tidyverse) # data wrangling

# using Martin Blauuw tutorials: https://maarten14c.github.io/

# First get post-bomb c14 ages

F14C.age(1.01168, sdev = c(0.004657), decimals = 0) # LKN-98

## calibrate 14C with rintcal

## LKN 

LKN192 <- caldist(age= 145, error = 19, "intcal20")

stat192 <- point.estimates(LKN192) %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  rename("value" = ".") %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate("Depth" = 192)

LKN192 <- LKN192 %>%
  as.data.frame() %>%
  rename(cal_BP = "cal BP") %>%
  summarise("Core" = "LKN", "Depth" = 192, min = min(cal_BP), max = max(cal_BP))

LKN98 <- caldist(age= -90, error = 40, postbomb = TRUE,
               "nh3", BCAD=FALSE)

stat98 <- point.estimates(LKN98) %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  rename("value" = ".") %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate("Depth" = 98)

LKN98 <- LKN98 %>%
  as.data.frame() %>%
  rename(cal_BP = "cal BP") %>%
  #names the other way around due to - numbers for min and max
  summarise("Core" = "LKN", "Depth" = 98, max = min(cal_BP), min = max(cal_BP))

LKN <- LKN192 %>%
  full_join(LKN98)

LKNstat <- stat192 %>%
  full_join(stat98)

LKN <- LKN %>%
  left_join(LKNstat, by = "Depth")


####### Modelling Likupang sites with pb-210 CRS dates and 14C#######

# following https://github.com/maquinolopez/BES-Palaeo-Sig---Plum

Bacon("LKN", postbomb = 3, d.max = 193)

agedepth()

accrate.depth.ghost()

accrate.age.ghost()

d0 <- mean(accrate.depth(0))
d2 <- mean(accrate.depth(2))
d5 <- mean(accrate.depth(5))
d10 <- mean(accrate.depth(10))
d15 <- mean(accrate.depth(15))
d20 <- mean(accrate.depth(20))

d25 <- mean(accrate.depth(25))
d36 <- mean(accrate.depth(36))
d40 <- mean(accrate.depth(40))
d45<- mean(accrate.depth(45))
d55<- mean(accrate.depth(55))
d59<- mean(accrate.depth(59))

d64<- mean(accrate.depth(64))
d74 <- mean(accrate.depth(74))
d85 <- mean(accrate.depth(85))
d89<- mean(accrate.depth(89))

d95<- mean(accrate.depth(95))

d106<- mean(accrate.depth(106))
d116<- mean(accrate.depth(116))
d121<- mean(accrate.depth(121))

d134<- mean(accrate.depth(134))
d143<- mean(accrate.depth(143))
d156<- mean(accrate.depth(156))

d164<- mean(accrate.depth(164))
d177<- mean(accrate.depth(177))
d181<- mean(accrate.depth(181))

d185<- mean(accrate.depth(185))
d188<- mean(accrate.depth(188))
d193<- mean(accrate.depth(193))

mean_acc <- vctrs::vec_c(d0,d2, d5,d10,d15,d20, d25, d36, d40,d45,d55,d59,d64,
        d74,d85,d89,d95, d106,d116,d121, d134, d143, d156, d164, d177, d181, d185,
        d188, d193)

age <- read.csv("LKN_crs_age.csv")%>%
  mutate(Depth = depth) %>%
  mutate(Age_BP = mean) %>%
  select(Depth, Age_BP)

dep <- read.csv("LKN_depth.csv") 
  
comb <- dep %>%
  left_join(age)

LKN_avg_sed_acc <- comb %>%
  arrange(Depth)%>%
  mutate(Mean_sed = mean_acc) 

#write.csv(LKN_avg_sed_acc, file = "LKN_avg_sed_acc.csv")

#LKRE

LKRE300 <- caldist(age= 708, error = 17, "intcal20")

stat300 <- point.estimates(LKRE300) %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  rename("value" = ".") %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate("Depth" = 300)

LKRE300 <- LKRE300 %>%
  as.data.frame() %>%
  rename(cal_BP = "cal BP") %>%
  summarise("Core" = "LKRE", "Depth" = 300, min = min(cal_BP), max = max(cal_BP))

LKRE220 <- caldist(age= 577, error = 37, "intcal20")

stat220 <- point.estimates(LKRE220) %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  rename("value" = ".") %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate("Depth" = 220)

LKRE220 <- LKRE220 %>%
  as.data.frame() %>%
  rename(cal_BP = "cal BP") %>%
  summarise("Core" = "LKRE", "Depth" = 220, min = min(cal_BP), max = max(cal_BP))

LKRE150 <- caldist(age= 260, error = 35, "intcal20")

stat150 <- point.estimates(LKRE150) %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  rename("value" = ".") %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate("Depth" = 150)

LKRE150 <- LKRE150 %>%
  as.data.frame() %>%
  rename(cal_BP = "cal BP") %>%
  summarise("Core" = "LKRE", "Depth" = 150, min = min(cal_BP), max = max(cal_BP))

LKRE100 <- caldist(age= 40, error = 30, "intcal20")

stat100 <- point.estimates(LKRE100) %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  rename("value" = ".") %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate("Depth" = 100)

LKRE100 <- LKRE100 %>%
  as.data.frame() %>%
  rename(cal_BP = "cal BP") %>%
  summarise("Core" = "LKRE", "Depth" = 100, min = min(cal_BP), max = max(cal_BP))

LKRE <- LKRE300 %>%
  full_join(LKRE220) %>%
  full_join(LKRE150) %>%
  full_join(LKRE100)

LKREstat <- stat300 %>%
  full_join(stat220) %>%
  full_join(stat150) %>%
  full_join(stat100)

LKRE <- LKRE %>%
  left_join(LKREstat, by = "Depth")

Bacon("LKRE")

agedepth()

d101<- mean(accrate.depth(101))
d125<- mean(accrate.depth(125))
d135<- mean(accrate.depth(135))
d149<- mean(accrate.depth(149))
d170<- mean(accrate.depth(170))
d180<- mean(accrate.depth(180))
d190<- mean(accrate.depth(190))
d200<- mean(accrate.depth(200))
d210<- mean(accrate.depth(210))
d220<- mean(accrate.depth(220))
d230<- mean(accrate.depth(230))
d270<- mean(accrate.depth(270))
d274<- mean(accrate.depth(274))
d280<- mean(accrate.depth(280))
d290<- mean(accrate.depth(290))
d295<- mean(accrate.depth(295))
d300<- mean(accrate.depth(300))


mean_acc <- vctrs::vec_c(d101,d125, d135, d149, d170, d180, d190,d200,d210,
                         d220,d230, d270,d274, d280,d290,d295, d300)

LKRE_age <- read.csv("LKRE_ages.csv")

LKRE_avg_sed_acc <- read.csv("LKRE_depth.csv") %>%
  left_join(LKRE_age)%>%
  arrange(Depth)%>%
  mutate(Mean_sed = mean_acc)%>%
  select(Depth, Age_BP, Mean_sed)

write.csv(LKRE_avg_sed_acc, file = "LKRE_avg_sed_acc.csv")


