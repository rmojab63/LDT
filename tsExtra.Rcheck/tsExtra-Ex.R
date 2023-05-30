pkgname <- "tsExtra"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "tsExtra-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('tsExtra')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("BindVariables")
### * BindVariables

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: BindVariables
### Title: Binds a of Variables
### Aliases: BindVariables

### ** Examples

v1 = ldt::Variable(c(1,2,3,2,3,4,5),"V1",F_Monthly(2022,12), list())
v2 = ldt::Variable(c(10,20,30,20,30,40,50),"V2",F_Monthly(2022,8), list())
L = ldt::BindVariables(list(v1,v2))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("BindVariables", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("F_CrossSection")
### * F_CrossSection

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: F_CrossSection
### Title: Creates a Cross-Section Frequency
### Aliases: F_CrossSection

### ** Examples


cs0 <- F_CrossSection(10) # this initializes a cross-section frequency

cs0_value_str <-  as.character(cs0) # this will be '10'.
cs0_class_str <- get.class.id(cs0) # this will be 'cs'.

cs_new <- as.ldtf("20", "cs") # this is a cross-section frequency. It points to position 20.




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("F_CrossSection", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("F_Daily")
### * F_Daily

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: F_Daily
### Title: Creates a Daily Frequency
### Aliases: F_Daily

### ** Examples


d0 <- F_Daily(2023, 1, 2) # This is 2/1/2023. Next observation belongs to 3/1/2023.

d0_value_str <-  as.character(d0) # this will be '20230102'.
d0_class_str <- get.class.id(d0) # this will be 'd'.

d_new <- as.ldtf("20230109", "d") # This is 9/1/2023.

## Not run: 
##D # Don't use invalid or unsupported dates.
##D 
##D d_invalid <- as.ldtf("1399109", "d") # this is a too old date and unsupported
##D d_invalid <- as.ldtf("20230132", "d") # invalid day in month
##D d_invalid <- as.ldtf("20231331", "d") # invalid month
##D 
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("F_Daily", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("F_DailyInWeek")
### * F_DailyInWeek

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: F_DailyInWeek
### Title: Creates an 'Daily-In-Week' Frequency
### Aliases: F_DailyInWeek

### ** Examples


dw0 <- F_DailyInWeek(2023, 5, 16, "mon", "fri") # This is 16/5/2023.
dw0_value_str <-  as.character(dw0) # this will be '20230516'.
dw0_class_str <- get.class.id(dw0) # this will be 'i:mon-fri'.

# Let's use the same date with another week definition:
dw1 <- F_DailyInWeek(2023, 5, 16, "wed", "sat") # This is NOT 16/5/2023. It is 17/5/2023. Since it was outside the week, we moved it forward.
dw2 <- F_DailyInWeek(2023, 5, 16, "wed", "sat", FALSE) # This is 13/5/2023. The original day was outside the week, but we moved backward too the end of the previous week (which is Saturday).

dw_new <- as.ldtf("20230519", "i:sat-wed") # This is 20/1/2023 (by default, it moves forward).

## Not run: 
##D # Don't use invalid or unsupported dates.
##D 
##D dw_invalid <- as.ldtf("1399109", "d3") # this is a too old date and unsupported
##D dw_invalid <- as.ldtf("20230132", "d4") # invalid day in month
##D dw_invalid <- as.ldtf("20231331", "d5") # invalid month
##D 
##D # don't use invalid week definitions:
##D F_DailyInWeek(2023, 5, 16, "Wednesday", "sat")
##D 
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("F_DailyInWeek", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("F_Hourly")
### * F_Hourly

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: F_Hourly
### Title: Creates an 'Hourly' Frequency
### Aliases: F_Hourly

### ** Examples


ho0 <- F_Hourly(F_Daily(2023,5,16),4)

ho0_value_str <-  as.character(ho0) # this will be '20230516:4'.
ho0_class_str <- get.class.id(ho0) # this will be 'ho|d'. The second part (i.e., 'd') shows that this frequency is defined in a 'Daily' frequency.

ho_new <- as.ldtf("20231101:3", "ho|i:wed-sat")

## Not run: 
##D # Don't make the following mistakes:
##D 
##D ho_invalid <- as.ldtf("20231101:3", "ho|j:wed-sat") # invalid format in day-based frequency
##D ho_invalid <- F_Hourly(F_Daily(2023,5,16),25) # invalid hour
##D 
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("F_Hourly", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("F_ListDate")
### * F_ListDate

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: F_ListDate
### Title: Creates an 'List-Date' Frequency
### Aliases: F_ListDate

### ** Examples


Ld0 <- F_ListDate(c("20231101","20220903","20200823","20230303"), "20200823")

Ld0_value_str <-  as.character(Ld0) # this will be '20200823'.
Ld0_class_str <- get.class.id(Ld0) # this will be 'Ld:20231101;20220903;20200823;20230303'.

Ld_new <- as.ldtf("20231101", "Ld:20231101;20220903;20200823;20230303")
Ld_new0 <- as.ldtf("20231101", "Ld") # compared to the previous one, its items will be empty

## Not run: 
##D # Don't make the following mistakes:
##D 
##D Ld_invalid <- as.ldtf("20231102", "Ld:20231101;20220903;20200823;20230303") # 'E' is not a member of the list
##D Ld_invalid <- F_ListDate(c("20231101","20220903","20200823","20230303"), "20231102")
##D 
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("F_ListDate", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("F_ListString")
### * F_ListString

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: F_ListString
### Title: Creates an 'List-String' Frequency
### Aliases: F_ListString

### ** Examples


L0 <- F_ListString(c("A","B","C","D"), "C")

L0_value_str <-  as.character(L0) # this will be 'C'.
L0_class_str <- get.class.id(L0) # this will be 'Ls:A;B;C;D'.

L_new <- as.ldtf("A", "Ls:A;B;C;D")
L_new0 <- as.ldtf("A", "Ls") # compared to the previous one, its items will be empty

## Not run: 
##D # Don't make the following mistakes:
##D 
##D L_invalid <- as.ldtf("E", "Ls:A;B;C;D") # 'E' is not a member of the list
##D L_invalid <- F_ListString(c("A","B","C","D"), "E")
##D 
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("F_ListString", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("F_Minutely")
### * F_Minutely

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: F_Minutely
### Title: Creates an 'Minute-ly' Frequency
### Aliases: F_Minutely

### ** Examples


mi0 <- F_Minutely(F_Daily(2023,5,16),1200)

mi0_value_str <-  as.character(mi0) # this will be '20230516:1200'.
mi0_class_str <- get.class.id(mi0) # this will be 'mi|d'. The second part (i.e., 'd') shows that this frequency is defined in a 'Daily' frequency.

mi_new <- as.ldtf("20231101:3", "mi|i:wed-sat")

## Not run: 
##D # Don't make the following mistakes:
##D 
##D mi_invalid <- as.ldtf("20231101:3", "mi|j:wed-sat") # invalid format in day-based frequency
##D mi_invalid <- F_Minutely(F_Daily(2023,5,16),2000) # invalid minute
##D 
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("F_Minutely", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("F_Monthly")
### * F_Monthly

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: F_Monthly
### Title: Creates a Monthly Frequency
### Aliases: F_Monthly

### ** Examples


m0 <- F_Monthly(2020, 2) # this is a monthly frequency that refers to the second month of the year 2020.

m0_value_str <-  as.character(m0) # this will be '2020M2'.
m0_class_str <- get.class.id(m0) # this will be 'm'.

m_new <- as.ldtf("2021m3", "m") # this is a monthly frequency that refers to the third month of the year 2021.

## Not run: 
##D # Don't make the following mistakes:
##D 
##D m_invalid <- F_Monthly(2020, 0)
##D m_invalid <- F_Monthly(2020, 5)
##D m_invalid <- as.ldtf("2021m0", "m")
##D m_invalid <- as.ldtf("2021m13", "m")
##D m_invalid <- as.ldtf("2021", "m")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("F_Monthly", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("F_MultiDaily")
### * F_MultiDaily

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: F_MultiDaily
### Title: Creates an Multi-Day Frequency
### Aliases: F_MultiDaily

### ** Examples


md0 <- F_MultiDaily(2023, 1, 2, 4) # This is 2/1/2023. Next observation belongs to 6/1/2023.

md0_value_str <-  as.character(md0) # this will be '20230102'.
md0_class_str <- get.class.id(md0) # this will be 'd4'.

md_new <- as.ldtf("20230109", "d") # This is 9/1/2023.

## Not run: 
##D # Don't use invalid or unsupported dates.
##D 
##D md_invalid <- as.ldtf("1399109", "d3") # this is a too old date and unsupported
##D md_invalid <- as.ldtf("20230132", "d4") # invalid day in month
##D md_invalid <- as.ldtf("20231331", "d5") # invalid month
##D 
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("F_MultiDaily", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("F_MultiWeekly")
### * F_MultiWeekly

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: F_MultiWeekly
### Title: Creates a Multi-Week Frequency
### Aliases: F_MultiWeekly

### ** Examples


mw0 <- F_MultiWeekly(2023, 1, 2,3) # This is 2/1/2023 which is Monday. Next observation belongs to 23/1/2023.

mw0_value_str <-  as.character(mw0) # this will be '20230102'.
mw0_class_str <- get.class.id(mw0) # this will be 'w3'.

mw_new <- as.ldtf("20230109", "w4") # This is 9/1/2023.

## Not run: 
##D # Don't use invalid or unsupported dates.
##D 
##D mw_invalid <- as.ldtf("1399109", "w4") # this is a too old date and unsupported
##D mw_invalid <- as.ldtf("20230132", "w5") # invalid day in month
##D mw_invalid <- as.ldtf("20231331", "w2") # invalid month
##D mw_invalid <- as.ldtf("20231012", "w0")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("F_MultiWeekly", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("F_MultiYearly")
### * F_MultiYearly

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: F_MultiYearly
### Title: Creates a Multi-Year Frequency
### Aliases: F_MultiYearly

### ** Examples


my0 <- F_MultiYearly(2020, 2) # this is a multi-year frequency that refers to the year 2020. The next observation is expected in 2022 (not the next year).

my0_value_str <-  as.character(my0) # this will be '2020'.
my0_class_str <- get.class.id(my0) # this will be 'z2'.

my_new <- as.ldtf("2020", "z3") # this is a multi-year frequency that refers to the year 2020. However, the next observation is expected in 2023.

## Not run: 
##D # Don't make the following mistakes:
##D 
##D my_invalid <- F_MultiYearly(2020, 0)
##D my_invalid <- F_MultiYearly(2020, -5)
##D my_invalid <- as.ldtf("2021", "z")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("F_MultiYearly", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("F_Quarterly")
### * F_Quarterly

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: F_Quarterly
### Title: Creates a Quarterly Frequency
### Aliases: F_Quarterly

### ** Examples


q0 <- F_Quarterly(2020, 2) # this is a quarterly frequency that refers to the second quarter of the year 2021.

q0_value_str <-  as.character(q0) # this will be '2020Q2'.
q0_class_str <- get.class.id(q0) # this will be 'q'.

q_new <- as.ldtf("2021q3", "q") # this is a quarterly frequency that refers to the third quarter of the year 2021.

## Not run: 
##D # Don't make the following mistakes:
##D 
##D q_invalid <- F_Quarterly(2020, 0)
##D q_invalid <- F_Quarterly(2020, 5)
##D q_invalid <- as.ldtf("2021q0", "q")
##D q_invalid <- as.ldtf("2021q5", "q")
##D q_invalid <- as.ldtf("2021", "q")
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("F_Quarterly", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("F_Secondly")
### * F_Secondly

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: F_Secondly
### Title: Creates an 'Second-ly' Frequency
### Aliases: F_Secondly

### ** Examples


se0 <- F_Secondly(F_Daily(2023,5,16),40032)

se0_value_str <-  as.character(se0) # this will be '20230516:40032'.
se0_class_str <- get.class.id(se0) # this will be 'se|d'. The second part (i.e., 'd') shows that this frequency is defined in a 'Daily' frequency.

se_new <- as.ldtf("20231101:3", "se|i:wed-sat")

## Not run: 
##D # Don't make the following mistakes:
##D 
##D mi_invalid <- as.ldtf("20231101:3", "se|j:wed-sat") # invalid format in day-based frequency
##D mi_invalid <- F_Secondly(F_Daily(2023,5,16),100000) # invalid second
##D 
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("F_Secondly", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("F_Weekly")
### * F_Weekly

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: F_Weekly
### Title: Creates a Weekly Frequency
### Aliases: F_Weekly

### ** Examples


w0 <- F_Weekly(2023, 1, 2) # This is 2/1/2023 which is Monday. Next observation belongs to 9/1/2023.

w0_value_str <-  as.character(w0) # this will be '20230102'.
w0_class_str <- get.class.id(w0) # this will be 'w'.

w_new <- as.ldtf("20230109", "w") # This is 9/1/2023.

## Not run: 
##D # Don't use invalid or unsupported dates:
##D 
##D w_invalid <- as.ldtf("1399109", "w") # this is a too old date and unsupported
##D w_invalid <- as.ldtf("20230132", "w") # invalid day in month
##D w_invalid <- as.ldtf("20231331", "w") # invalid month
##D 
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("F_Weekly", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("F_XTimesADay")
### * F_XTimesADay

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: F_XTimesADay
### Title: Creates an 'X-Times-A-Day' Frequency
### Aliases: F_XTimesADay

### ** Examples


xd0 <- F_XTimesADay(F_Daily(2023,5,16),13, 12)

xd0_value_str <-  as.character(xd0) # this will be '20230516:12'.
xd0_class_str <- get.class.id(xd0) # this will be 'da13|d'. The second part (i.e., 'd') shows that this frequency is defined in a 'Daily' frequency.

xd_new <- as.ldtf("20231101:3", "da3|i:wed-sat")

## Not run: 
##D # Don't make the following mistakes:
##D 
##D xd_invalid <- as.ldtf("20231101:3", "da|i:wed-sat") # invalid format in day-based frequency
##D xd_invalid <- F_XTimesADay(F_Daily(2023,5,16),4,0) # invalid position
##D 
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("F_XTimesADay", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("F_XTimesAYear")
### * F_XTimesAYear

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: F_XTimesAYear
### Title: Creates an 'X-Times-A-Year' Frequency
### Aliases: F_XTimesAYear

### ** Examples


xty0 <- F_XTimesAYear(2020, 3, 1) # this frequency divides the year 2020 into 3 partitions and refers to the first partition.

xty_value_str <-  as.character(xty0) # this will be '2020:1'.
xty_class_str <- get.class.id(xty0) # this will be 'y3'.

xty_new <- as.ldtf("2021:24", "z24") # this frequency divides the year 2021 into 24 partitions and refers to the last partition.

## Not run: 
##D # Don't make the following mistakes:
##D 
##D xty_invalid <- F_XTimesAYear(2020, 3, 0)
##D xty_invalid <- F_XTimesAYear(2020, 24, 25)
##D xty_invalid <- as.ldtf("2021:13", "y12")
##D xty_invalid <- as.ldtf("2021:0", "y1")
##D xty_invalid <- as.ldtf("2021", "y1")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("F_XTimesAYear", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("F_XTimesZYears")
### * F_XTimesZYears

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: F_XTimesZYears
### Title: Creates an 'X-Times-Z-Years' Frequency
### Aliases: F_XTimesZYears

### ** Examples


xtzy0 <- F_XTimesZYears(2020, 3, 2, 3) # this frequency divides the year 2020 into 3 partitions and refers to the last partition. The next observation belongs to 2022 (not the next year).

xtzy_value_str <-  as.character(xtzy0) # this will be '2020:3'.
xtzy_class_str <- get.class.id(xtzy0) # this will be 'x3z2'.

xtzy_new <- as.ldtf("2021:3", "x3z4") # this frequency divides the year 2021 into 3 partitions and refers to the last partition. The next observation occurs after 4 years.

## Not run: 
##D # Don't make the following mistakes:
##D 
##D xtzy_invalid <- F_XTimesZYears(2020, 3, 5, 0)
##D xtzy_invalid <- F_XTimesZYears(2020, 3, 0, 1)
##D xtzy_invalid <- as.ldtf("2021:25", "x24y2")
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("F_XTimesZYears", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("F_Yearly")
### * F_Yearly

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: F_Yearly
### Title: Creates an Annual Frequency
### Aliases: F_Yearly

### ** Examples


y0 <- F_Yearly(2020) # this initializes a 'yearly' frequency

y0_value_str <-  as.character(y0) # this will be '2020'.
y0_class_str <- get.class.id(y0) # this will be 'y'.

y_new <- as.ldtf("2021", "y") # this is a yearly frequency. It points to year 2021.




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("F_Yearly", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Variable")
### * Variable

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Variable
### Title: Creates a Variable
### Aliases: Variable

### ** Examples

v1 = ldt::Variable(c(1,2,3,2,3,4,5),"V1",F_Monthly(2022,12),
     list(c("key1","value1"), c("key2", "value2")))



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Variable", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("get.character.info")
### * get.character.info

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get.character.info
### Title: Return Value and Class as String
### Aliases: get.character.info

### ** Examples


freq <- F_XTimesADay(F_Daily(2023,5,16),13, 12)
freq_class_id <- get.character.info(freq) # this will be 'da13|d'.




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get.character.info", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("get.class.id")
### * get.class.id

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: get.class.id
### Title: Gets Class 'Id' of a Frequency
### Aliases: get.class.id

### ** Examples


freq <- F_XTimesADay(F_Daily(2023,5,16),13, 12)
freq_class_id <- get.class.id(freq) # this will be 'da13|d'.




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("get.class.id", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("seq")
### * seq

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: seq
### Title: Generates a Sequence from a Range of Frequency
### Aliases: seq

### ** Examples

from <- F_Monthly(2020,1)
to <- F_Monthly(2021,12)
sequence1 <- seq(from, to, 1) # this will be '2020M1', '2020M2', ..., '2021M12'
sequence2 <- seq(from, to, 2) # this will be '2020M1', '2020M3', ..., '2021M11'
sequence3 <- seq(from, to, 3) # this will be '2020M1', '2020M4', ..., '2021M10'

# backward:
sequence4 <- seq(to, from, -1) # this will be '2021M12', '2021M11', ..., '2020M1'




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("seq", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("seq0")
### * seq0

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: seq0
### Title: Generates a Sequence from a Range of Frequency
### Aliases: seq0

### ** Examples

start <- F_Monthly(2020,1)
sequence1 <- seq0(start, 24, 1) # this will be '2020M1', '2020M2', ..., '2021M12'
sequence2 <- seq0(start, 24, 2) # this will be '2020M1', '2020M3', ..., '2023M11'
sequence3 <- seq0(start, 24, 3) # this will be '2020M1', '2020M4', ..., '2025M10'

# backward:
sequence4 <- seq0(start, 24, -1) # this will be '2020M1', '2019M12', ..., '2018M2'




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("seq0", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
