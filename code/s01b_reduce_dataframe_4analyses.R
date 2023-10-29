# Make file smaller for downstream analyses

tidypath = here("data", "tidy")
tidy_data = read.csv(file.path(tidypath,"tidy_H_original_revalsubgroups_global_01.02.2023.csv"))

data2save = tidy_data %>%
  filter(dx1_short == "ASD" | dx1_short == "TD") %>%
  filter(task == "RestingState" & 
           exclude == "No" & 
           sex == "Male")
write.csv(data2save, file = file.path(tidypath,"tidy_H_revalsubgroups_global.csv"))