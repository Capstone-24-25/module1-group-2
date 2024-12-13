library(conflicted)
library(tidyverse)
library(tidymodels)
library(vip)
library(ranger)
tidymodels_prefer()
set.seed(1000)

load('data/biomarker-clean.RData')

biomarker = biomarker_clean %>% 
  select(-c(ados)) %>%
  mutate(group = as.factor(group))
set.seed(1000) #splitting data and folds
biomarker_split = initial_split(biomarker, prop = 0.8, strata = group)
biomarker_train = training(biomarker_split)
biomarker_test = training(biomarker_split)
biomarker_fold = vfold_cv(biomarker_train, v = 10, strata = group)

#recipe for variable selection
biomarker_recipe = recipe(group ~ ., data = biomarker_train) %>%
  step_normalize(all_predictors())

prep(biomarker_recipe) %>%
  bake(new_data = biomarker_train)

#setting tuning variables
random_forest = rand_forest(mtry = tune(),
                            trees = tune(),
                            min_n = tune()) %>%
  set_engine('ranger', importance = 'impurity') %>%
  set_mode('classification')

elastic_net = logistic_reg(penalty = tune(), 
                         mixture = tune()) %>% 
  set_mode("classification") %>% 
  set_engine("glmnet")

#workflow for models
rf_wf = workflow() %>%
  add_model(random_forest) %>%
  add_recipe(biomarker_recipe)

elastic_wf = workflow() %>%
  add_model(elastic_net) %>%
  add_recipe(biomarker_recipe)

#setting range of tune values
rf_grid = grid_regular(mtry(range = c(1, 5)), 
                       trees(range = c(100,500)), 
                       min_n(range = c(5,20)), 
                       levels = 10)

elastic_grid = grid_regular(penalty(), 
                            mixture(range = c(0,1)), 
                            levels = 10)

#fitting and tuning models
#rf_tune = tune_grid(
#  rf_wf,
#  resamples = biomarker_fold,
#  grid = rf_grid)
#save(rf_tune, file = 'data/rf_tune.rda')

#elastic_tune = tune_grid(
#  elastic_wf,
#  resamples = biomarker_fold,
#  grid = elastic_grid)
#save(elastic_tune, file = 'data/elastic_tune.rda')

load('data/elastic_tune.rda')
load('data/rf_tune.rda')

#finding best model by accuracy
set.seed(1000)
best_rf_train = select_best(rf_tune, metric = 'accuracy')
rf_final_wf_train = finalize_workflow(rf_wf, best_rf_train)
rf_final_fit = fit(rf_final_wf_train, data = biomarker_train)

best_elastic_train = select_best(elastic_tune, metric = 'accuracy')
elastic_final_wf_train = finalize_workflow(elastic_wf, best_elastic_train)
elastic_final_fit = fit(elastic_final_wf_train, data = biomarker_train)

#looking at variables that matter
elastic_estimate = elastic_final_fit %>% 
  tidy() %>%
  filter(estimate != 0)

rf_final_fit %>%
  extract_fit_engine() %>%
  vip(num_feature = 50, aesthetics = list(fill = 'red4'))

#putting important variables into data frame
rf_var = data.frame(term = names(rf_final_fit$fit$fit$fit$variable.importance), 
                    imp = rf_final_fit$fit$fit$fit$variable.importance, 
                    row.names = NULL) %>%
  arrange(-imp) %>%
  select(-imp) %>%
  slice(1:10)

elastic_var = elastic_estimate %>%
  select(term)

#finding intersection between the two models
var_select = merge(rf_var, elastic_var) %>%
  list()

#new recipe with common important variables only
biomarker_small_recipe = recipe(group ~ `Calcineurin` + `Cystatin C` + 
                                  FSTL1 + IgD + 
                                  MAPK14 + MAPK2 + 
                                  PPID + RELT, 
                                data = biomarker_train) %>%
  step_normalize(all_predictors())

prep(biomarker_small_recipe) %>%
  bake(new_data = biomarker_train)

#setting model type
log_reg = logistic_reg() %>%
  set_engine('glm') %>%
  set_mode('classification')

#setting up workflow
log_wf = workflow() %>%
  add_model(log_reg) %>%
  add_recipe(biomarker_small_recipe)

#fitting logistic model
log_fit = fit(log_wf, biomarker_train)

#looking at estimates
log_estimate = log_fit %>%
  tidy()

#extracting accuracy
acc = augment(log_fit, biomarker_test) %>% accuracy(group, .pred_class)

# The model achieves an accuracy of 0.779 with 8 predictor proteins 
# List of predictor proteins is stored as `var_select`