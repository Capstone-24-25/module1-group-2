library(glmnet)
library(tidyverse)
library(tidymodels)
library(modelr)
library(yardstick)
load("data/biomarker-clean.RData")

# Split testing and training data and formating data
set.seed(101422)

biomarker_split <- biomarker_clean %>%
  initial_split(prop = 0.8)

training_biomarker <- training(biomarker_split)
testing_biomarker <- testing(biomarker_split)

x <- training_biomarker %>%
  select(-ados, -group) %>%
  as.matrix()
y <- training_biomarker %>%
  select(group) %>%
  mutate(group = (group == "ASD")) %>%
  pull(group)

# Lasso regression model
model <- glmnet(x, y, family = "binomial", alpha = 1)
cv.model <- cv.glmnet(x, y, family = "binomial", alpha = 1)

best_lambda <- cv.model$lambda.min

final_model <- glmnet(x, y, family = "binomial", lambda = best_lambda)

coefficients <- coef(final_model, s = best_lambda)
selected_variables <- rownames(coefficients)[which(coefficients != 0)]

# Evaluating on test data

x_test <- testing_biomarker %>%
  select(-ados, -group) %>%
  as.matrix()

predicted_prob <- predict(final_model, s = best_lambda, newx = x_test, type = "response")

predicted_group <- ifelse(predicted_prob > 0.5, "ASD", "TD")

testing_biomarker <- testing_biomarker %>%
  mutate(predicted_group = predicted_group)

confusion_matrix <- table(testing_biomarker$group, testing_biomarker$predicted_group)
