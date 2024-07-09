## -------------------------------------------------------------------------------------------------------------
#| echo: false
#| warning: false
#| message: false

library(gt)
library(gtsummary)
library(tidyverse)
library(tidymodels)
library(tidytext)
#library(textrecipes)
library(here)
library(aplore3)


mean_age = mean(glow500$age) %>% round()
glow1 = glow500 %>% mutate(age_c = age - mean_age)

set.seed(513613)


## -------------------------------------------------------------------------------------------------------------
# model
model = logistic_reg()
# recipe
recipe = recipe(fracture ~ priorfrac + age_c, data = glow1) %>%
            step_dummy(priorfrac) %>%
            step_interact(terms = ~ age_c:starts_with("priorfrac"))
# workflow
workflow = workflow() %>% add_model(model) %>% add_recipe(recipe)

fit = workflow %>% fit(data = glow1)


## -------------------------------------------------------------------------------------------------------------
#| echo: false
#| code-fold: true

tidy(fit, conf.int = T) %>% gt() %>% 
  tab_options(table.font.size = 35) %>%
  fmt_number(decimals = 3)


## -------------------------------------------------------------------------------------------------------------
glow_m3 = glm(fracture ~ priorfrac + age_c + priorfrac*age_c, 
              data = glow1, family = binomial)


## -------------------------------------------------------------------------------------------------------------
tidy(glow_m3, conf.int = T) %>% gt() %>% 
  tab_options(table.font.size = 35) %>%
  fmt_number(decimals = 3)


## -------------------------------------------------------------------------------------------------------------
#| fig-width: 3
#| fig-height: 3

ggplot(glow1, aes(x = fracture)) + geom_bar()


## -------------------------------------------------------------------------------------------------------------
glow = glow1 %>%
    dplyr::select(-sub_id, -site_id, -phy_id, -age, -bmi, -weight)


## -------------------------------------------------------------------------------------------------------------
glow_split = initial_split(glow, strata = fracture, prop = 0.8)
glow_split


## -------------------------------------------------------------------------------------------------------------
glow_train = training(glow_split)
glow_test = testing(glow_split)


## -------------------------------------------------------------------------------------------------------------
glimpse(glow_train)


## -------------------------------------------------------------------------------------------------------------
glimpse(glow_test)


## -------------------------------------------------------------------------------------------------------------
lasso_mod = logistic_reg(penalty = 0.001, mixture = 1) %>%

            set_engine("glmnet")


## -------------------------------------------------------------------------------------------------------------
glow_rec_main = recipe(fracture ~ ., data = glow_train) %>%

  step_dummy(priorfrac, premeno, momfrac, armassist, smoke, raterisk)

glow_workflow_main = workflow() %>%

      add_model(lasso_mod) %>% add_recipe(glow_rec_main)
  
glow_fit_main = glow_workflow_main %>% fit(glow_train)

tidy.glmnet(glow_fit_main, conf.int = T) %>% gt() %>% 
  tab_options(table.font.size = 35) %>%
  fmt_number(decimals = 3)
## -------------------------------------------------------------------------------------------------------------
library(vip)  
vi_data_main = glow_fit_main %>% 
    pull_workflow_fit() %>%
    vi(lambda = 0.001) %>%
    filter(Importance != 0)
vi_data_main





####


log_model = logistic_reg()

glow_workflow_int2 = workflow() %>%
      add_model(log_model) %>% add_recipe(glow_rec_int2)
  
glow_fit_int2 = glow_workflow_int2 %>% fit(glow_train)


## -------------------------------------------------------------------------------------------------------------
print(tidy(glow_fit_int2), n=60)


## -------------------------------------------------------------------------------------------------------------
glow_test_pred = predict(glow_fit_int2, new_data = glow_test, type = "prob") %>%
    bind_cols(glow_test)


## -------------------------------------------------------------------------------------------------------------
glow_test_pred %>% 
  roc_auc(truth = fracture, 
                  .pred_No)


## -------------------------------------------------------------------------------------------------------------
#| fig-width: 7
#| fig-height: 7
#| fig-align: center

glow_test_pred %>% 
  roc_curve(truth = fracture, .pred_No) %>%
  autoplot()


## -------------------------------------------------------------------------------------------------------------
glow_test_pred = predict(glow_fit_int2, new_data = glow_test, type = "prob") %>%
    bind_cols(glow_test)


## -------------------------------------------------------------------------------------------------------------
glow_test_pred %>% 
  roc_auc(truth = fracture, 
                  .pred_No)


## -------------------------------------------------------------------------------------------------------------
#| fig-width: 7
#| fig-height: 7
#| fig-align: center

glow_test_pred %>% 
  roc_curve(truth = fracture, .pred_No) %>%
  autoplot()

