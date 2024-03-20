library(ompr)
library(gurobi)
library(tidyverse)
library(ompr.roi)
library(ROI)            # load package
##remotes::install_github("FlorianSchwendinger/ROI.plugin.gurobi", INSTALL_opts="--no-multiarch", force = T)
library(roi.plugin)
ROI_installed_solvers() 
ROI_available_solvers() # check available solvers
ROI_installed_solvers() # check installed solvers

set.seed(1234)
grid_size <- 10000
n <- 1000
customer_locations <- data.frame(
  id = 1:n,
  x = round(runif(n) * grid_size),
  y = round(runif(n) * grid_size)
)

m <- 200
warehouse_locations <- data.frame(
  id = 1:m,
  x = round(runif(m) * grid_size),
  y = round(runif(m) * grid_size)
)
fixedcost <- round(rnorm(m, mean = grid_size * 10, sd = grid_size * 5))

transportcost <- function(i, j) {
  customer <- customer_locations[i, ]
  warehouse <- warehouse_locations[j, ]
  round(sqrt((customer$x - warehouse$x)^2 + (customer$y - warehouse$y)^2))
}
transportcost(1, 3)

model <- MILPModel() %>%
  # 1 iff i gets assigned to warehouse j
  add_variable(x[i, j], i = 1:n, j = 1:m, type = "binary") %>%
  
  # 1 iff warehouse j is built
  add_variable(y[j], j = 1:m, type = "binary") %>%
  
  # maximize the preferences
  set_objective(sum_expr(transportcost(i, j) * x[i, j], i = 1:n, j = 1:m) + 
                  sum_expr(fixedcost[j] * y[j], j = 1:m), "min") %>%
  
  # every customer needs to be assigned to a warehouse
  add_constraint(sum_expr(x[i, j], j = 1:m) == 1, i = 1:n) %>% 
  
  # if a customer is assigned to a warehouse, then this warehouse must be built
  add_constraint(x[i,j] <= y[j], i = 1:n, j = 1:m)
model

result <- model %>% 
  solve_model(with_ROI(solver = "gurobi", verbose = TRUE))
result

matching <- result %>% 
  get_solution(x[i,j]) %>%
  filter(value > .9) %>%  
  select(i, j)

