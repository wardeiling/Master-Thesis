# This script generates a figure that shows how marginal and conditional effects may not 
# be equivalent under certain conditions (GM3) in a multilevel linear model.

### Example without Problem ----

### Step 1. Obtain data for 10 different individuals with increasing random effects

# create a vector with different values for the random intercept, mean = 0, sigma_u0 = 1
sample10 <- data.frame(ri_vals = seq(from = -2.5, to = 2.5, by = 0.4),
                             intercept = 2,
                             slope = 0.8)
sample10$intercept_i <- sample10$intercept - sample10$ri_vals

par(mfrow = c(1, 2))

### Step 2. Obtain the conditional average: the subject with b_i = 0, the mean individual

# when we fill in sigma_u = 0, the random effect is 0 for all individuals.
# The estimated effect should therefore represent the conditional average.

cond_slope <- 0.8
cond_intercept <- 2

### Step 3. Obtain the population average: the average Y for every level of X over all individuals

marg_slope <- 0.8
marg_intercept <- 2.10

### Step 4. Create Figure

# X is from -5 to 15
# Y is from -5 to 15

y_range <- x_range <- seq(-5, 10, length.out = 100)

plot(x_range, y_range, type = "n", xlab = "X", ylab = "Y")

for (i in 1:nrow(sample10)) {
  abline(a = sample10$intercept_i[i], b = 0.8, col = "grey", lwd = 2)
}

abline(a = cond_intercept, b = cond_slope, col = "black", lwd = 4)
abline(a = marg_intercept, b = marg_slope, col = "red", , lwd = 4)

### Example with Problem ----

### Step 2. Obtain the conditional average: the subject with b_i = 0, the mean individual

# when we fill in sigma_u = 0, the random effect is 0 for all individuals.
# The estimated effect should therefore represent the conditional average.

cond_slope <- 0.8
cond_intercept <- 2.0

### Step 3. Obtain the population average: the average Y for every level of X over all individuals

marg_slope <- 0.98
marg_intercept <- 1.82

### Step 4. Create Figure

# X is from -5 to 15
# Y is from -5 to 15

y_range <- x_range <- seq(-5, 10, length.out = 100)

plot(x_range, y_range, type = "n", xlab = "X", ylab = "Y")

for (i in 1:nrow(sample10)) {
  abline(a = sample10$intercept_i[i], b = 0.8, col = "grey", lwd = 2)
}

abline(a = cond_intercept, b = cond_slope, col = "black", lwd = 4)
abline(a = marg_intercept, b = marg_slope, col = "red", , lwd = 4)

