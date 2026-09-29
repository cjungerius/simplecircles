# ==========================================
# 1. Base Conversion Helpers
# ==========================================

#' Conversion between degrees and radians
#'
#' @description
#' Helper functions to convert between degrees and radians.
#'
#' @param x A numeric vector of angles.
#'
#' @returns A numeric vector of converted values.
#'
#' @examples
#' deg2rad(47)
#' deg2rad(c(10, 30, 40))
#' rad2deg(pi / 7)
#'
#' @name degrad
NULL

#' @rdname degrad
#' @export
deg2rad <- function(x) { x * (pi / 180) }

#' @rdname degrad
#' @export
rad2deg <- function(x) { x * (180 / pi) }

# ==========================================
# 2. Circular Difference & Distance
# ==========================================

#' Shortest Directed Circular Difference
#'
#' Calculates the shortest element-wise directed difference between two angles.
#' Unlike simple algebraic subtraction, this function correctly accounts for
#' circular wrapping.
#'
#' @param x A numeric vector of angles.
#' @param y A numeric vector of reference angles to subtract from `x`.
#' @param units A character string specifying the unit of measurement.
#'   Must be one of `"degrees"` (default) or `"radians"`.
#'
#' @returns A numeric vector of shortest directed differences. The output range
#'   is centered around zero: `[-180, 180)` for degrees, or `[-pi, pi)` for radians.
#'   Positive values indicate that `x` is counter-clockwise from `y`.
#'
#' @examples
#' # Simple difference crossing the 360/0 boundary
#' circ_diff(x = 10, y = 350)
#'
#' # Vectorized operations
#' headings <- c(80, 10, 350)
#' target <- 45
#' circ_diff(headings, target)
#'
#' @export
circ_diff <- function(x, y, units = c("degrees", "radians")) {
  units <- match.arg(units)

  if (units == "degrees") {
    ((x - y + 180) %% 360) - 180
  } else {
    ((x - y + pi) %% (2 * pi)) - pi
  }
}

#' Absolute Circular Distance
#'
#' Calculates the absolute shortest distance between two angles.
#'
#' @param x A numeric vector of angles.
#' @param y A numeric vector of reference angles.
#' @param units A character string specifying the unit of measurement.
#'   Must be one of `"degrees"` (default) or `"radians"`.
#'
#' @returns A numeric vector of absolute distances `[0, 180]` or `[0, pi]`.
#'
#' @examples
#' circ_dist(10, 350) # Returns 20, not -20 or 340
#'
#' @export
circ_dist <- function(x, y, units = c("degrees", "radians")) {
  abs(circ_diff(x, y, units = units))
}

# ==========================================
# 3. Circular Wrapping
# ==========================================

#' Wrap Angles to standard circular boundaries
#'
#' Constrains angles into a standard single rotation.
#'
#' @param x A numeric vector of angles.
#' @param bounds A character string specifying the target range.
#'   `"positive"` constrains to `[0, 360)` or `[0, 2pi)`.
#'   `"centered"` constrains to `[-180, 180)` or `[-pi, pi)`.
#' @param units A character string specifying the unit of measurement.
#'   Must be one of `"degrees"` (default) or `"radians"`.
#'
#' @returns A numeric vector of normalized angles.
#'
#' @examples
#' circ_wrap(c(370, -45), bounds = "positive")
#' circ_wrap(c(370, 315), bounds = "centered")
#'
#' @export
circ_wrap <- function(x, bounds = c("positive", "centered"), units = c("degrees", "radians")) {
  bounds <- match.arg(bounds)
  units <- match.arg(units)

  if (bounds == "positive") {
    if (units == "degrees") x %% 360 else x %% (2 * pi)
  } else {
    circ_diff(x, 0, units = units)
  }
}


# ==========================================
# 4. Circular Central Tendency & Dispersion
# ==========================================

#' Mean Resultant Length
#'
#' Computes the mean resultant length (often denoted as *R*) of circular data.
#' This is a measure of concentration for circular distributions.
#'
#' @param x A numeric vector of angles.
#' @param na.rm Logical value indicating whether `NA` values should be stripped.
#' @param units A character string specifying the unit of measurement.
#'   Must be one of `"degrees"` (default) or `"radians"`.
#'
#' @returns A numeric value between 0 (uniformly distributed) and 1 (perfectly concentrated).
#'
#' @examples
#' circ_resultant_length(c(90, 90, 90)) # 1
#' circ_resultant_length(c(0, 90, 180, 270)) # 0
#'
#' @export
circ_resultant_length <- function(x, na.rm = FALSE, units = c("degrees", "radians")) {
  units <- match.arg(units)

  if (units == "degrees") {
    c <- mean(cospi(x / 180), na.rm = na.rm)
    s <- mean(sinpi(x / 180), na.rm = na.rm)
  } else {
    c <- mean(cos(x), na.rm = na.rm)
    s <- mean(sin(x), na.rm = na.rm)
  }

  sqrt(c^2 + s^2)
}

#' Circular Mean
#'
#' Computes the average angle of a vector of circular data.
#'
#' @param x A numeric vector of angles.
#' @param na.rm Logical value indicating whether `NA` values should be stripped.
#' @param units A character string specifying the unit of measurement.
#'   Must be one of `"degrees"` (default) or `"radians"`.
#'
#' @returns The mean angle in the specified units.
#'
#' @examples
#' circ_mean(c(350, 10))
#' circ_mean(c(0, 90, 180, NA), na.rm = TRUE)
#'
#' @export
circ_mean <- function(x, na.rm = FALSE, units = c("degrees", "radians")) {
  units <- match.arg(units)

  if (units == "degrees") {
    s <- mean(sinpi(x / 180), na.rm = na.rm)
    c <- mean(cospi(x / 180), na.rm = na.rm)
    atan2(s, c) * (180 / pi)
  } else {
    s <- mean(sin(x), na.rm = na.rm)
    c <- mean(cos(x), na.rm = na.rm)
    atan2(s, c)
  }
}

#' Circular Variance
#'
#' Computes the Mardia circular variance, calculated as 1 minus the mean
#' resultant length.
#'
#' @param x A numeric vector of angles.
#' @param na.rm Logical value indicating whether `NA` values should be stripped.
#' @param units A character string specifying the unit of measurement.
#'   Must be one of `"degrees"` (default) or `"radians"`.
#'
#' @returns A numeric value between 0 (no variance) and 1 (maximum variance).
#'
#' @examples
#' circ_var(c(5, 10, 15))
#' circ_var(c(0, 90, 180, 270))
#'
#' @export
circ_var <- function(x, na.rm = FALSE, units = c("degrees", "radians")) {
  1 - circ_resultant_length(x, na.rm = na.rm, units = units)
}

#' Circular Standard Deviation
#'
#' Computes the circular standard deviation of a numeric vector of angles.
#'
#' @param x A numeric vector of angles.
#' @param na.rm Logical value indicating whether `NA` values should be stripped.
#' @param units A character string specifying the unit of measurement.
#'   Must be one of `"degrees"` (default) or `"radians"`.
#'
#' @returns The circular standard deviation in the specified units. Returns
#'   `Inf` if the data is perfectly uniformly distributed (variance = 1).
#'
#' @examples
#' circ_sd(c(350, 5, 10))
#'
#' @export
circ_sd <- function(x, na.rm = FALSE, units = c("degrees", "radians")) {
  units <- match.arg(units)

  r <- circ_resultant_length(x, na.rm = na.rm, units = units)
  sd_rad <- sqrt(-2 * log(r))

  if (units == "degrees") {
    return(rad2deg(sd_rad))
  }

  sd_rad
}


# ==========================================
# 5. Circular Correlation & Scaling
# ==========================================

#' Circular Correlation
#'
#' Computes the Jammalamadaka circular correlation coefficient between two
#' vectors of angles.
#'
#' @param x A numeric vector of angles.
#' @param y A numeric vector of angles of the same length as `x`.
#' @param na.rm Logical value indicating whether pairs with `NA` values
#'   should be stripped.
#' @param units A character string specifying the unit of measurement.
#'   Must be one of `"degrees"` (default) or `"radians"`.
#'
#' @returns A numeric value between -1 and 1.
#'
#' @examples
#' wind_dir_day1 <- c(10, 45, 350, 20)
#' wind_dir_day2 <- c(15, 50, 355, 25)
#' circ_cor(wind_dir_day1, wind_dir_day2)
#'
#' @export
circ_cor <- function(x, y, na.rm = FALSE, units = c("degrees", "radians")) {
  units <- match.arg(units)

  if (na.rm) {
    valid <- complete.cases(x, y)
    x <- x[valid]
    y <- y[valid]
  }

  mean_x <- circ_mean(x, units = units)
  mean_y <- circ_mean(y, units = units)

  diff_x <- circ_diff(x, mean_x, units = units)
  diff_y <- circ_diff(y, mean_y, units = units)

  if (units == "degrees") {
    sin_x <- sinpi(diff_x / 180)
    sin_y <- sinpi(diff_y / 180)
  } else {
    sin_x <- sin(diff_x)
    sin_y <- sin(diff_y)
  }

  sum(sin_x * sin_y) / sqrt(sum(sin_x^2) * sum(sin_y^2))
}

#' Scale a Linear Variable to a Circular Range
#'
#' Maps a linear sequence bounded by a maximum value (e.g., 24-hour time or
#' 365-day years) onto a circular scale.
#'
#' @param x A numeric vector of linear values.
#' @param max_val The maximum value of the linear scale that represents a full cycle.
#' @param units A character string specifying the target unit of measurement.
#'   Must be one of `"degrees"` (default) or `"radians"`.
#'
#' @returns A numeric vector of angles.
#'
#' @examples
#' # Convert noon and midnight to degrees
#' scale_to_circ(c(12, 24), max_val = 24)
#'
#' # Convert day 182 of the year to radians
#' scale_to_circ(182.5, max_val = 365, units = "radians")
#'
#' @export
scale_to_circ <- function(x, max_val, units = c("degrees", "radians")) {
  units <- match.arg(units)

  proportion <- (x %% max_val) / max_val

  if (units == "degrees") {
    proportion * 360
  } else {
    proportion * (2 * pi)
  }
}
