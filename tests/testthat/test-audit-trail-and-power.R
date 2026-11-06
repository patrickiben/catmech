test_that("catmech_audit_trail generates auditable reports in markdown, json, and html", {
  d <- simulate_catmech_lab(120, seed = 42)
  d <- d[setdiff(names(d), c("lab_true", "corrupted"))]
  fit <- catmech_model(d, "lab_value", c("age", "sex", "outcome"), flag = "outlier_flag", m = 3, seed = 12)

  # Pooled model
  f <- outcome ~ lab_value + age
  pooled <- catmech_pool_lm(fit, f)

  # Markdown format
  trail_md <- catmech_audit_trail(fit, pooled = pooled, format = "markdown")
  expect_s3_class(trail_md, "catmech_audit_trail")
  txt_md <- as.character(trail_md)
  expect_true(grepl("# Catmech Audit Trail & Compliance Report", txt_md))
  expect_true(grepl("Original Input Data Hash:", txt_md))
  expect_true(grepl("Completion `m=1`:", txt_md))
  expect_true(grepl("Downstream Pooled Inference", txt_md))
  expect_true(grepl("lab_value", txt_md))

  # Test print method
  pr_txt <- utils::capture.output(val <- print(trail_md))
  expect_identical(val, trail_md)
  expect_true(any(grepl("Regulatory Audit & Compliance Certification", pr_txt)))

  # JSON format
  trail_json <- catmech_audit_trail(fit, pooled = pooled, format = "json")
  expect_s3_class(trail_json, "catmech_audit_trail")
  parsed <- jsonlite::fromJSON(as.character(trail_json))
  expect_identical(parsed$specification$engine, "gaussian_model")
  expect_identical(parsed$specification$completions_m, 3L)
  expect_true(nzchar(parsed$integrity$original_data_sha256))
  expect_equal(length(parsed$integrity$completed_datasets_sha256), 3L)

  # HTML format
  trail_html <- catmech_audit_trail(fit, pooled = pooled, format = "html")
  expect_s3_class(trail_html, "catmech_audit_trail")
  txt_html <- as.character(trail_html)
  expect_true(grepl("<!DOCTYPE html>", txt_html))
  expect_true(grepl("<title>Catmech Audit Trail & Compliance Report</title>", txt_html))

  # File writing
  tmp_file <- tempfile(fileext = ".md")
  on.exit(unlink(tmp_file), add = TRUE)
  catmech_audit_trail(fit, pooled = pooled, format = "markdown", path = tmp_file)
  expect_true(file.exists(tmp_file))
  read_back <- paste(readLines(tmp_file), collapse = "\n")
  expect_true(grepl("Original Input Data Hash:", read_back))

  # Sensitivity object inclusion
  sens <- catmech_sensitivity(list(default = fit), formula = f, term = "lab_value")
  trail_sens <- catmech_audit_trail(fit, pooled = pooled, sensitivity = sens, format = "markdown")
  expect_true(grepl("## 4. Sensitivity Envelope", as.character(trail_sens)))
})

test_that("catmech_audit_power computes prospective and plan-based power and plots curve", {
  # Prospective calculation
  pwr_prosp <- catmech_audit_power(
    N = 1000,
    target_ratio = 1.25,
    power_target = 0.80,
    error_rate = 0.15,
    correlation = 0.85,
    fractions = c(0.05, 0.10, 0.15, 0.20, 0.25)
  )

  expect_s3_class(pwr_prosp, "catmech_audit_power")
  expect_identical(pwr_prosp$N, 1000L)
  expect_identical(nrow(pwr_prosp$table), 5L)
  expect_true(all(pwr_prosp$table$power >= 0 & pwr_prosp$table$power <= 1))
  expect_true(all(pwr_prosp$table$expected_ratio >= 1))
  # Monotonicity: power should generally non-decrease with larger audit fraction
  expect_true(pwr_prosp$table$power[5L] >= pwr_prosp$table$power[1L])

  # Test print method
  pwr_txt <- utils::capture.output(val_pwr <- print(pwr_prosp))
  expect_identical(val_pwr, pwr_prosp)
  expect_true(any(grepl("Catmech Audit Power & Sample Size Calculator", pwr_txt)))
  expect_true(any(grepl("Recommended Audit Size:", pwr_txt)))

  # Test plot method
  tf <- tempfile(fileext = ".png")
  on.exit(unlink(tf), add = TRUE)
  grDevices::png(tf, width = 600, height = 400)
  plot(pwr_prosp)
  grDevices::dev.off()
  expect_true(file.exists(tf) && file.size(tf) > 0)

  # Calculation from an existing audit plan object
  # Mock a minimal audit plan
  plan_mock <- structure(
    list(
      k = data.frame(strategy = "recorded", estimate = 0.40, lower = 0.20, upper = 0.70, stringsAsFactors = FALSE),
      pilot = list(fit = list(original = data.frame(x = 1:500)))
    ),
    class = "catmech_audit_plan"
  )

  pwr_plan <- catmech_audit_power(plan = plan_mock, target_ratio = 1.25, power_target = 0.80)
  expect_s3_class(pwr_plan, "catmech_audit_power")
  expect_identical(pwr_plan$N, 500L)
  expect_identical(pwr_plan$k_summary$source, "empirical_plan")

  # Error handling
  expect_error(catmech_audit_power(target_ratio = 0.9), "target_ratio must be a single number strictly greater than 1")
  expect_error(catmech_audit_power(power_target = 1.5), "power_target must be strictly between 0 and 1")
  expect_error(catmech_audit_power(N = NULL, plan = NULL), "N must be a single integer >= 10")
})
