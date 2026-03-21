---
title: "Area Quantile EDA"
author: "Codex"
date: "2026-03-21"
output: html_document
---



## Purpose

This notebook explores the area-augmented `stan_ready_data_with_area` object. The goals are:

1. visualize size-distribution trajectories alongside counts,
2. quantify within-line ploidy contrasts over time,
3. test whether early morphology shifts relate to later failure or survival,
4. relate morphology to glucose depletion,
5. port the size-to-cost scaling idea from `workflow/volume_aware_anal.Rmd` onto the row-aligned area dataset.

## Setup



Loaded area-aware Stan data from dev_data/stan_ready_data_with_area.Rds.

The object contains 76 conditions, 5258 count observations, and 11 aligned area quantiles per count row.

## Coverage Check


|metric                          | value|
|:-------------------------------|-----:|
|count rows                      |  5258|
|rows with any area summary      |  5257|
|rows missing all area quantiles |     1|

The area join is effectively complete: only one of the 5,258 count observations lacks an aligned area summary, so downstream patterns should not be driven by missingness.

## Size Trajectories

These plots treat area as an additional phenotype channel parallel to live/dead counts. `q50` captures the middle of the alive-cell area distribution, while `q99` and `tail_ratio = q99 / q50` track the upper tail.

![plot of chunk size-trajectories](figure/size-trajectories-1.png)


|cellLine   |ploidy_label |    G0| peak_live| peak_q50| peak_tail_ratio|
|:----------|:------------|-----:|---------:|--------:|---------------:|
|MCF10A     |baseline     |  0.00|    1356.0|   759.75|            3.14|
|MCF10A     |elevated     |  0.00|    1967.5|  1361.50|            5.42|
|MCF10A     |baseline     |  0.10|    4099.5|   873.25|            6.96|
|MCF10A     |elevated     |  0.10|    2617.5|  1584.50|            6.31|
|MCF10A     |baseline     |  0.25|    4936.5|   852.50|            5.60|
|MCF10A     |elevated     |  0.25|    3108.5|  1436.75|            6.22|
|MCF10A     |baseline     |  0.50|    6694.5|   826.00|            4.37|
|MCF10A     |elevated     |  0.50|    4294.5|  1404.50|            6.62|
|MCF10A     |baseline     |  1.00|    8085.5|   810.00|            4.23|
|MCF10A     |elevated     |  1.00|    5312.0|  1356.25|            6.33|
|MCF10A     |baseline     |  5.00|   10304.0|   810.50|            5.46|
|MCF10A     |elevated     |  5.00|   11792.5|  1318.75|            6.07|
|MCF10A     |baseline     | 25.00|   11392.5|   796.00|            3.91|
|MCF10A     |elevated     | 25.00|   13018.5|  1369.00|            5.61|
|MDA-MB-231 |baseline     |  0.00|     382.0|   893.50|           16.29|
|MDA-MB-231 |elevated     |  0.00|     252.5|   989.50|           14.10|
|MDA-MB-231 |baseline     |  0.10|     814.0|  1565.50|           21.54|
|MDA-MB-231 |elevated     |  0.10|     681.0|  1828.00|           11.27|

The count and area trajectories already suggest that morphology is not just a proxy for cell number. In several conditions, `q50` and `tail_ratio` keep moving after the main expansion phase has slowed, which supports looking at morphology as a state variable rather than just a confluence readout.

## Ploidy Contrasts Through Time

These contrasts show the elevated-ploidy condition minus the baseline-ploidy condition within each cell line, glucose condition, and time point.

![plot of chunk ploidy-contrast-trajectories](figure/ploidy-contrast-trajectories-1.png)


|cellLine     |    G0| max_abs_delta_live| max_abs_delta_q50| max_abs_delta_tail_ratio|
|:------------|-----:|------------------:|-----------------:|------------------------:|
|MCF10A       |  0.10|             1941.5|           1069.25|                     2.63|
|MCF10A       |  0.25|             2254.0|            847.00|                     2.44|
|SUM-159-fuse |  1.00|             1273.0|            846.25|                     2.64|
|SNU668       |  5.00|             3518.5|            836.50|                     2.30|
|MCF10A       |  0.00|             1266.0|            815.25|                     3.05|
|SNU668       | 25.00|             3028.0|            791.00|                     2.14|
|SNU668       |  0.10|              187.0|            784.75|                     2.70|
|SUM-159-fuse |  0.50|              634.5|            760.75|                     1.08|
|SNU668       |  0.25|              259.0|            752.25|                     3.23|
|SNU668       |  1.00|              731.0|            750.50|                     4.75|
|SUM-159-fuse | 25.00|             4518.5|            739.25|                     0.80|
|SNU668       |  0.50|              379.0|            719.50|                     2.92|
|SUM-159-chem |  0.10|              318.5|            713.25|                     1.02|
|SUM-159-fuse |  5.00|             3916.5|            679.50|                     1.55|
|SUM-159-chem |  0.25|              276.5|            673.25|                     0.83|

The contrast view is useful because it makes it obvious where ploidy effects are primarily morphological versus primarily count-based. The larger `|delta_q50|` and `|delta_tail_ratio|` cases are the ones most worth following up in a mechanistic branch.

## Area-Rate and Death-Rate Coupling

This replaces the earlier early-versus-late scatter with a rate-based view. For each replicate, a smoothing spline is fit to dead-cell counts and to area `q50`. We then examine:

- `death_rate_per_cell = max(dD/dt, 0) / N`
- `area_rate = d(q50)/dt`

and compute lagged correlations between the two signals.

![plot of chunk lag-coupling-plots](figure/lag-coupling-plots-1.png)


|cellLine     | median_best_corr| median_best_lag_hours|
|:------------|----------------:|---------------------:|
|MCF10A       |            0.012|                     4|
|MDA-MB-231   |           -0.772|                     0|
|SNU668       |           -0.911|                     8|
|SUM-159-chem |           -0.901|                    24|
|SUM-159-fuse |           -0.786|                     8|

The lag summaries suggest that for most lines the area-change and death-rate signals are coupled with negative correlation and modest positive lag, meaning shrinking median area often coincides with or slightly precedes rising death pressure. `MCF10A` is the main exception here.

## Morphology and Glucose Depletion

These plots address confluence by dropping count observations above 70% of the maximum total-cell count for each `cellLine x ploidy` combination. They then compare median alive-cell area to current glucose, glucose lagged by 24h, and glucose lagged by 48h. Glucose values are shown on a log scale and clipped to `0.02` only for plotting.

Within each panel, points connected by a line are the median time trajectory for one `cellLine x ploidy x G0` condition after the confluence filter.

![plot of chunk morphology-glucose](figure/morphology-glucose-1.png)


|cellLine     | cor_q50_vs_glucose| cor_q50_vs_glucose_lag24| cor_q50_vs_glucose_lag48|
|:------------|------------------:|------------------------:|------------------------:|
|MCF10A       |             -0.111|                   -0.293|                   -0.380|
|MDA-MB-231   |              0.158|                    0.447|                    0.498|
|SNU668       |             -0.109|                    0.018|                    0.158|
|SUM-159-chem |             -0.094|                   -0.043|                   -0.263|
|SUM-159-fuse |              0.004|                    0.066|                   -0.043|

The useful summary here is whether area tracks current glucose or past glucose more strongly. If the 24h- or 48h-lag correlation is stronger, that supports a delayed morphology response to nutrient history rather than an instantaneous glucose effect.

In the current summaries that delayed-history pattern is most obvious for `MCF10A` and `MDA-MB-231`, where lagged glucose correlates with median area more strongly than contemporaneous glucose. The other lines look weaker or less consistent, so this should be treated as line-specific rather than universal.

## Size-To-Cost Scaling

This ports the most interesting idea from `workflow/volume_aware_anal.Rmd`: relate glucose cost per unit cell production to a size proxy. Here the calculation is rebuilt from the row-aligned area dataset and the cost model is explicitly fit as:

`glucose consumed ~ a * biomass built + m * live-cell AUC`

with `a >= 0` and `m >= 0`.

Here `glucose consumed` is not the nominal `G0`. It is the estimated glucose at the first glucose measurement minus the estimated glucose at the last glucose measurement for each well.


|cellLine     | ploidy_metric| ploidy_abs|ploidy_label |     a|  m|  rmse|  n| area_24h| vol_flat| vol_sphere| a_base| vol_flat_base| vol_sphere_base| ratio_cost| ratio_vol_flat| ratio_vol_sphere| cost_per_flat_vol| cost_per_sphere_vol|
|:------------|-------------:|----------:|:------------|-----:|--:|-----:|--:|--------:|--------:|----------:|------:|-------------:|---------------:|----------:|--------------:|----------------:|-----------------:|-------------------:|
|MCF10A       |          0.00|        2.0|baseline     | 0.000|  0| 0.714| 24|   813.00|   813.00|   23181.25|  0.000|        813.00|        23181.25|      1.000|          1.000|            1.000|                 0|                   0|
|MCF10A       |          1.00|        4.0|elevated     | 0.000|  0| 0.491| 24|  1274.50|  1274.50|   45499.83|  0.000|        813.00|        23181.25|      1.223|          1.568|            1.963|                 0|                   0|
|MDA-MB-231   |          0.00|        3.5|baseline     | 0.001|  0| 0.135| 24|  1523.00|  1523.00|   59436.05|  0.001|       1523.00|        59436.05|      1.000|          1.000|            1.000|                 0|                   0|
|MDA-MB-231   |          0.33|        4.4|elevated     | 0.001|  0| 0.092| 24|  1531.50|  1531.50|   59934.60|  0.001|       1523.00|        59436.05|      0.884|          1.006|            1.008|                 0|                   0|
|SNU668       |          0.00|        2.6|baseline     | 0.001|  0| 0.322| 24|   703.00|   703.00|   18639.67|  0.001|        703.00|        18639.67|      1.000|          1.000|            1.000|                 0|                   0|
|SNU668       |          1.00|        5.2|elevated     | 0.004|  0| 0.672| 24|  1422.75|  1422.75|   53665.22|  0.001|        703.00|        18639.67|      3.036|          2.024|            2.879|                 0|                   0|
|SUM-159-chem |          0.00|        2.0|baseline     | 0.001|  0| 0.263| 24|  1070.75|  1070.75|   35037.54|  0.001|       1070.75|        35037.54|      1.000|          1.000|            1.000|                 0|                   0|
|SUM-159-chem |          1.00|        4.0|elevated     | 0.002|  0| 0.160| 24|  1537.50|  1537.50|   60287.92|  0.001|       1070.75|        35037.54|      1.341|          1.436|            1.721|                 0|                   0|
|SUM-159-fuse |          0.00|        2.0|baseline     | 0.001|  0| 0.351| 32|  1680.75|  1680.75|   68906.04|  0.001|       1680.75|        68906.04|      1.000|          1.000|            1.000|                 0|                   0|
|SUM-159-fuse |          1.00|        4.0|elevated     | 0.001|  0| 0.303| 32|  1500.00|  1500.00|   58094.79|  0.001|       1680.75|        68906.04|      0.547|          0.892|            0.843|                 0|                   0|

The very small `m` estimates and fairly tight `p_build` relationships imply the current accounting fit is behaving almost like a one-parameter biomass-yield model, especially outside `MCF10A`. That is exactly why an implied biomass-trajectory view is worth plotting.

![plot of chunk cost-input-diagnostics](figure/cost-input-diagnostics-1.png)

![plot of chunk cost-scaling-plots](figure/cost-scaling-plots-1.png)

![plot of chunk cost-fit-quality](figure/cost-fit-quality-1.png)

![plot of chunk cost-implied-timecourse](figure/cost-implied-timecourse-1.png)![plot of chunk cost-implied-timecourse](figure/cost-implied-timecourse-2.png)![plot of chunk cost-implied-timecourse](figure/cost-implied-timecourse-3.png)![plot of chunk cost-implied-timecourse](figure/cost-implied-timecourse-4.png)![plot of chunk cost-implied-timecourse](figure/cost-implied-timecourse-5.png)


|cellLine   |ploidy_label |    G0| rmse_timecourse| cor_timecourse|
|:----------|:------------|-----:|---------------:|--------------:|
|MCF10A     |baseline     |  0.00|         547.999|         -0.962|
|MCF10A     |elevated     |  0.00|         403.032|             NA|
|MCF10A     |baseline     |  0.10|        1394.941|          0.892|
|MCF10A     |elevated     |  0.10|         261.675|          0.398|
|MCF10A     |baseline     |  0.25|        1616.337|          0.886|
|MCF10A     |elevated     |  0.25|         485.591|          0.881|
|MCF10A     |baseline     |  0.50|        2168.486|          0.992|
|MCF10A     |elevated     |  0.50|         525.570|          0.974|
|MCF10A     |baseline     |  1.00|        2372.712|          0.928|
|MCF10A     |elevated     |  1.00|         782.245|          0.949|
|MCF10A     |baseline     |  5.00|        1526.918|          0.889|
|MCF10A     |elevated     |  5.00|        1958.273|          0.788|
|MCF10A     |baseline     | 25.00|        7502.202|          0.870|
|MCF10A     |elevated     | 25.00|        3719.268|          0.916|
|MDA-MB-231 |baseline     |  0.00|          54.726|             NA|
|MDA-MB-231 |elevated     |  0.00|          66.852|             NA|
|MDA-MB-231 |baseline     |  0.10|         115.416|         -0.586|
|MDA-MB-231 |elevated     |  0.10|          60.015|         -0.797|
|MDA-MB-231 |baseline     |  0.25|         123.347|          0.565|
|MDA-MB-231 |elevated     |  0.25|         102.856|          0.330|

The timecourse summary shows this accounting model is not uniformly good. It looks promising in several non-MCF10A settings, but some conditions, especially parts of `MCF10A`, still have large trajectory error or even pathological correlation. So the result supports a basic biomass ODE as a plausible next branch, but not yet as a globally adequate description.

![plot of chunk cost-efficiency](figure/cost-efficiency-1.png)

## Notes

- `q50`, `q99`, and `tail_ratio = q99 / q50` are used as the main morphology summaries to avoid redundant quantile-derived features.
- The rate-coupling section is exploratory: it uses smoothed derivatives to estimate lead-lag structure between area change and death burden.
- The cost model still works with replicate-level count rows but well-level glucose summaries, matching the spirit of the original notebook while keeping the area data aligned to the actual count observations.
- `glucose_consumed` is estimated from the first and last glucose measurements within each well, not from nominal `G0`.
- Any strong effect seen here should still be treated as exploratory until promoted into a scripted analysis branch.
