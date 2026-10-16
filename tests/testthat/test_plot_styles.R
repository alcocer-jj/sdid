# sdid_plot has two styles. "minimal", the default, draws black trajectories with the synthetic control
# dashed, outlined stretches of time weight with a right-hand weight axis, a red dashed treatment line and dotted
# gridlines. "effects" adds the diff-in-diff diagram, with the trajectories at half opacity.

style.setup = function() {
  data(california_prop99)
  s = panel.matrices(california_prop99)
  e = synthdid_estimate(s$Y, s$N0, s$T0)
  list(s = s, e = e, ests = list(sdid = e, sc = sc_estimate(s$Y, s$N0, s$T0), did = did_estimate(s$Y, s$N0, s$T0)))
}

geoms.of = function(p) { unique(vapply(p$layers, function(l) class(l$geom)[1], "")) }
layer.data = function(p, geom) {
  b = ggplot2::ggplot_build(p)
  b$data[[which(vapply(b$plot$layers, function(l) class(l$geom)[1], "") == geom)[1]]]
}
weight.axis = function(p) {
  ax = ggplot2::ggplot_build(p)$plot$scales$get_scales("y")$secondary.axis
  if (is.null(ax) || inherits(ax, "waiver")) NULL else ax
}
# TRUE if a right-hand axis is actually drawn, which is what a reader sees, as opposed to what the scale declares.
right.axis.drawn = function(p) {
  g = ggplot2::ggplotGrob(p)
  any(vapply(g$grobs[grepl("^axis-r", g$layout$name)], function(x) !inherits(x, "zeroGrob"), logical(1)))
}
# TRUE for each panel whose band starts on one of the labeled values of its left axis.
bands.on.breaks = function(p) {
  b = ggplot2::ggplot_build(p)
  r = b$data[[which(vapply(b$plot$layers, function(l) class(l$geom)[1], "") == "GeomRibbon")[1]]]
  vapply(split(r, as.character(r$PANEL)), function(d) {
    breaks = b$layout$panel_params[[as.integer(as.character(d$PANEL[1]))]]$y$breaks
    any(abs(breaks - min(d$ymin)) < 1e-9, na.rm = TRUE)
  }, logical(1))
}

test_that("minimal is the default and draws no diff-in-diff diagram", {
  st = style.setup()
  expect_identical(ggplot2::ggplot_build(sdid_plot(st$e))$data, ggplot2::ggplot_build(sdid_plot(st$e, style = "minimal"))$data)
  expect_true(all(geoms.of(sdid_plot(st$e)) %in% c("GeomRibbon", "GeomLine", "GeomVline", "GeomCustomAnn", "GeomPath")))
  expect_true(all(c("GeomPoint", "GeomSegment", "GeomCurve") %in% geoms.of(sdid_plot(st$e, style = "effects"))))
  expect_false("GeomCurve" %in% geoms.of(plot(st$e)))
  expect_true("GeomCurve" %in% geoms.of(plot(st$e, style = "effects")))
  expect_error(sdid_plot(st$e, style = "classic"))
})

test_that("the effects style draws the diagram in black, the counterfactual dotted and the arrow red", {
  st = style.setup()
  p = sdid_plot(st$e, style = "effects")
  lines = layer.data(p, "GeomLine")
  expect_equal(unique(lines$alpha), .5)
  expect_setequal(unique(lines$linetype), c("solid", "dashed"))
  b = ggplot2::ggplot_build(p)
  geoms = vapply(b$plot$layers, function(l) class(l$geom)[1], "")
  segments = do.call(rbind, lapply(b$data[geoms == "GeomSegment"], function(d) d[, c("x", "xend", "linetype", "colour")]))
  # Two trend segments, solid, and the counterfactual, dotted. No vertical segments close the quadrilateral.
  expect_equal(nrow(segments), 3)
  expect_true(all(segments$x != segments$xend))
  expect_equal(sort(segments$linetype), c("dotted", "solid", "solid"))
  expect_true(all(segments$colour == "black"))
  points = do.call(rbind, lapply(b$data[geoms == "GeomPoint"], function(d) d[, c("colour", "alpha")]))
  expect_true(all(points$colour == "black"))
  expect_true(all(points$alpha == 1))
  arrow = layer.data(p, "GeomCurve")
  expect_equal(nrow(arrow), 1)
  expect_identical(arrow$colour, "red")
  # The arrow runs from the counterfactual to the actual post-treatment average, a drop of the estimated effect.
  expect_equal(arrow$yend - arrow$y, c(st$e), tolerance = 1e-8)
  # An SC estimate has no time weights, so its segment layers draw nothing and only the arrow and its endpoints remain.
  sc = ggplot2::ggplot_build(sdid_plot(st$ests$sc, style = "effects"))
  sc.geoms = vapply(sc$plot$layers, function(l) class(l$geom)[1], "")
  expect_equal(sum(vapply(sc$data[sc.geoms == "GeomSegment"], nrow, 1)), 0)
  expect_equal(nrow(sc$data[[which(sc.geoms == "GeomCurve")]]), 1)
})

test_that("the minimal style draws black trajectories, outlined bands and a red dashed treatment line", {
  st = style.setup()
  p = sdid_plot(st$e)
  lines = layer.data(p, "GeomLine")
  expect_true(all(lines$colour == "black"))
  expect_true(all(lines$alpha == 1))
  expect_setequal(unique(lines$linetype), c("solid", "dashed"))
  ribbon = layer.data(p, "GeomRibbon")
  expect_true(all(ribbon$fill == "black"))
  expect_equal(unique(ribbon$alpha), .34)
  ribbon.layer = p$layers[[which(vapply(p$layers, function(l) class(l$geom)[1], "") == "GeomRibbon")]]
  expect_identical(ribbon.layer$geom_params$outline.type, "full")
  vline = layer.data(p, "GeomVline")
  expect_identical(unique(vline$colour), "red")
  expect_identical(unique(vline$linetype), "dashed")
  expect_equal(unique(vline$alpha), .8)
  expect_equal(unique(vline$xintercept), mean(as.numeric(colnames(st$s$Y)[st$s$T0 + 0:1])))
  expect_equal(unique(layer.data(sdid_plot(st$e, onset.alpha = .5), "GeomVline")$alpha), .5)
})

test_that("bands are drawn only over positive weights and the stretches beside them", {
  st = style.setup()
  lambda = attr(st$e, "weights")$lambda
  years = as.numeric(colnames(st$s$Y))[seq_along(lambda)]
  ribbon = layer.data(sdid_plot(st$e), "GeomRibbon")
  positive = years[lambda > 0]
  beside = years[intersect(c(which(lambda > 0) - 1, which(lambda > 0) + 1), seq_along(years))]
  expect_true(all(ribbon$x %in% c(positive, beside)))
  expect_true(all(positive %in% ribbon$x))
  for (stretch in split(ribbon, ribbon$group)) { expect_true(any(stretch$ymax > stretch$ymin)) }
})

test_that("the weight axis maps each ribbon back to its time weight", {
  st = style.setup()
  lambda = attr(st$e, "weights")$lambda
  years = as.numeric(colnames(st$s$Y))[seq_along(lambda)]
  for (comparable in c(FALSE, TRUE)) {
    p = sdid_plot(st$e, lambda.comparable = comparable)
    ax = weight.axis(p)
    expect_false(is.null(ax), info = paste("comparable", comparable))
    to.weight = if (!is.null(ax$transform)) ax$transform else ax$trans
    ribbon = layer.data(p, "GeomRibbon")
    expect_equal(to.weight(ribbon$ymax), lambda[match(ribbon$x, years)], tolerance = 1e-10, info = paste("comparable", comparable))
    expect_equal(unique(to.weight(ribbon$ymin)), 0, tolerance = 1e-10)
    expect_true(all(ax$breaks >= 0 & ax$breaks <= max(lambda) + 1e-9))
  }
})

test_that("each band starts on a labeled value of the left axis", {
  st = style.setup()
  data(CPS)
  cps = panel.matrices(CPS, unit = "state", time = "year", outcome = "log_wage", treatment = "min_wage")
  plots = list(single = sdid_plot(st$e), comparable = sdid_plot(st$e, lambda.comparable = TRUE),
               cps = sdid_plot(synthdid_estimate(cps$Y, cps$N0, cps$T0)), facets = sdid_plot(st$ests),
               placebo = sdid_placebo_plot(st$e), placebo.1975 = sdid_placebo_plot(st$e, placebo.time = 1975))
  for (name in names(plots)) { expect_true(all(bands.on.breaks(plots[[name]])), info = name) }
  expect_equal(min(layer.data(plots$single, "GeomRibbon")$ymin), 0)
})

test_that("the axes have loose ticks, as the Stata command chooses them", {
  st = style.setup()
  data(CPS)
  s = panel.matrices(CPS, unit = "state", time = "year", outcome = "log_wage", treatment = "abort_ban")
  cps = synthdid_estimate(s$Y, s$N0, s$T0)
  # The ticks and labels drawn, and the range of the data, leaving out the range frame's lines.
  drawn = function(p, axis) {
    b = ggplot2::ggplot_build(p)
    pp = b$layout$panel_params[[1]]
    data = b$data[vapply(b$plot$layers, function(l) class(l$geom)[1], "") != "GeomPath"]
    values = unlist(lapply(data, function(d) if (axis == "x") c(d$x, d$xend) else c(d$y, d$ymin, d$ymax, d$yend)))
    keep = !is.na(pp[[axis]]$breaks)
    list(at = pp[[axis]]$breaks[keep], labels = pp[[axis]]$get_labels()[keep], data = range(values[is.finite(values)]),
         ends = range(pp[[axis]]$breaks, na.rm = TRUE))
  }
  for (estimate in list(st$e, cps)) {
    for (axis in c("x", "y")) {
      loose = drawn(sdid_plot(estimate), axis)
      # Loose ticks enclose the data, are those of the loose form of ggplot2's algorithm, and are evenly spaced.
      expect_lte(loose$at[1], loose$data[1])
      expect_gte(loose$at[length(loose$at)], loose$data[2])
      # The time axis, whose periods are years, leaves out the step 2.5.
      steps = if (axis == "x") { c(1, 5, 2, 4, 3) } else { c(1, 5, 2, 2.5, 4, 3) }
      expect_equal(loose$at, scales::extended_breaks(Q = steps, only.loose = TRUE)(loose$data))
      expect_equal(length(unique(round(diff(loose$at), 10))), 1)
    }
  }
  # Prop 99 and CPS in particular.
  expect_equal(drawn(sdid_plot(st$e), "y")$at, c(0, 40, 80, 120, 160))
  expect_equal(drawn(sdid_plot(cps), "x")$at, seq(1970, 2020, by = 10))
  expect_equal(range(drawn(sdid_placebo_plot(st$e), "x")$at), c(1970, 2000))
  # The choice of tick styles, facets and their direction is gone.
  for (removed in c("ticks", "facet", "facet.vertical")) {
    expect_error(do.call(sdid_plot, stats::setNames(list(st$e, 1), c("estimates", removed))), "unused argument", info = removed)
  }
  expect_error(sdid_placebo_plot(st$e, ticks = "loose"), "unused argument")
})

test_that("stacked panels each get an outcome axis and a weight axis, and only the bottom one a time axis", {
  st = style.setup()
  p = sdid_plot(st$ests)
  b = ggplot2::ggplot_build(p)
  geoms = vapply(b$plot$layers, function(l) class(l$geom)[1], "")
  layout = b$layout$layout
  # Panels follow the order of the estimates, with a title above each, as they are not a staggered design's cohorts.
  expect_identical(as.character(layout$facet), names(st$ests))
  expect_identical(p$theme$strip.text$hjust, 0)
  frame = b$data[[which(geoms == "GeomPath")[1]]]
  for (panel in layout$PANEL) {
    in.panel = frame[frame$PANEL == panel, ]
    expect_true(any(in.panel$x == -Inf), info = panel)
    expect_equal(any(in.panel$y == -Inf), panel == layout$PANEL[which.max(layout$ROW)], info = panel)
  }
  # A weight axis for each panel with a band, the synthetic control panel having none.
  axes = b$data[[which(geoms == "GeomSdidWeightAxis")]]
  banded = unique(layer.data(p, "GeomRibbon")$PANEL)
  expect_setequal(unique(axes$PANEL), banded)
  for (panel in unique(axes$PANEL)) {
    a = axes[axes$PANEL == panel, ]
    expect_equal(a$label[1], "0.00")
    expect_equal(min(a$y), min(a$ymin))
  }
  expect_false("GeomSdidWeightAxis" %in% geoms.of(sdid_plot(st$ests, lambda.axis = FALSE)))
})

test_that("a staggered design's cohorts are labeled at their adoption periods under each panel", {
  set.seed(1)
  N = 23; periods = 15; adoption = c(rep(NA, 14), rep(c(6, 9, 12), each = 3))
  d = expand.grid(time = seq_len(periods), unit = sprintf("u%02d", seq_len(N)), stringsAsFactors = FALSE)
  i = match(d$unit, sprintf("u%02d", seq_len(N)))
  d$treat = as.integer(!is.na(adoption[i]) & d$time >= adoption[i])
  d$y = rnorm(N)[i] + .3 * d$time + rnorm(nrow(d), sd = .2) + 2 * d$treat
  fit = sdid(d, "y", "unit", "time", "treat", vce = "noinference", graph = FALSE, quietly = TRUE)
  p = sdid_plot(fit)
  b = ggplot2::ggplot_build(p)
  geoms = vapply(b$plot$layers, function(l) class(l$geom)[1], "")
  expect_identical(as.character(b$layout$layout$facet), c("6", "9", "12"))
  expect_s3_class(p$theme$strip.text, "element_blank")
  labels = b$data[[which(geoms == "GeomSdidCohortLabel")]]
  expect_equal(labels$x, c(6, 9, 12))
  expect_identical(labels$label, c("6", "9", "12"))
  expect_identical(labels$tick, c(FALSE, FALSE, TRUE))
  # Each treatment line is midway between the last period before adoption and the adoption period.
  expect_equal(sort(unique(layer.data(p, "GeomVline")$xintercept)), c(5.5, 8.5, 11.5))
  # Every time axis label is kept, the bottom panel's label moving instead if it would overlap one.
  x = b$layout$panel_params[[3]]$x
  shown = x$get_labels()[!is.na(x$breaks)]
  expect_equal(x$breaks[!is.na(x$breaks)], c(0, 5, 10, 15))
  expect_identical(shown, c("0", "5", "10", "15"))
})

test_that("the secondary weight axis and its full title appear only for a single panel", {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  st = style.setup()
  single = sdid_plot(st$e)
  expect_false(is.null(weight.axis(single)))
  expect_true("GeomCustomAnn" %in% geoms.of(single))
  expect_identical(single$coordinates$clip, "off")
  expect_true(right.axis.drawn(single))
  expect_true(right.axis.drawn(sdid_plot(st$e, style = "effects")))
  # Without a weight axis no right-hand axis may be drawn. ggplot2 3.5 and later once drew a copy of the left axis here.
  for (p in list(sdid_plot(st$ests), sdid_plot(st$e, lambda.axis = FALSE), sdid_plot(st$ests$sc),
                 sdid_plot(st$e, overlay = 1), sdid_placebo_plot(st$e), sdid_plot(st$ests, style = "effects"))) {
    expect_null(weight.axis(p))
    expect_false("GeomCustomAnn" %in% geoms.of(p))
    expect_false(right.axis.drawn(p))
  }
})

test_that("the minimal style sets the documented theme and frame", {
  st = style.setup()
  single = sdid_plot(st$e)
  expect_s3_class(single$theme$panel.grid.major, "element_line")
  expect_identical(single$theme$panel.grid.major$linetype, "dotted")
  expect_identical(single$theme$panel.grid.major$colour, "lightgray")
  expect_equal(single$theme$panel.grid.major$linewidth, .6)
  expect_s3_class(single$theme$panel.grid.minor, "element_blank")
  expect_equal(grid::convertUnit(single$theme$axis.ticks.length, "cm", valueOnly = TRUE), .25)
  expect_equal(single$theme$axis.text$size, 10)
  expect_equal(single$theme$legend.text$size, 12)
  expect_s3_class(single$theme$panel.border, "element_blank")
  expect_s3_class(single$theme$axis.line, "element_line")
  # Each axis ends at its outermost ticks, which enclose the data, padded by 2% of the axis on each side.
  for (p in list(single, sdid_plot(st$e, style = "effects"), sdid_plot(st$e, lambda.axis = FALSE))) {
    pp = ggplot2::ggplot_build(p)$layout$panel_params[[1]]
    for (axis in c("x", "y")) {
      ticks = range(pp[[axis]]$breaks, na.rm = TRUE)
      expect_equal(pp[[paste0(axis, ".range")]], ticks + c(-1, 1) * .02 * diff(ticks), info = axis)
    }
  }
  # Panels stacked in a column are open.
  expect_s3_class(sdid_plot(st$ests)$theme$panel.border, "element_blank")
})

test_that("a single panel's axis lines span the data and the outermost ticks", {
  # The frame lines, drawn at the bottom and left edges of the panel.
  frame = function(p) {
    b = ggplot2::ggplot_build(p)
    paths = b$data[vapply(b$plot$layers, function(l) class(l$geom)[1], "") == "GeomPath"]
    bottom = Filter(function(d) all(d$y == -Inf), paths)
    left = Filter(function(d) all(d$x == -Inf), paths)
    list(x = if (length(bottom)) range(bottom[[1]]$x), y = if (length(left)) range(left[[1]]$y), b = b)
  }
  span = function(b, axis) {
    values = unlist(lapply(b$data, function(d) if (axis == "x") c(d$x, d$xend, d$xintercept) else c(d$y, d$ymin, d$ymax, d$yend)))
    breaks = b$layout$panel_params[[1]][[axis]]$breaks
    range(c(values, breaks)[is.finite(c(values, breaks))])
  }
  st = style.setup()
  data(CPS)
  cps = sdid:::panel.matrices(CPS, unit = "state", time = "year", outcome = "log_wage", treatment = "abort_ban")
  for (p in list(sdid_plot(st$e), sdid_plot(st$e, style = "effects"),
                 sdid_plot(sdid:::synthdid_estimate(cps$Y, cps$N0, cps$T0)))) {
    f = frame(p)
    expect_equal(f$x, span(f$b, "x"))
    expect_equal(f$y, span(f$b, "y"))
    expect_s3_class(p$theme$axis.line.x.bottom, "element_blank")
    expect_s3_class(p$theme$axis.line.y.left, "element_blank")
  }
  # On CPS the ticks run from 1970 to 2020, enclosing the data, so the frame runs between them.
  f = frame(sdid_plot(sdid:::synthdid_estimate(cps$Y, cps$N0, cps$T0)))
  expect_equal(f$x, c(1970, 2020))
  # Faceted plots keep their box and draw no frame.
  expect_null(frame(sdid_plot(st$ests))$x)
})

test_that("a time axis of whole-number periods gets whole-number ticks", {
  panel = function(times, adoption, N = 30, seed = 4) {
    set.seed(seed)
    d = expand.grid(time = times, unit = sprintf("s%02d", 1:N), stringsAsFactors = FALSE)
    a = c(rep(NA, N - length(adoption)), adoption)[match(d$unit, sprintf("s%02d", 1:N))]
    d$treat = as.integer(!is.na(a) & d$time >= a)
    d$y = .3 + .002 * (d$time - times[1]) + rnorm(nrow(d), sd = .02) + .08 * d$treat
    sdid(d, "y", "unit", "time", "treat", vce = "noinference", graph = FALSE, quietly = TRUE)
  }
  ticks = function(p, axis = "x") {
    pp = ggplot2::ggplot_build(p)$layout$panel_params
    b = pp[[length(pp)]][[axis]]$breaks
    b[!is.na(b)]
  }
  # Years from 2003 to 2022, which the step 2.5 would label 2002.5, 2007.5 and so on.
  years = panel(2003:2022, c(2009, 2013))
  expect_equal(scales::extended_breaks(only.loose = TRUE)(c(2003, 2022)), seq(2002.5, 2022.5, by = 5))
  expect_equal(ticks(sdid_plot(years)), seq(2002, 2022, by = 4))
  at = ticks(sdid_plot(years))
  expect_true(all(at == round(at)))
  # A short panel whose ticks the algorithm could put at half periods takes whole-number steps instead.
  short = panel(2017:2019, 2019, N = 12)
  expect_equal(ticks(sdid_plot(short)), 2017:2019)
  # The outcome axis keeps every step, and periods that are not whole numbers keep 2.5 on the time axis.
  expect_true(any(abs(ticks(sdid_plot(years), "y") * 100 - round(ticks(sdid_plot(years), "y") * 100)) < 1e-9))
  halves = panel(seq(0.5, 10, by = .5), 6, N = 12)
  expect_equal(ticks(sdid_plot(halves)), scales::extended_breaks(only.loose = TRUE)(c(.5, 10)))
})

test_that("a cohort label moves just enough to clear the time axis's labels, or below them if it cannot", {
  # Draws the label in a panel of the given width, in points, and returns where its text lands and how far down.
  placed = function(at, neighbors, labels, width = 400) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off())
    grid::pushViewport(grid::viewport(width = grid::unit(width, "pt"), height = grid::unit(200, "pt")))
    g = grid::makeContent(sdid:::cohort_label_grob("2013", at, TRUE, neighbors, labels))
    text = Filter(function(child) inherits(child, "text"), g$children)
    if (length(text) == 0) { return(NULL) }
    half = function(s) { grid::convertWidth(grid::grobWidth(grid::textGrob(s, gp = grid::gpar(fontsize = 10))), "pt", valueOnly = TRUE) / 2 }
    list(x = grid::convertX(text[[1]]$x, "pt", valueOnly = TRUE), y = grid::convertY(text[[1]]$y, "pt", valueOnly = TRUE),
         room = half("2013") + half("2014") + 5.5)
  }
  ticks = c(2002, 2006, 2010, 2014, 2018, 2022)
  npc = (ticks - 2001.6) / (2022.4 - 2001.6)
  at = (2013 - 2001.6) / (2022.4 - 2001.6)
  labels = as.character(ticks)
  # At 450 points the label would overlap 2014 but fits before it, so it moves away by just enough, its tick staying
  # at 2013. At 600 points it already clears 2014 and stays.
  wide = placed(at, npc, labels, width = 450)
  expect_lt(wide$x, at * 450)
  expect_equal(abs(wide$x - npc[4] * 450), wide$room, tolerance = 1e-6)
  expect_equal(placed(at, npc, labels, width = 600)$x, at * 600)
  # Far from any label, it stays at its tick.
  alone = placed(at, npc[c(1, 6)], labels[c(1, 6)], width = 600)
  expect_equal(alone$x, at * 600)
  # Where there is no room between 2010 and 2014, it goes below them, under its tick.
  narrow = placed(at, npc, labels, width = 200)
  expect_equal(narrow$x, at * 200)
  expect_lt(narrow$y, wide$y - 11)
  # A tick at the adoption period already labels it, so the cohort's label is not repeated.
  expect_null(placed(npc[3], npc, labels))
})

test_that("user ticks replace the automatic ones, with the design following them", {
  st = style.setup()
  axes = function(p) { b = ggplot2::ggplot_build(p); list(pp = b$layout$panel_params[[1]], b = b) }
  given = seq(1965, 2005, by = 5)
  a = axes(sdid_plot(st$e, x.breaks = given, y.breaks = c(0, 50, 100, 150, 175)))
  expect_equal(a$pp$x$breaks[!is.na(a$pp$x$breaks)], given)
  expect_equal(a$pp$x.range, range(given) + c(-1, 1) * .02 * diff(range(given)))
  expect_equal(a$pp$y$breaks[!is.na(a$pp$y$breaks)], c(0, 50, 100, 150, 175))
  # The frame spans the given ticks, and the band sits on the given tick of 0.
  geoms = vapply(a$b$plot$layers, function(l) class(l$geom)[1], "")
  bottom = Filter(function(d) all(d$y == -Inf), a$b$data[geoms == "GeomPath"])[[1]]
  expect_equal(range(bottom$x), range(given))
  expect_equal(min(a$b$data[[which(geoms == "GeomRibbon")[1]]]$ymin), 0)
  # Ticks inside the data still give an axis covering it.
  inside = axes(sdid_plot(st$e, x.breaks = c(1980, 1990)))
  expect_equal(inside$pp$x$breaks[!is.na(inside$pp$x$breaks)], c(1980, 1990))
  expect_equal(inside$pp$x.range, c(1970, 2000) + c(-1, 1) * .02 * 30)
  # The placebo plot passes them on, and other values are rejected.
  expect_equal(axes(sdid_placebo_plot(st$e, x.breaks = given))$pp$x$breaks[!is.na(axes(sdid_placebo_plot(st$e, x.breaks = given))$pp$x$breaks)], given)
  expect_error(sdid_plot(st$e, x.breaks = c(1970, NA)), "finite numbers")
  expect_error(sdid_plot(st$e, y.breaks = "a"), "finite numbers")
  expect_error(sdid_plot(st$e, y.breaks = as.Date("2000-01-01")), "y.breaks must be numbers")
})

test_that("the weight axis stops at its ticks", {
  st = style.setup()
  stored = function(p) { g = p$guides; if (is.environment(g) || inherits(g, "ggproto")) { g$guides } else { g } }
  single = stored(sdid_plot(st$e))
  if (utils::packageVersion("ggplot2") >= "3.5.0") {
    # Only the weight axis has a guide, which caps it. The outcome and time axes keep ggplot2's default, uncapped.
    expect_identical(names(single), "y.sec")
    expect_identical(single$y.sec$params$cap, "both")
  } else {
    expect_length(single, 0)
  }
  expect_length(stored(sdid_plot(st$e, lambda.axis = FALSE)), 0)
  expect_length(stored(sdid_plot(st$ests)), 0)
})

test_that("each estimate gets its own panel, labeled by its name, in either style", {
  st = style.setup()
  for (style in c("minimal", "effects")) {
    expect_error(sdid_plot(list(a = st$e, a = st$e), style = style), "names of the list must not repeat")
  }
  expect_error(sdid_placebo_plot(st$e, overlay = TRUE))
  # The list's names label the panels, in the order given.
  renamed = stats::setNames(st$ests, c("c", "a", "b"))
  expect_identical(as.character(ggplot2::ggplot_build(sdid_plot(renamed))$layout$layout$facet), c("c", "a", "b"))
  # The time weights are on a common scale only when asked for.
  expect_false(formals(sdid_plot)$lambda.comparable)
})

test_that("a named list of sdid results is plotted as a list of their estimates", {
  data(california_prop99)
  fit = function(method) {
    sdid(california_prop99, "PacksPerCapita", "State", "Year", "treated", method = method, vce = "noinference",
         graph = FALSE, quietly = TRUE)
  }
  a = fit("sdid"); b = fit("sc")
  for (plot in list(sdid_plot, sdid_units_plot)) {
    from.fits = ggplot2::ggplot_build(plot(list(SDID = a, SC = b)))
    from.estimates = ggplot2::ggplot_build(plot(list(SDID = a$estimates[[1]], SC = b$estimates[[1]])))
    expect_identical(from.fits$data, from.estimates$data)
    expect_equal(nrow(from.fits$layout$layout), 2)
    # A fit and an estimate can be mixed.
    expect_identical(ggplot2::ggplot_build(plot(list(SDID = a, SC = b$estimates[[1]])))$data, from.estimates$data)
  }
  # A staggered result in a list stops, as its cohorts would each need a panel.
  set.seed(1)
  N = 23; periods = 15; adoption = c(rep(NA, 14), rep(c(6, 9, 12), each = 3))
  d = expand.grid(time = seq_len(periods), unit = sprintf("u%02d", seq_len(N)), stringsAsFactors = FALSE)
  i = match(d$unit, sprintf("u%02d", seq_len(N)))
  d$treat = as.integer(!is.na(adoption[i]) & d$time >= adoption[i])
  d$y = rnorm(N)[i] + .3 * d$time + rnorm(nrow(d), sd = .2) + 2 * d$treat
  staggered = sdid(d, "y", "unit", "time", "treat", vce = "noinference", graph = FALSE, quietly = TRUE)
  for (plot in list(sdid_plot, sdid_units_plot)) {
    expect_error(plot(list(A = a, B = staggered)), "single adoption period")
  }
})

test_that("the RMSE plot has been removed", {
  expect_false(exists("synthdid_rmse_plot", envir = asNamespace("sdid")))
})

test_that("plots no longer take standard error or classic-style arguments", {
  st = style.setup()
  removed = list(se.method = "none", trajectory.linetype = 1, trajectory.alpha = .5, guide.linetype = 2,
                 diagram.alpha = .95, effect.alpha = .95, ci.alpha = .3)
  for (name in names(removed)) {
    expect_error(do.call(sdid_plot, c(list(st$e), removed[name])), "unused argument", info = name)
  }
  expect_error(sdid_units_plot(st$e, se.method = "none"), "unused argument")
  for (removed in c("negligible.threshold", "negligible.alpha", "units")) {
    expect_error(do.call(sdid_units_plot, stats::setNames(list(st$e, 1), c("estimates", removed))), "unused argument", info = removed)
  }
})

test_that("the units plot shows each control's own estimate, whose weighted average is the estimate", {
  st = style.setup()
  # The controls in the scatter's point layers, weighted and unweighted.
  points = function(p) {
    b = ggplot2::ggplot_build(p)
    geoms = vapply(b$plot$layers, function(l) class(l$geom)[1], "")
    layers = Filter(function(d) nrow(d) > 0, b$data[geoms == "GeomPoint"])
    do.call(rbind, lapply(layers, function(d) d[, c("x", "y", "PANEL")]))
  }
  for (estimate in list(st$e, st$ests$sc, st$ests$did)) {
    pts = points(sdid_units_plot(estimate))
    expect_equal(nrow(pts), st$s$N0)
    expect_equal(sum(pts$x * pts$y), c(estimate), tolerance = 1e-10)
    expect_setequal(pts$x, attr(estimate, "weights")$omega)
  }
  # Zero-weight controls are open circles, and weighted.only leaves them out.
  b = ggplot2::ggplot_build(sdid_units_plot(st$e))
  geoms = vapply(b$plot$layers, function(l) class(l$geom)[1], "")
  shapes = vapply(b$data[geoms == "GeomPoint"], function(d) d$shape[1], 0)
  expect_true(1 %in% shapes)
  expect_equal(nrow(points(sdid_units_plot(st$e, weighted.only = TRUE))), sum(attr(st$e, "weights")$omega > 0))
  # The effect is a red line, with its label on it.
  expect_equal(unique(layer.data(sdid_units_plot(st$e), "GeomHline")$yintercept), c(st$e))
  expect_identical(layer.data(sdid_units_plot(st$e), "GeomLabel")$label, sprintf("ATT %.1f", c(st$e)))
})

test_that("label.units chooses which controls are named", {
  st = style.setup()
  omega = attr(st$e, "weights")$omega
  named = function(p, by.panel = FALSE) {
    b = ggplot2::ggplot_build(p)
    geoms = vapply(b$plot$layers, function(l) class(l$geom)[1], "")
    d = if ("GeomTextRepel" %in% geoms) { b$data[[which(geoms == "GeomTextRepel")]] } else { b$data[[which(geoms == "GeomText")[1]]] }
    d = d[d$label != "", ]
    if (by.panel) { stats::setNames(as.vector(table(factor(d$PANEL, levels = b$layout$layout$PANEL))), as.character(b$layout$layout$panel)) } else { d$label }
  }
  expect_length(named(sdid_units_plot(st$e)), st$s$N0)
  expect_length(named(sdid_units_plot(st$e, label.units = "weighted")), sum(omega > 0))
  top = rownames(st$s$Y)[order(omega, decreasing = TRUE)[1:5]]
  expect_setequal(named(sdid_units_plot(st$e, label.units = 5)), top)
  expect_length(named(sdid_units_plot(st$e, label.units = 500)), sum(omega > 0))
  expect_length(named(sdid_units_plot(st$e, label.units = 0)), 0)
  # Several panels default to ten names each, any tied at the cutoff included, as all of DiD's equal weights are.
  per.panel = named(sdid_units_plot(st$ests), by.panel = TRUE)
  expect_equal(per.panel[["sdid"]], 10)
  expect_equal(per.panel[["sc"]], min(10, sum(attr(st$ests$sc, "weights")$omega > 0)))
  expect_equal(per.panel[["did"]], st$s$N0)
  for (bad in list("some", -1, 2.5, c(1, 2), NA)) { expect_error(sdid_units_plot(st$e, label.units = bad), "label.units must be") }
})

test_that("the bars run from the effect to each weighted control's estimate, largest weight at the top", {
  st = style.setup()
  omega = attr(st$e, "weights")$omega
  p = sdid_units_plot(st$e, type = "bars")
  b = ggplot2::ggplot_build(p)
  segments = layer.data(p, "GeomSegment")
  expect_equal(nrow(segments), sum(omega > 0))
  expect_true(all(abs(segments$x - c(st$e)) < 1e-9))
  # From the top of the plot down, the weights decrease.
  units = sub("^.*\r", "", levels(b$plot$data$id))
  shown = stats::setNames(omega, rownames(st$s$Y)[1:st$s$N0])[rev(units)]
  expect_false(is.unsorted(rev(shown)))
  expect_equal(unname(shown[1]), max(omega))
  # The weights are written at the bars' ends.
  expect_setequal(layer.data(p, "GeomText")$label[-nrow(layer.data(p, "GeomText"))],
                  formatC(omega[omega > 0], format = "f", digits = 3))
  # The effect has its own unlabeled tick, and the ticks cover every bar.
  x = b$layout$panel_params[[1]]$x
  expect_true(any(abs(x$breaks - c(st$e)) < 1e-9, na.rm = TRUE))
  expect_identical(x$get_labels()[which(abs(x$breaks - c(st$e)) < 1e-9)], "")
  ticks = x$breaks[!is.na(x$breaks) & abs(x$breaks - c(st$e)) > 1e-9]
  expect_lte(min(ticks), min(segments$xend))
  expect_gte(max(ticks), max(segments$xend))
  # label.units chooses the bars drawn.
  expect_equal(nrow(layer.data(sdid_units_plot(st$e, type = "bars", label.units = 6), "GeomSegment")), 6)
})

test_that("a staggered fit is plotted cohort by cohort, or one cohort alone", {
  set.seed(1)
  N = 23; periods = 15; adoption = c(rep(NA, 14), rep(c(6, 9, 12), each = 3))
  d = expand.grid(time = seq_len(periods), unit = sprintf("u%02d", seq_len(N)), stringsAsFactors = FALSE)
  i = match(d$unit, sprintf("u%02d", seq_len(N)))
  d$treat = as.integer(!is.na(adoption[i]) & d$time >= adoption[i])
  d$y = rnorm(N)[i] + .3 * d$time + rnorm(nrow(d), sd = .2) + 2 * d$treat
  fit = sdid(d, "y", "unit", "time", "treat", vce = "noinference", graph = FALSE, quietly = TRUE)
  for (type in c("scatter", "bars")) {
    b = ggplot2::ggplot_build(sdid_units_plot(fit, type = type))
    expect_identical(as.character(b$layout$layout$panel), c("6", "9", "12"), info = type)
    one = ggplot2::ggplot_build(sdid_units_plot(fit, type = type, cohort = 9))
    expect_equal(nrow(one$layout$layout), 1)
  }
  expect_equal(unique(layer.data(sdid_units_plot(fit, cohort = "9"), "GeomHline")$yintercept), c(fit$estimates[["9"]]))
  expect_error(sdid_units_plot(fit, cohort = 7), "one of the adoption periods, 6, 9, 12")
  data(california_prop99)
  block = sdid(california_prop99, "PacksPerCapita", "State", "Year", "treated", vce = "noinference", graph = FALSE, quietly = TRUE)
  expect_error(sdid_units_plot(block, cohort = 1989), "staggered design")
  expect_equal(nrow(ggplot2::ggplot_build(sdid_units_plot(list(a = block, b = block)))$layout$layout), 2)
  expect_error(sdid_units_plot(block, type = "lines"))
})

test_that("building a plot does not open a graphics device", {
  st = style.setup()
  before = grDevices::dev.list()
  invisible(ggplot2::ggplot_build(sdid_plot(st$e)))
  expect_identical(grDevices::dev.list(), before)
})

test_that("minimal plots render without warnings or messages", {
  # Printing to the PDF device checks that every label can be drawn there, including the time-weight axis title.
  # The null file leaves no Rplots.pdf behind. ggplot2 4.0 and later also report titles set for unused aesthetics.
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  st = style.setup()
  cp = california_prop99
  cp$date = as.Date(sprintf("%04d/01/01", cp$Year))
  sd = panel.matrices(cp, time = "date")
  plots = list(sdid_plot(st$e), sdid_plot(st$ests), sdid_plot(st$e, lambda.comparable = TRUE),
               sdid_plot(st$e, spaghetti.units = rownames(st$s$Y)[1:3]), sdid_plot(st$ests$sc),
               sdid_plot(st$ests$did), sdid_plot(synthdid_estimate(sd$Y, sd$N0, sd$T0)),
               sdid_placebo_plot(st$e), sdid_placebo_plot(st$e, placebo.time = 1975),
               sdid_plot(st$e, style = "effects"), sdid_plot(st$ests, style = "effects"),
               sdid_plot(st$ests$sc, style = "effects"), sdid_plot(st$e, style = "effects", overlay = 1),
               sdid_placebo_plot(st$e, style = "effects", placebo.time = 1975),
               sdid_units_plot(st$e), sdid_units_plot(st$ests), sdid_units_plot(st$ests, weighted.only = TRUE),
               sdid_units_plot(st$e, type = "bars"), sdid_units_plot(st$ests, type = "bars"))
  for (p in plots) { expect_no_message(expect_no_warning(print(p))) }
})

test_that("the placebo plot passes the style and placebo time through", {
  st = style.setup()
  expect_false("GeomCurve" %in% geoms.of(sdid_placebo_plot(st$e)))
  expect_true("GeomCurve" %in% geoms.of(sdid_placebo_plot(st$e, style = "effects")))
  vlines = layer.data(sdid_placebo_plot(st$e, placebo.time = 1975), "GeomVline")
  # Each line is midway between the last period before treatment and the first period treated.
  expect_setequal(vlines$xintercept, c(1988.5, 1974.5))
})
