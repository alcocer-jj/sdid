# Column names used inside ggplot2::aes(). They are data frame columns, not
# variables, but R CMD check cannot tell the difference.
utils::globalVariables(c(
  "color",
  "direction",
  "effect",
  "estimate",
  "facet",
  "frame",
  "group",
  "id",
  "label",
  "name",
  "own",
  "panel",
  "right",
  "shade",
  "size",
  "text",
  "tick",
  "iteration",
  "method",
  "role",
  "se",
  "show",
  "stretch",
  "unit",
  "weight",
  "x",
  "xend",
  "xintercept",
  "y",
  "yend",
  "ymax",
  "ymin"
))

# Stops with an informative message if ggplot2 is not installed. Plotting
# functions call ggplot2 through its namespace, so it is never attached to the
# user's search path.
require_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Plotting requires the package `ggplot2`. Install it to use this function."
    )
  }
}

# A layer drawing each panel's time-weight axis just outside its right edge,
# for panels stacked in a column, where one secondary axis cannot describe
# bands drawn on different scales. Each row gives a tick's position and label
# in y, and the span of the axis line in ymin and ymax. The line, ticks,
# labels and title are placed in absolute units from the panel's edge, as an
# axis is, with ticks 0.25 cm long and labels 2.2 points beyond them.
weight_axis_layer <- function(data) {
  geom <- ggplot2::ggproto(
    "GeomSdidWeightAxis",
    ggplot2::Geom,
    required_aes = c("y", "ymin", "ymax", "label"),
    draw_panel = function(data, panel_params, coord) {
      d <- coord$transform(data, panel_params)
      right <- grid::unit(1, "npc")
      tick <- grid::unit(.25, "cm")
      pad <- grid::unit(2.2, "pt")
      line <- grid::gpar(col = "black", lwd = .4 * ggplot2::.pt)
      labels <- grid::textGrob(
        d$label,
        x = right + tick + pad,
        y = grid::unit(d$y, "npc"),
        hjust = 0,
        gp = grid::gpar(fontsize = 10, col = "grey30")
      )
      grid::grobTree(
        grid::segmentsGrob(
          right,
          grid::unit(d$ymin[1], "npc"),
          right,
          grid::unit(d$ymax[1], "npc"),
          gp = line
        ),
        grid::segmentsGrob(
          right,
          grid::unit(d$y, "npc"),
          right + tick,
          grid::unit(d$y, "npc"),
          gp = line
        ),
        labels,
        grid::textGrob(
          expression(lambda),
          x = right +
            tick +
            pad +
            grid::grobWidth(labels) +
            grid::unit(5.5, "pt"),
          y = grid::unit(mean(c(d$ymin[1], d$ymax[1])), "npc"),
          rot = -90,
          vjust = 0,
          gp = grid::gpar(fontsize = 11)
        )
      )
    }
  )
  ggplot2::layer(
    geom = geom,
    stat = "identity",
    position = "identity",
    data = data,
    mapping = ggplot2::aes(y = y, ymin = ymin, ymax = ymax, label = label),
    inherit.aes = FALSE,
    show.legend = FALSE
  )
}

# A layer labeling each cohort's adoption period under its panel, where a time
# axis's labels would be, with a tick on the panel that has the time axis.
# That panel's label is moved, if it would overlap the time axis's labels in
# `neighbors`, a data frame of their positions and labels, just far enough to
# clear them, keeping its tick at the adoption period, or, if there is no room
# between them, to a line below them. The overlap depends on the size the plot
# is drawn at, so it is resolved when the label is drawn.
cohort_label_layer <- function(data, neighbors = NULL) {
  geom <- ggplot2::ggproto(
    "GeomSdidCohortLabel",
    ggplot2::Geom,
    required_aes = c("x", "label", "tick"),
    draw_panel = function(data, panel_params, coord) {
      d <- coord$transform(data, panel_params)
      near <- if (
        isTRUE(d$tick[1]) && !is.null(neighbors) && nrow(neighbors) > 0
      ) {
        coord$transform(data.frame(x = neighbors$x), panel_params)$x
      } else {
        numeric(0)
      }
      cohort_label_grob(
        d$label[1],
        d$x[1],
        isTRUE(d$tick[1]),
        near,
        if (length(near) > 0) {
          neighbors$label
        } else {
          character(0)
        }
      )
    }
  )
  ggplot2::layer(
    geom = geom,
    stat = "identity",
    position = "identity",
    data = data,
    mapping = ggplot2::aes(x = x, label = label, tick = tick),
    inherit.aes = FALSE,
    show.legend = FALSE
  )
}

# The label of a cohort's adoption period at `at`, a proportion of the panel's
# width, with a tick if `tick` is TRUE, and the positions and labels of the
# time axis's labels it must not overlap.
cohort_label_grob <- function(
  label,
  at,
  tick,
  neighbors = numeric(0),
  neighbor.labels = character(0)
) {
  grid::gTree(
    label = label,
    at = at,
    tick = tick,
    neighbors = neighbors,
    neighbor.labels = neighbor.labels,
    cl = "sdidcohortlabel"
  )
}

#' @exportS3Method grid::makeContent
makeContent.sdidcohortlabel <- function(x) {
  gp <- grid::gpar(fontsize = 10, col = "grey30")
  tick <- grid::unit(.25, "cm")
  y <- grid::unit(0, "npc") - tick - grid::unit(2.2, "pt")
  children <- if (isTRUE(x$tick)) {
    grid::gList(grid::segmentsGrob(
      grid::unit(x$at, "npc"),
      grid::unit(0, "npc"),
      grid::unit(x$at, "npc"),
      grid::unit(0, "npc") - tick,
      gp = grid::gpar(col = "black")
    ))
  } else {
    grid::gList()
  }
  keep <- !is.na(x$neighbors) &
    !is.na(x$neighbor.labels) &
    nzchar(x$neighbor.labels)
  neighbors <- x$neighbors[keep]
  neighbor.labels <- x$neighbor.labels[keep]
  # A time axis label at the adoption period already gives it, so the cohort's
  # label is not repeated.
  if (any(abs(neighbors - x$at) < 1e-9)) {
    return(grid::setChildren(x, children))
  }
  # Positions in points from the panel's left edge, and the room each pair of
  # labels needs, half of each width and a gap of 5.5 points.
  panel <- grid::convertWidth(grid::unit(1, "npc"), "pt", valueOnly = TRUE)
  width <- function(s) {
    grid::convertWidth(
      grid::grobWidth(grid::textGrob(s, gp = gp)),
      "pt",
      valueOnly = TRUE
    )
  }
  position <- x$at * panel
  if (length(neighbors) > 0) {
    half <- width(x$label) / 2
    room <- vapply(
      neighbor.labels,
      function(s) {
        half + width(s) / 2 + 5.5
      },
      0
    )
    left <- which(neighbors < x$at)
    right <- which(neighbors > x$at)
    left <- if (length(left) > 0) {
      left[which.max(neighbors[left])]
    } else {
      NA
    }
    right <- if (length(right) > 0) {
      right[which.min(neighbors[right])]
    } else {
      NA
    }
    lowest <- if (is.na(left)) {
      -Inf
    } else {
      neighbors[left] * panel + room[left]
    }
    highest <- if (is.na(right)) {
      Inf
    } else {
      neighbors[right] * panel - room[right]
    }
    if (lowest <= highest) {
      # Room between the two neighbors, so move as little as clears both.
      position <- min(max(position, lowest), highest)
    } else {
      # No room between them, so the label goes on a line below theirs, under
      # its tick.
      y <- y - grid::unit(12, "pt")
    }
  }
  children <- grid::gList(
    children,
    grid::textGrob(
      x$label,
      x = grid::unit(position, "pt"),
      y = y,
      vjust = 1,
      gp = gp
    )
  )
  grid::setChildren(x, children)
}

# The frame of panels stacked in a column. Each panel gets an outcome axis line
# as in add_range_frame, the bottom panel a time axis line, each panel with a
# band its own weight axis, and, for a staggered design's cohorts, each panel a
# label of its adoption period under it, the bottom panel's moved if needed to
# clear the time axis's labels.
add_stacked_frames <- function(p, ribbons, weights, vlines) {
  built <- ggplot2::ggplot_build(p)
  layout <- built$layout$layout
  finite <- function(v) {
    v[is.finite(v)]
  }
  bottom.row <- which.max(layout$ROW)
  frames <- list()
  for (i in seq_len(nrow(layout))) {
    in.panel <- lapply(built$data, function(d) {
      d[d$PANEL == layout$PANEL[i], , drop = FALSE]
    })
    ys <- finite(unlist(lapply(in.panel, function(d) {
      c(d$y, d$ymin, d$ymax, d$yend)
    })))
    y <- range(c(ys, finite(built$layout$panel_params[[i]]$y$breaks)))
    frames[[length(frames) + 1]] <- data.frame(
      facet = layout$facet[i],
      x = -Inf,
      y = y,
      group = paste("y", i)
    )
    if (i == bottom.row) {
      xs <- finite(unlist(lapply(in.panel, function(d) {
        c(d$x, d$xend, d$xintercept)
      })))
      x <- range(c(xs, finite(built$layout$panel_params[[i]]$x$breaks)))
      frames[[length(frames) + 1]] <- data.frame(
        facet = layout$facet[i],
        x = x,
        y = -Inf,
        group = "x"
      )
    }
  }
  p <- p +
    ggplot2::geom_path(
      ggplot2::aes(x = x, y = y, group = group),
      data = do.call(rbind, frames),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = .4
    )
  # Weight axes, with ticks chosen as for the single panel's axis, except that
  # labels are kept at least 12% of the panel's height apart rather than 5%,
  # since stacked panels are shorter than a single panel.
  axes <- list()
  for (f in unique(weights$facet)) {
    band <- weights[weights$facet == f, ]
    bottom <- band$band.bottom[1]
    height <- band$band.height[1]
    top <- band$band.top[1]
    largest <- max((band$ymax - bottom) / height * top)
    if (!is.finite(largest) || largest <= 0) {
      next
    }
    i <- which(as.character(layout$facet) == as.character(f))
    span <- diff(built$layout$panel_params[[i]]$y.range)
    n.ticks <- max(1, min(4, floor((largest / top * height) / span / .12)))
    breaks <- pretty(c(0, largest), n = 4)
    breaks <- breaks[breaks <= largest * (1 + 1e-9)]
    if (n.ticks < 4) {
      breaks <- breaks[seq(
        1,
        length(breaks),
        by = max(1, ceiling((length(breaks) - 1) / n.ticks))
      )]
    }
    digits <- if (length(breaks) > 1) {
      max(2, -floor(log10(min(diff(breaks)))))
    } else {
      2
    }
    at <- bottom + height * breaks / top
    axes[[length(axes) + 1]] <- data.frame(
      facet = f,
      y = at,
      ymin = min(at),
      ymax = max(at),
      label = formatC(breaks, format = "f", digits = digits)
    )
  }
  if (length(axes) > 0) {
    p <- p + weight_axis_layer(do.call(rbind, axes))
  }
  if (!is.null(vlines) && nrow(vlines) > 0) {
    labels <- unique(vlines[, c("facet", "adoption")])
    bottom.facet <- as.character(layout$facet[bottom.row])
    labels <- data.frame(
      facet = labels$facet,
      x = labels$adoption,
      label = as.character(labels$facet),
      tick = as.character(labels$facet) == bottom.facet
    )
    # The time axis's labels, which the bottom panel's label is moved to clear.
    x <- built$layout$panel_params[[bottom.row]]$x
    shown <- !is.na(x$breaks)
    p <- p +
      cohort_label_layer(
        labels,
        data.frame(
          x = x$breaks[shown],
          label = as.character(x$get_labels()[shown])
        )
      )
  }
  p
}

# The choice of ticks shared by the plot functions, loose ticks as the Stata
# command chooses them. `padding` is the proportion each axis is padded by
# beyond its ends. Returns choose.ticks(),
# which chooses the ticks of a range, and axis.functions(), which gives an
# axis's limits, breaks and labels functions, remembering each choice so
# that they agree.
tick_machinery <- function(padding = .02) {
  # The algorithm's steps are 1, 5, 2, 2.5, 4 and 3 times a power of ten.
  # A time axis whose periods are whole numbers, such as years, leaves out
  # 2.5, which would put ticks between periods, as at 2002.5, and if a short
  # panel still gets a fractional tick, takes steps of whole numbers instead,
  # so that every tick is a period.
  steps <- function(whole) {
    if (whole) {
      c(1, 5, 2, 4, 3)
    } else {
      c(1, 5, 2, 2.5, 4, 3)
    }
  }
  fractional <- function(at, whole) {
    whole && any(abs(at - round(at)) > 1e-9)
  }
  whole.ticks <- function(low, high, step, loose) {
    step <- max(1, ceiling(step))
    if (loose) {
      seq(floor(low / step) * step, ceiling(high / step) * step, by = step)
    } else {
      seq(ceiling(low / step) * step, floor(high / step) * step, by = step)
    }
  }
  choose.ticks <- function(low, high, whole = FALSE, given = NULL) {
    # Breaks given by the user are used as they are, the axis spanning them and
    # the data.
    if (!is.null(given)) {
      return(list(at = sort(unique(given)), labels = NULL))
    }
    if (!(is.finite(low) && is.finite(high) && high > low)) {
      return(list(at = unique(c(low, high)), labels = NULL))
    }
    Q <- steps(whole)
    at <- scales::extended_breaks(n = 5, Q = Q, only.loose = TRUE)(c(
      low,
      high
    ))
    if (fractional(at, whole)) {
      at <- whole.ticks(low, high, diff(at)[1], TRUE)
    }
    list(at = at, labels = NULL)
  }
  # Each axis's ticks are chosen once and remembered with the range they set,
  # so that the range, the ticks and their labels always come from the same
  # choice. ggplot2 asks for an axis's limits, breaks and labels separately,
  # calls the breaks and labels functions with the range after padding, and
  # in a later build sees the range frame's lines, which span the limits, as
  # data. So each request is matched to a remembered choice by its range,
  # allowing for rounding, and a new choice is made only for a range not
  # seen before. Choices are kept by axis, so that the time and outcome axes
  # never reuse each other's, even for the same range.
  chosen <- new.env(parent = emptyenv())
  chosen$entries <- list()
  same.range <- function(a, b) {
    isTRUE(all.equal(a, b, tolerance = 1e-8, scale = max(1, abs(diff(b)))))
  }
  find.choice <- function(limits, axis) {
    for (entry in chosen$entries) {
      if (entry$axis == axis && same.range(entry$limits, limits)) {
        return(entry)
      }
    }
    NULL
  }
  remember <- function(low, high, axis, whole, given) {
    choice <- choose.ticks(low, high, whole, given)
    choice$limits <- base::range(c(choice$at, low, high))
    choice$axis <- axis
    chosen$entries <- c(chosen$entries, list(choice))
    choice
  }
  # The limits, breaks and labels functions of an axis, whole-number or not,
  # with any breaks the user gave.
  axis.functions <- function(axis, whole, given) {
    choice <- function(range) {
      width <- diff(range) / (1 + 2 * padding)
      unpadded <- range + c(1, -1) * padding * width
      for (limits in list(unpadded, range)) {
        entry <- find.choice(limits, axis)
        if (!is.null(entry)) {
          return(entry)
        }
      }
      remember(unpadded[1], unpadded[2], axis, whole, given)
    }
    list(
      limits = function(range) {
        if (!is.null(find.choice(range, axis))) {
          return(range)
        }
        remember(range[1], range[2], axis, whole, given)$limits
      },
      breaks = function(range) {
        choice(range)$at
      },
      labels = ggplot2::waiver()
    )
  }
  list(choose.ticks = choose.ticks, axis.functions = axis.functions)
}

# Draws the outcome and time axis lines of a single panel as a range frame.
# Each line runs from the lowest to the highest of the plotted data and the
# axis's breaks, so it covers the data even when the outermost breaks fall
# inside it, and stops where the data or the outermost break ends rather than
# at the edge of the panel. The weight axis on the right keeps its own line.
add_range_frame <- function(p) {
  built <- ggplot2::ggplot_build(p)
  scales <- built$layout$panel_params[[1]]
  finite <- function(v) {
    v[is.finite(v)]
  }
  xs <- finite(unlist(lapply(built$data, function(d) {
    c(d$x, d$xend, d$xintercept)
  })))
  ys <- finite(unlist(lapply(built$data, function(d) {
    c(d$y, d$ymin, d$ymax, d$yend)
  })))
  x <- range(c(xs, finite(scales$x$breaks)))
  y <- range(c(ys, finite(scales$y$breaks)))
  p +
    ggplot2::theme(
      axis.line.x.bottom = ggplot2::element_blank(),
      axis.line.y.left = ggplot2::element_blank()
    ) +
    ggplot2::annotate(
      "path",
      x = x,
      y = c(-Inf, -Inf),
      colour = "black",
      linewidth = .4
    ) +
    ggplot2::annotate(
      "path",
      x = c(-Inf, -Inf),
      y = y,
      colour = "black",
      linewidth = .4
    )
}


#' Takes an `sdid()` object and plots the treated and synthetic control
#' trajectories along with the lambda weights conveying which pre-treatment
#' periods were most informative for predicting post-treatment outcomes among
#' the controls.
#'
#' Two plotting styles are available: "minimal" and "effects". Its default,
#' `style = 'minimal'`, plots both trajectories in black, with the treatment
#' line being solid and synthetic line dashed. The time weights are plotted as
#' a gray band with their scale on a right-hand axis, while the last
#' pre-treatment period is plotted a vertical red dashed line. The alternative,
#' `style = 'effects'` adds a diff-in-diff diagram conveying the treatment
#' effect while reducing the trajectory lines opacity.
#'
#' Feeding a staggered treatment `sdid()` object plots cohort panels stacked in
#' a column and are identified by their adoption periods, labeled under each
#' panel at that period. Adding your own x- or y-scale to the returned plot
#' replaces the one this function sets, and with it the time-weight axis or the
#' x padding. Adding facets or a coordinate system misplaces the time-weight
#' axes of stacked panels. If a list of estimates is passed, plots each of them
#' in stacked in its own panel. Plotting several estimates in one shared panel
#' is not supported.
#'
#' Requires ggplot2 3.5.0.
#'
#' @param estimates The result of an sdid object or a list of sdid estimates.
#'        For a staggered sdid result, each adoption period is plotted in its
#'        a stacked panel.
#' @param treated.name A string conveying the name of the treated curve that
#'        appears in the legend. Defaults to "treated".
#' @param control.name A string conveying the name of the control curve that
#'        appears in the legend. Defaults to 'synthetic control'.
#' @param spaghetti.units A list of units to plot individually alongside the
#'        treatment and control curves. Defaults to the empty list.
#' @param spaghetti.matrices A list of matrices, one for each element of
#'        estimates, of trajectories to plot individually. The rows of these
#'        matrices should be the same length as the trajectories in Y and they
#'        must be named so trajectories can be labeled in the plot.
#' @param lambda.comparable Logical argument where if TRUE, the lambda
#'        weight axis is scaled to range from 0 to 1. Useful for when comparing
#'        across stacked panels. Defaults to FALSE.
#' @param overlay A number between 0 and 1 that shifts the synthetic control
#'        toward the treated trajectory by that fraction of the difference in
#'        their time-weighted pre-treatment averages. At 0, the default, both
#'        are drawn as estimated. Values near 1 bring them close, which makes
#'        it easier to judge whether they move in parallel, and at 1 they
#'        overlap and the effects diagram is left out, as for a synthetic
#'        control estimate. An estimate's own value can be set as its
#'        attribute "overlay". A vector draws the plot at each value, marked
#'        by the frame aesthetic, which plotly shows as an animation.
#' @param lambda.plot.scale A number setting the height of the time weights'
#'        band, the range of the trajectories divided by this number.
#'        Defaults to 3.
#' @param effect.curvature A number setting the curvature of the arrow that
#'        marks the effect in the effects style. Defaults to 0.3.
#' @param line.width A number setting the width of the trajectory lines.
#'        Defaults to 0.5.
#' @param point.size A number setting the size of the points in the effects
#'        style's diagram. Defaults to 1.
#' @param onset.alpha A number between 0 and 1 setting the opacity of the red
#'        dashed line marking the start of treatment, drawn midway between
#'        the last period before treatment and the first treated period.
#'        Defaults to 0.8.
#' @param spaghetti.line.width A number setting the width of the trajectories
#'        drawn by spaghetti.units and spaghetti.matrices. Defaults to 0.2.
#' @param spaghetti.label.size A number setting the size of their labels.
#'        Defaults to 2.
#' @param spaghetti.line.alpha A number between 0 and 1 setting the opacity
#'        of those trajectories. Defaults to 0.3.
#' @param spaghetti.label.alpha A number between 0 and 1 setting the opacity
#'        of their labels. Defaults to 0.5.
#' @param alpha.multiplier A vector with one number per estimate, multiplying
#'        the opacity of everything drawn for that estimate, to emphasize
#'        some estimates over others. Defaults to 1 for every estimate.
#' @param style The plot's style, "minimal" (the default) or "effects". See
#'        Details.
#' @param x.breaks,y.breaks Optional ticks for the time and outcome axes, as
#'        numbers, or as dates for a time axis of dates. A numeric axis runs
#'        between its outermost ticks, widened to cover the data if needed,
#'        and the frame, padding, cohort labels and time weights' band follow
#'        them. The band's baseline sits on the highest outcome tick below
#'        it, if there is one, and with stacked panels the outcome ticks
#'        apply to every panel. A time axis of dates uses the given dates as
#'        ticks and keeps the range of the data. Prefer these to adding a
#'        scale to the returned plot, which replaces the one it sets.
#' @param lambda.axis Logical argument where if TRUE, each panel with a time
#'        weights' band shows the band's scale on a right-hand axis. Defaults
#'        to TRUE.
#' @return A ggplot object, which draws the plot when printed and can be
#'   extended with `+`, for example with labels, themes or further layers.
#'   Scales, facets and coordinate systems added this way replace those the plot
#'   sets, as described above.
#' @export sdid_plot
sdid_plot <- function(
  estimates,
  treated.name = "treated",
  control.name = "synthetic control",
  spaghetti.units = c(),
  spaghetti.matrices = NULL,
  lambda.comparable = FALSE,
  overlay = 0,
  lambda.plot.scale = 3,
  effect.curvature = .3,
  line.width = .5,
  point.size = 1,
  onset.alpha = .8,
  spaghetti.line.width = .2,
  spaghetti.label.size = 2,
  spaghetti.line.alpha = .3,
  spaghetti.label.alpha = .5,
  alpha.multiplier = NULL,
  style = c("minimal", "effects"),
  lambda.axis = TRUE,
  x.breaks = NULL,
  y.breaks = NULL
) {
  estimates <- plot_estimates(estimates)
  cohorts <- isTRUE(attr(estimates, "sdid.cohorts"))
  require_ggplot2()
  style <- match.arg(style)
  for (breaks in list(x.breaks = x.breaks, y.breaks = y.breaks)) {
    if (
      !is.null(breaks) &&
        (!(is.numeric(breaks) || inherits(breaks, "Date")) ||
          length(breaks) == 0 ||
          any(!is.finite(as.numeric(breaks))))
    ) {
      stop(
        "x.breaks and y.breaks must be finite numbers, or dates for a time axis of dates."
      )
    }
  }
  if (inherits(y.breaks, "Date")) {
    stop("y.breaks must be numbers.")
  }
  if (!is.null(x.breaks)) {
    x.breaks <- as.numeric(x.breaks)
  }
  if (inherits(estimates, "sdid_estimate")) {
    estimates <- list(estimates)
  }
  if (is.null(names(estimates))) {
    names(estimates) <- sprintf("estimate %d", 1:length(estimates))
  }
  if (is.null(alpha.multiplier)) {
    alpha.multiplier <- rep(1, length(estimates))
  }
  if (
    !is.null(spaghetti.matrices) &&
      length(spaghetti.matrices) != length(estimates)
  ) {
    stop("spaghetti.matrices must be the same length as estimates")
  }
  multiple.frames <- length(overlay) > 1
  treated <- 1
  control <- 2
  groups <- factor(c(control, treated), labels = c(control.name, treated.name))
  estimate.factors <- factor(
    1:(length(estimates) + 1),
    labels = c(treated.name, names(estimates))
  )
  facet_factors <- factor(seq_along(estimates), labels = names(estimates))
  grid <- expand.grid(
    estimate = 1:length(estimates),
    overlay = 1:length(overlay)
  )
  plot.descriptions <- lapply(1:nrow(grid), function(row) {
    est <- estimates[[grid$estimate[row]]]
    over <- overlay[grid$overlay[row]]
    setup <- attr(est, "setup")
    weights <- attr(est, "weights")
    Y <- setup$Y - contract3(setup$X, weights$beta)
    N0 <- setup$N0
    N1 <- nrow(Y) - N0
    T0 <- setup$T0
    T1 <- ncol(Y) - T0

    lambda.synth <- c(weights$lambda, rep(0, T1))
    lambda.target <- c(rep(0, T0), rep(1 / T1, T1))
    omega.synth <- c(weights$omega, rep(0, N1))
    omega.target <- c(rep(0, N0), rep(1 / N1, N1))

    # pull estimate-specific overlay from attribute if present if we're given a
    # synthetic control estimate or overlay is one, take note: we'll plot it
    # differently
    if (!is.null(attr(est, "overlay"))) {
      over <- attr(est, "overlay")
    }
    is.sc <- all(weights$lambda == 0) || over == 1

    intercept.offset <- over *
      c((omega.target - omega.synth) %*% Y %*% lambda.synth)
    obs.trajectory <- as.numeric(omega.target %*% Y)
    syn.trajectory <- as.numeric(omega.synth %*% Y) + intercept.offset
    spaghetti.trajectories <- Y[
      rownames(Y) %in% spaghetti.units,
      ,
      drop = FALSE
    ]
    if (!is.null(spaghetti.matrices)) {
      more.spaghetti.trajectories <- spaghetti.matrices[[grid$estimate[row]]]
      if (ncol(more.spaghetti.trajectories) != ncol(Y)) {
        stop(
          "The elements of spaghetti.matrices must be matrices with the same number of columns as Y"
        )
      }
      if (is.null(rownames(more.spaghetti.trajectories))) {
        stop("The elements of the list spaghetti.matrices must have named rows")
      }
      spaghetti.trajectories <- rbind(
        spaghetti.trajectories,
        more.spaghetti.trajectories
      )
    }

    treated.post <- omega.target %*% Y %*% lambda.target
    treated.pre <- omega.target %*% Y %*% lambda.synth
    control.post <- omega.synth %*% Y %*% lambda.target + intercept.offset
    control.pre <- omega.synth %*% Y %*% lambda.synth + intercept.offset
    sdid.post <- as.numeric(control.post + treated.pre - control.pre)

    time <- as.numeric(timesteps(Y))
    if (length(time) == 0 || !all(is.finite(time))) {
      time <- 1:(T0 + T1)
    }
    pre.time <- lambda.synth %*% time
    post.time <- lambda.target %*% time

    # construct objects on graph
    lines <- data.frame(
      x = rep(time, 2),
      y = c(obs.trajectory, syn.trajectory),
      color = rep(groups[c(treated, control)], each = length(time))
    )
    points <- data.frame(
      x = c(post.time, post.time),
      y = c(treated.post, sdid.post),
      color = groups[c(treated, control)]
    )
    did.points <- data.frame(
      x = c(pre.time, pre.time, post.time, post.time),
      y = c(treated.pre, control.pre, control.post, treated.post),
      color = groups[c(treated, control, control, treated)]
    )
    did.segments <- data.frame(
      x = c(pre.time, pre.time),
      xend = c(post.time, post.time),
      y = c(control.pre, treated.pre),
      yend = c(control.post, treated.post),
      color = groups[c(control, treated)]
    )
    hallucinated.segments <- data.frame(
      x = pre.time,
      xend = post.time,
      y = treated.pre,
      yend = sdid.post
    )
    arrows <- data.frame(
      x = post.time,
      xend = post.time,
      y = sdid.post,
      yend = treated.post,
      xscale = max(time) - post.time,
      color = groups[control]
    )
    spaghetti.lines <- data.frame(
      x = rep(time, nrow(spaghetti.trajectories)),
      y = as.vector(t(spaghetti.trajectories)),
      unit = rep(rownames(spaghetti.trajectories), each = length(time))
    )
    spaghetti.labels <- data.frame(
      x = rep(time[1], nrow(spaghetti.trajectories)),
      y = as.vector(spaghetti.trajectories[, 1]),
      unit = rownames(spaghetti.trajectories)
    )

    # The treatment line is drawn midway between the last period before
    # treatment and the first period treated, the adoption period, which is
    # recorded with it.
    T0s <- attr(est, "T0s")
    if (is.null(T0s)) {
      T0s <- T0
    }
    vlines <- data.frame(
      xintercept = (time[T0s] + time[T0s + 1]) / 2,
      adoption = time[T0s + 1]
    )

    if (lambda.comparable) {
      height <- (max(c(obs.trajectory)) - min(c(obs.trajectory))) /
        lambda.plot.scale
      bottom <- min(c(obs.trajectory)) - height
      ribbons <- data.frame(
        x = time[1:T0],
        ymin = rep(bottom, T0),
        ymax = bottom + height * lambda.synth[1:T0],
        color = groups[control],
        band.bottom = bottom,
        band.height = height,
        band.top = 1
      )
    } else {
      height <- (max(c(obs.trajectory, syn.trajectory)) -
        min(c(obs.trajectory, syn.trajectory))) /
        lambda.plot.scale
      bottom <- min(c(obs.trajectory, syn.trajectory)) - height
      ribbons <- data.frame(
        x = time[1:T0],
        ymin = rep(bottom, T0),
        ymax = bottom + height * lambda.synth[1:T0] / max(lambda.synth),
        color = groups[control],
        band.bottom = bottom,
        band.height = height,
        band.top = max(lambda.synth)
      )
    }
    elements <- list(
      lines = lines,
      points = points,
      did.segments = did.segments,
      did.points = did.points,
      hallucinated.segments = hallucinated.segments,
      arrows = arrows,
      spaghetti.lines = spaghetti.lines,
      spaghetti.labels = spaghetti.labels,
      vlines = vlines,
      ribbons = ribbons
    )
    lapply(elements, function(x) {
      if (nrow(x) > 0) {
        x$frame <- over
        x$is.sc <- is.sc
        # offset because the treated pseudo-estimate factor is first
        x$estimate <- estimate.factors[grid$estimate[row] + 1]
      }
      x
    })
  })

  # Each estimate has its own panel, labeled by its name, so names must not
  # repeat.
  if (anyDuplicated(names(estimates))) {
    stop(
      "Each estimate is plotted in its own panel, labeled by its name, ",
      "so the names of the list must not repeat."
    )
  }
  concatenate.field <- function(field) {
    do.call(
      rbind,
      lapply(plot.descriptions, function(desc) {
        element <- desc[[field]]
        estimate.factor <- element$estimate[1]
        # offset because the treated pseudo-estimate factor is first
        element$facet <- facet_factors[as.integer(estimate.factor) - 1]
        element$show <- alpha.multiplier[as.integer(element$estimate) - 1] # "
        # show treated observations
        element$show[element$color == groups[treated]] <- 1
        if ("color" %in% colnames(element)) {
          element$role <- factor(
            ifelse(
              element$color == groups[treated],
              treated.name,
              control.name
            ),
            levels = c(treated.name, control.name)
          )
        }
        element
      })
    )
  }
  conc <- lapply(names(plot.descriptions[[1]]), concatenate.field)
  names(conc) <- names(plot.descriptions[[1]])
  no.sc <- function(x) {
    x[!x$is.sc, ]
  }

  # invoke the geom with 'frame=frame' included in the aesthetic only if there
  # are multiple frames This cuts down on warnings without restricting function.
  # ggplotly understands frame and uses it to create animations if there are
  # multiple ggplot2's display doesn't and warns 'ignoring unknown aesthetic' if
  # you include it returns geom(aes, data=data, ...) where aes=base.aes if
  # there's a single frame
  #					    aes=(base.aes, frame=frame) if there's more than one
  with.frame <- function(geom, base.aes, data, ...) {
    new.aes <- if (multiple.frames) {
      modifyList(base.aes, ggplot2::aes(frame = frame))
    } else {
      base.aes
    }
    do.call(geom, c(list(new.aes, data = data), list(...)))
  }

  # Opacities are actual opacities, so alpha uses an identity scale. The effects
  # style draws the trajectories at half opacity and adds the diff-in-diff
  # diagram.
  ribbon.alpha <- .34
  trajectory.opacity <- if (style == "effects") {
    .5
  } else {
    1
  }
  ribbons <- no.sc(conc$ribbons)
  # Ticks on the numeric axes. Each axis is labeled as the Stata command labels
  # its axes by default, with about five nice ticks covering the data, and its
  # range is widened to its outermost ticks. The ticks come from the algorithm
  # ggplot2 uses by default (Talbot, Lin and Hanrahan 2010) in its loose form,
  # which places the first tick at or below the smallest value and the last at
  # or above the largest. Both axes are padded by 2% of their range beyond their
  # ends, so that the content sits close to the axis lines.
  padding <- .02
  machinery <- tick_machinery(padding)
  choose.ticks <- machinery$choose.ticks
  axis.functions <- machinery$axis.functions
  y.axis <- axis.functions("y", FALSE, y.breaks)
  # The time axis is a whole-number axis when every period is a whole number.
  x.axis <- axis.functions(
    "x",
    all(abs(conc$lines$x - round(conc$lines$x)) < 1e-9),
    x.breaks
  )
  # Move each panel's band down onto the nearest tick below it, so that a weight
  # of zero lines up with a labeled value on the left axis. Moving the band
  # changes the panel's range, which can change the ticks, so the step repeats
  # until the band sits on one of the ticks of its final range.
  for (f in unique(ribbons$facet)) {
    in.panel <- ribbons$facet == f
    bottom <- min(ribbons$ymin[in.panel])
    top <- max(c(
      conc$lines$y[conc$lines$facet == f],
      conc$spaghetti.lines$y[conc$spaghetti.lines$facet == f],
      ribbons$ymax[in.panel]
    ))
    target <- bottom
    for (attempt in 1:5) {
      breaks <- choose.ticks(target, top, FALSE, y.breaks)$at
      if (any(abs(breaks - target) <= 1e-9 * (top - target))) {
        break
      }
      below <- breaks[breaks < target]
      target <- if (length(below) > 0) {
        max(below)
      } else {
        floor(target / diff(breaks)[1]) * diff(breaks)[1]
      }
    }
    aligned <- any(
      abs(choose.ticks(target, top, FALSE, y.breaks)$at - target) <=
        1e-9 * (top - target)
    )
    if (aligned) {
      shift <- bottom - target
      ribbons$ymin[in.panel] <- ribbons$ymin[in.panel] - shift
      ribbons$ymax[in.panel] <- ribbons$ymax[in.panel] - shift
      ribbons$band.bottom[in.panel] <- ribbons$band.bottom[in.panel] - shift
    }
  }
  # Keep only the stretches of each band where the weights are positive, with
  # the zero-weight period on either side so each stretch rises from the
  # baseline. Each stretch is its own group, so nothing is drawn between them.
  ribbons <- do.call(
    rbind,
    lapply(
      split(ribbons, list(ribbons$facet, ribbons$estimate), drop = TRUE),
      function(band) {
        band <- band[order(band$x), ]
        positive <- band$ymax - band$ymin > 1e-12 * band$band.height
        keep <- positive |
          c(positive[-1], FALSE) |
          c(FALSE, positive[-length(positive)])
        band$stretch <- paste(
          band$estimate,
          cumsum(keep & !c(FALSE, keep[-length(keep)]))
        )
        band[keep, ]
      }
    )
  )
  if (is.null(ribbons)) {
    # No band at all, as for an SC estimate, whose time weights are all zero.
    ribbons <- no.sc(conc$ribbons)[0, ]
    ribbons$stretch <- character(0)
  }
  # Each positive stretch of the band is outlined all the way around, so the
  # outline also runs beneath it.
  p <- ggplot2::ggplot() +
    ggplot2::geom_ribbon(
      ggplot2::aes(
        x = x,
        ymin = ymin,
        ymax = ymax,
        group = stretch,
        alpha = ribbon.alpha * show
      ),
      data = ribbons,
      fill = "black",
      color = "black",
      linewidth = line.width,
      outline.type = "full",
      show.legend = FALSE
    ) +
    with.frame(
      ggplot2::geom_line,
      ggplot2::aes(
        x = x,
        y = y,
        group = color,
        linetype = role,
        alpha = trajectory.opacity * show
      ),
      data = conc$lines,
      color = "black",
      linewidth = line.width
    )
  if (style == "effects") {
    # The diff-in-diff diagram in black. The treated and synthetic control trend
    # segments are solid, the treated unit's counterfactual, which runs parallel
    # to the synthetic control's trend, is dotted, and a red arrow shows the
    # effect.
    p <- p +
      with.frame(
        ggplot2::geom_segment,
        ggplot2::aes(
          x = x,
          xend = xend,
          y = y,
          yend = yend,
          group = color,
          alpha = show
        ),
        data = no.sc(conc$did.segments),
        color = "black",
        linetype = "solid",
        linewidth = line.width
      ) +
      with.frame(
        ggplot2::geom_segment,
        ggplot2::aes(
          x = x,
          xend = xend,
          y = y,
          yend = yend,
          group = estimate,
          alpha = show
        ),
        data = no.sc(conc$hallucinated.segments),
        color = "black",
        linetype = "dotted",
        linewidth = line.width
      ) +
      with.frame(
        ggplot2::geom_point,
        ggplot2::aes(x = x, y = y, alpha = show),
        data = no.sc(conc$did.points),
        color = "black",
        size = point.size
      ) +
      with.frame(
        ggplot2::geom_point,
        ggplot2::aes(x = x, y = y, alpha = show),
        data = conc$points,
        color = "black",
        shape = 21,
        size = point.size
      ) +
      ggplot2::geom_curve(
        ggplot2::aes(x = x, xend = xend, y = y, yend = yend, alpha = show),
        data = conc$arrows,
        curvature = effect.curvature,
        color = "red",
        linewidth = line.width,
        arrow = ggplot2::arrow(length = ggplot2::unit(.2, "cm"))
      )
  }
  p <- p +
    ggplot2::geom_vline(
      ggplot2::aes(xintercept = xintercept, alpha = onset.alpha * show),
      data = conc$vlines,
      linewidth = line.width,
      color = "red",
      linetype = "dashed"
    )
  if (nrow(conc$spaghetti.labels) > 0) {
    p <- p +
      ggplot2::geom_text(
        ggplot2::aes(
          x = x,
          y = y,
          label = unit,
          alpha = spaghetti.label.alpha * show
        ),
        data = conc$spaghetti.labels,
        color = "black",
        size = spaghetti.label.size
      ) +
      ggplot2::geom_line(
        ggplot2::aes(
          x = x,
          y = y,
          group = unit,
          alpha = spaghetti.line.alpha * show
        ),
        data = conc$spaghetti.lines,
        color = "black",
        linewidth = spaghetti.line.width
      )
  }
  single.panel <- all(conc$lines$facet == conc$lines$facet[1])
  stacked <- !single.panel
  if (stacked) {
    # Panels are stacked in a column, in the order the estimates were given.
    facet.levels <- unique(as.character(conc$lines$facet))
    p <- p +
      ggplot2::facet_wrap(
        ggplot2::vars(facet = factor(facet, levels = facet.levels)),
        ncol = 1,
        scales = "free_y"
      )
  }
  p <- p +
    tryCatch(
      {
        as.Date(colnames(attr(estimates[[1]], "setup")$Y))
        ggplot2::scale_x_continuous(
          labels = function(time) {
            # Dates are stored as days since R's origin, 1970-01-01, so the
            # axis is numeric and its labels convert the days back to dates.
            # R before 4.3.0 requires the origin to be given.
            as.Date(time, origin = "1970-01-01")
          },
          breaks = if (is.null(x.breaks)) {
            ggplot2::waiver()
          } else {
            x.breaks
          },
          expand = ggplot2::expansion(mult = padding)
        )
      },
      error = function(e) {
        # Date axes keep ggplot2's ticks, since nice numbers of days are not
        # nice dates.
        ggplot2::scale_x_continuous(
          limits = x.axis$limits,
          breaks = x.axis$breaks,
          labels = x.axis$labels,
          expand = ggplot2::expansion(mult = padding)
        )
      }
    )

  # Right-hand axis for the time weights. A ribbon reaching height h above the
  # band's bottom represents the weight h / band.height * band.top. It is drawn
  # only for a single band in a single panel, since one secondary axis cannot
  # describe several bands drawn on different scales.
  bands <- unique(ribbons[,
    c("band.bottom", "band.height", "band.top"),
    drop = FALSE
  ])
  weight.title.y <- NULL
  if (lambda.axis && single.panel && nrow(bands) == 1) {
    band.bottom <- bands$band.bottom
    band.height <- bands$band.height
    band.top <- bands$band.top
    largest <- max((ribbons$ymax - band.bottom) / band.height * band.top)
    if (is.finite(largest) && largest > 0) {
      # Fewer ticks when the weights occupy little of the panel, keeping labels
      # about 5% of its height apart.
      y.span <- diff(range(c(conc$lines$y, ribbons$ymin, ribbons$ymax)))
      occupied <- (largest / band.top * band.height) / y.span
      n.ticks <- max(1, min(4, floor(occupied / .05)))
      breaks <- pretty(c(0, largest), n = 4)
      breaks <- breaks[breaks <= largest * (1 + 1e-9)]
      if (n.ticks < 4) {
        step <- max(1, ceiling((length(breaks) - 1) / n.ticks))
        breaks <- breaks[seq(1, length(breaks), by = step)]
      }
      digits <- if (length(breaks) > 1) {
        max(2, -floor(log10(min(diff(breaks)))))
      } else {
        2
      }
      p <- p +
        ggplot2::scale_y_continuous(
          limits = y.axis$limits,
          breaks = y.axis$breaks,
          labels = y.axis$labels,
          expand = ggplot2::expansion(mult = padding),
          sec.axis = ggplot2::sec_axis(
            ~ (. - band.bottom) / band.height * band.top,
            name = NULL,
            breaks = breaks,
            labels = function(w) {
              formatC(w, format = "f", digits = digits)
            }
          )
        )
      weight.title.y <- band.bottom + band.height * max(breaks) / band.top
    }
  }
  # Plots without a weight axis still get the same vertical padding.
  if (is.null(weight.title.y)) {
    p <- p +
      ggplot2::scale_y_continuous(
        limits = y.axis$limits,
        breaks = y.axis$breaks,
        labels = y.axis$labels,
        expand = ggplot2::expansion(mult = padding)
      )
  }

  p <- p +
    ggplot2::scale_linetype_manual(
      values = stats::setNames(
        c("solid", "dashed"),
        c(treated.name, control.name)
      ),
      name = NULL
    ) +
    ggplot2::xlab("") +
    ggplot2::ylab("") +
    ggplot2::scale_alpha_identity() +
    ggplot2::theme_bw() +
    ggplot2::theme(
      legend.direction = "horizontal",
      legend.position = "top",
      legend.key.width = ggplot2::unit(1.5, "cm"),
      panel.grid.major = ggplot2::element_line(
        color = "lightgray",
        linewidth = .6,
        linetype = "dotted"
      ),
      panel.grid.minor = ggplot2::element_blank(),
      axis.ticks.length = ggplot2::unit(.25, "cm"),
      axis.text = ggplot2::element_text(size = 10),
      legend.text = ggplot2::element_text(size = 12)
    )
  # A single panel gets an open frame. Its outcome and time axis lines are drawn
  # as a range frame once the plot is complete, below. The weight axis stops at
  # its outermost ticks, as it describes only the band. Panels stacked in a
  # column each get the single panel's open frame, drawn once the plot is
  # complete, below, with the time axis only under the bottom panel, so that the
  # gridlines run through. A staggered design's cohorts are identified by their
  # adoption periods, labeled under each panel, and other panels by a title
  # above them. Panels side by side keep the box around each panel.
  if (stacked) {
    p <- p +
      ggplot2::theme(
        panel.border = ggplot2::element_blank(),
        axis.line = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_line(color = "black"),
        strip.background = ggplot2::element_blank(),
        strip.text = if (cohorts) {
          ggplot2::element_blank()
        } else {
          ggplot2::element_text(
            size = 11,
            colour = "black",
            hjust = 0,
            margin = ggplot2::margin(4, 0, 4, 0)
          )
        },
        panel.spacing = ggplot2::unit(
          if (cohorts) {
            .9
          } else {
            .4
          },
          "cm"
        )
      ) +
      ggplot2::coord_cartesian(clip = "off")
    # Room on the right for the weight axes, and below for a cohort label that
    # must go under the time axis's labels.
    right <- if (lambda.axis && nrow(ribbons) > 0) {
      5.5 + 48
    } else {
      5.5
    }
    bottom <- if (cohorts) {
      5.5 + 12
    } else {
      5.5
    }
    p <- p +
      ggplot2::theme(plot.margin = ggplot2::margin(5.5, right, bottom, 5.5))
  }
  if (single.panel) {
    p <- p +
      ggplot2::theme(
        panel.border = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(color = "black", linewidth = .4),
        axis.ticks = ggplot2::element_line(color = "black")
      )
    # Only when the weight axis exists. A guide for a secondary axis that does
    # not exist makes ggplot2 3.5 and later draw a copy of the left axis on the
    # right.
    if (
      utils::packageVersion("ggplot2") >= "3.5.0" && !is.null(weight.title.y)
    ) {
      p <- p + ggplot2::guides(y.sec = ggplot2::guide_axis(cap = "both"))
    }
  }
  # The weight axis title runs vertically beside the axis labels, as an axis
  # title would, but is centered on the span of the weight axis (zero to its top
  # tick) rather than on the whole panel. It is placed in data coordinates at
  # the panel's right edge, offset past the ticks and labels, so the panel stops
  # clipping and the right margin widens to hold it. Widths are left to grid to
  # measure when the plot is drawn, on the device in use, rather than measured
  # here, which would open a graphics device as a side effect.
  if (!is.null(weight.title.y)) {
    label.widths <- lapply(
      formatC(breaks, format = "f", digits = digits),
      function(label) {
        grid::grobWidth(grid::textGrob(label, gp = grid::gpar(fontsize = 10)))
      }
    )
    gap <- grid::unit(5.5, "pt")
    offset <- grid::unit(.25, "cm") +
      grid::unit(2.2, "pt") +
      do.call(grid::unit.pmax, label.widths) +
      gap
    title.grob <- grid::textGrob(
      expression("Time weight" ~ (lambda)),
      x = grid::unit(0, "npc") + offset,
      y = grid::unit(.5, "npc"),
      rot = -90,
      hjust = .5,
      vjust = 0,
      gp = grid::gpar(fontsize = 11)
    )
    p <- p +
      ggplot2::annotation_custom(
        title.grob,
        xmin = Inf,
        xmax = Inf,
        ymin = band.bottom,
        ymax = weight.title.y
      ) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::theme(
        plot.margin = grid::unit.c(
          grid::unit(5.5, "pt"),
          grid::unit(5.5, "pt") + gap + grid::grobWidth(title.grob),
          grid::unit(5.5, "pt"),
          grid::unit(5.5, "pt")
        )
      )
  }
  if (single.panel) {
    p <- add_range_frame(p)
  }
  if (stacked) {
    p <- add_stacked_frames(
      p,
      ribbons,
      if (lambda.axis) {
        ribbons
      } else {
        ribbons[0, ]
      },
      if (cohorts) {
        conc$vlines
      } else {
        NULL
      }
    )
  }
  p
}

#' Takes an `sdid()` object and plots the treated and synthetic control
#' trajectories along with a placebo plot stacked in a column. Requires ggplot2.
#' @param estimate The result of `sdid()` with a single adoption period.
#' @param treated.fraction The share of the periods before treatment that
#'        the placebo treats, as in `sdid_placebo()`. Defaults to NULL.
#' @param style The plot's style, "minimal" (the default) or "effects", as
#'        in `sdid_plot()`.
#' @param placebo.time The period at which the placebo treatment begins, as
#'        in `sdid_placebo()`. Defaults to "auto".
#' @param x.breaks,y.breaks Optional ticks for the time and outcome axes, as
#'        in `sdid_plot()`. The outcome ticks apply to both panels.
#' @return A ggplot object, which draws the plot when printed and can be
#'   extended with `+`, for example with labels, themes or further layers.
#'   Scales, facets and coordinate systems added this way replace those the plot
#'   sets, as for [sdid_plot()].
#' @export sdid_placebo_plot
sdid_placebo_plot <- function(
  estimate,
  treated.fraction = NULL,
  style = c("minimal", "effects"),
  placebo.time = "auto",
  x.breaks = NULL,
  y.breaks = NULL
) {
  require_ggplot2()
  style <- match.arg(style)
  estimate <- placebo_estimate(estimate)
  estimates <- list(
    estimate = estimate,
    placebo = sdid_placebo(
      estimate,
      treated.fraction = treated.fraction,
      placebo.time = placebo.time
    )
  )
  sdid_plot(
    estimates,
    style = style,
    x.breaks = x.breaks,
    y.breaks = y.breaks
  )
}

#' Plots each control unit's weight and its own estimate
#'
#' Synthetic difference-in-differences estimates the effect as a weighted
#' average, over the control units, of the difference-in-differences comparison
#' of the treated units with each control alone, using the time weights. This
#' plot shows those comparisons with their unit weights, for an estimate or each
#' cohort of a staggered design.
#'
#' `type = "scatter"` (the default) plots each control's weight against its own
#' estimate, with the effect as a red line, so that the controls that make up
#' the synthetic control are at the right and whether they agree with each other
#' is seen in their spread about the line. Every control is drawn, those with no
#' weight as open circles, so that the plot also shows whether the weighting
#' changed the answer, and `label.units` chooses which are named. `type =
#' "bars"` draws a bar for each control with weight chosen by `label.units`,
#' ordered from the largest weight down, running from the effect to the
#' control's own estimate, with its weight at the end. A control pulls the
#' estimate up if its bar runs right of the line and down if left, by its weight
#' times the bar's length.
#'
#' Requires ggplot2, and the ggrepel package to place the scatter's labels
#' without overlap. Without ggrepel, labels that would overlap are left out.
#' @param estimates The result of `sdid()`, or a named list of results with
#'        a single adoption period each, drawn in panels labeled by their
#'        names. A staggered result is drawn with a panel for each cohort,
#'        unless `cohort` chooses one.
#' @param type The kind of plot. "scatter" (the default) plots each
#'        control's weight against its own estimate. "bars" draws a bar from
#'        the effect to each weighted control's own estimate, ordered by
#'        weight. See Details.
#' @param cohort For a staggered result, the adoption period of the one
#'        cohort to plot, as it appears in the result's `cohorts`. Defaults
#'        to NULL, plotting every cohort in its own panel.
#' @param label.units Which controls to name. "all" names every control,
#'        "weighted" those with positive weight, and a whole number n the n
#'        with the largest weights, with any tied at the cutoff included.
#'        Defaults to "all" for a single panel and 10 for several. The
#'        scatter still draws the controls it leaves unnamed, while the bars
#'        draw only the controls they name.
#' @param weighted.only TRUE to leave controls with no weight out of the
#'        scatter. Defaults to FALSE, drawing them as open circles.
#' @return A ggplot object, which draws the plot when printed and can be
#'   extended with `+`, for example with labels, themes or further layers.
#'   Scales, facets and coordinate systems added this way replace those the plot
#'   sets.
#' @export sdid_units_plot
sdid_units_plot <- function(
  estimates,
  type = c("scatter", "bars"),
  cohort = NULL,
  label.units = NULL,
  weighted.only = FALSE
) {
  require_ggplot2()
  type <- match.arg(type)
  labelled.given <- !is.null(label.units)
  if (!is.null(cohort)) {
    if (!inherits(estimates, "sdid") || estimates$design != "staggered") {
      stop(
        "cohort chooses a cohort of the result of sdid for a staggered design."
      )
    }
    cohort <- as.character(cohort)
    if (length(cohort) != 1 || !cohort %in% names(estimates$estimates)) {
      stop(sprintf(
        "cohort must be one of the adoption periods, %s.",
        paste(names(estimates$estimates), collapse = ", ")
      ))
    }
    estimates <- estimates$estimates[cohort]
  } else {
    estimates <- plot_estimates(estimates)
  }
  if (inherits(estimates, "sdid_estimate")) {
    estimates <- list(estimates)
  }
  if (is.null(names(estimates))) {
    names(estimates) <- sprintf("estimate %d", seq_along(estimates))
  }
  single <- length(estimates) == 1
  if (!labelled.given) {
    label.units <- if (single) {
      "all"
    } else {
      10
    }
  }
  valid <- (is.character(label.units) &&
    length(label.units) == 1 &&
    label.units %in% c("all", "weighted")) ||
    (is.numeric(label.units) &&
      length(label.units) == 1 &&
      is.finite(label.units) &&
      label.units >= 0 &&
      label.units == round(label.units))
  if (!valid) {
    stop('label.units must be "all", "weighted" or a whole number of units.')
  }

  # Each control's own estimate, the difference-in-differences of the treated
  # units with that control alone, whose average with the unit weights is the
  # estimate.
  data <- do.call(
    rbind,
    lapply(seq_along(estimates), function(e) {
      estimate <- estimates[[e]]
      setup <- attr(estimate, "setup")
      weights <- attr(estimate, "weights")
      Y <- setup$Y - contract3(setup$X, weights$beta)
      N0 <- setup$N0
      N1 <- nrow(Y) - N0
      T0 <- setup$T0
      T1 <- ncol(Y) - T0
      change <- c(rep(0, T0), rep(1 / T1, T1)) - c(weights$lambda, rep(0, T1))
      own <- sum(colMeans(Y[(N0 + 1):nrow(Y), , drop = FALSE]) * change) -
        as.vector(Y[1:N0, , drop = FALSE] %*% change)
      data.frame(
        unit = rownames(Y)[1:N0],
        weight = weights$omega,
        own = own,
        effect = c(estimate),
        panel = names(estimates)[e],
        stringsAsFactors = FALSE
      )
    })
  )
  data$panel <- factor(data$panel, levels = names(estimates))
  effects <- data[!duplicated(data$panel), c("panel", "effect")]
  effects$label <- sprintf(
    "ATT %s",
    formatC(effects$effect, format = "f", digits = 1)
  )

  padding <- .02
  machinery <- tick_machinery(padding)
  look <- ggplot2::theme_bw() +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(
        colour = "grey85",
        linetype = "dotted"
      ),
      axis.ticks = ggplot2::element_line(colour = "black"),
      axis.ticks.length = ggplot2::unit(.25, "cm"),
      axis.text = ggplot2::element_text(size = 10),
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 12),
      legend.key = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        size = 11,
        colour = "black",
        hjust = 0
      ),
      panel.spacing = ggplot2::unit(.6, "cm")
    )
  facet <- if (single) {
    NULL
  } else {
    ggplot2::facet_wrap(ggplot2::vars(panel), ncol = 1, scales = "free")
  }

  # The controls chosen by label.units in each panel, which the scatter names
  # and the bars draw.
  data$named <- unlist(
    lapply(split(data$weight, data$panel), function(w) {
      positive <- w > 0
      if (identical(label.units, "all")) {
        return(rep(TRUE, length(w)))
      }
      if (identical(label.units, "weighted") || sum(positive) <= label.units) {
        return(positive)
      }
      if (label.units == 0) {
        return(rep(FALSE, length(w)))
      }
      positive & w >= sort(w[positive], decreasing = TRUE)[label.units]
    }),
    use.names = FALSE
  )[order(order(data$panel))]

  if (type == "scatter") {
    if (weighted.only) {
      data <- data[data$weight > 0, ]
    }
    data$name <- ifelse(data$named, data$unit, "")
    # The effect's label sits midway along its line. Points spread under it keep
    # the controls' labels off it.
    ranges <- tapply(
      data$weight,
      data$panel,
      function(w) {
        range(c(0, w))
      },
      simplify = FALSE
    )
    effects$x <- vapply(
      as.character(effects$panel),
      function(f) {
        mean(ranges[[f]])
      },
      0
    )
    guard <- do.call(
      rbind,
      lapply(seq_len(nrow(effects)), function(k) {
        width <- diff(ranges[[as.character(effects$panel[k])]])
        data.frame(
          panel = effects$panel[k],
          weight = effects$x[k] + seq(-.07, .07, length.out = 7) * width,
          own = effects$effect[k],
          name = "",
          size = 5
        )
      })
    )
    repel <- rbind(
      data.frame(
        panel = data$panel,
        weight = data$weight,
        own = data$own,
        name = data$name,
        size = 2,
        shade = ifelse(data$weight > 0, "black", "grey50")
      ),
      cbind(guard, shade = "black")
    )
    weighted <- data[data$weight > 0, ]
    labels <- if (requireNamespace("ggrepel", quietly = TRUE)) {
      list(
        ggrepel::geom_text_repel(
          ggplot2::aes(
            x = weight,
            y = own,
            label = name,
            point.size = size,
            colour = shade
          ),
          data = repel,
          size = 3,
          max.overlaps = Inf,
          min.segment.length = 0,
          segment.colour = "grey60",
          segment.size = .25,
          box.padding = .25,
          point.padding = .1,
          seed = 1
        ),
        ggplot2::scale_colour_identity()
      )
    } else {
      list(ggplot2::geom_text(
        ggplot2::aes(x = weight, y = own, label = name),
        data = data,
        size = 3,
        vjust = -.8,
        check_overlap = TRUE,
        colour = "grey20"
      ))
    }
    # No border round the effect's label. ggplot2 4.0.0 renamed the argument,
    # from label.size to linewidth, which it takes through ..., so the version
    # decides which is used.
    border <- if (utils::packageVersion("ggplot2") >= "4.0.0") {
      list(linewidth = 0)
    } else {
      list(label.size = 0)
    }
    x.axis <- machinery$axis.functions("x", FALSE, NULL)
    y.axis <- machinery$axis.functions("y", FALSE, NULL)
    p <- ggplot2::ggplot(data) +
      ggplot2::geom_hline(
        ggplot2::aes(yintercept = effect),
        data = effects,
        colour = "red",
        linetype = "dashed"
      ) +
      ggplot2::geom_point(
        ggplot2::aes(x = weight, y = own),
        data = weighted,
        size = 2
      ) +
      ggplot2::geom_point(
        ggplot2::aes(x = weight, y = own),
        data = data[data$weight == 0, ],
        shape = 1,
        size = 2,
        colour = "grey55"
      ) +
      labels +
      do.call(
        ggplot2::geom_label,
        c(
          list(
            mapping = ggplot2::aes(x = x, y = effect, label = label),
            data = effects,
            colour = "red",
            fill = "white",
            size = 3.4,
            inherit.aes = FALSE
          ),
          border
        )
      ) +
      ggplot2::scale_x_continuous(
        limits = x.axis$limits,
        breaks = x.axis$breaks,
        labels = x.axis$labels,
        expand = ggplot2::expansion(mult = padding),
        oob = scales::oob_keep
      ) +
      ggplot2::scale_y_continuous(
        limits = y.axis$limits,
        breaks = y.axis$breaks,
        labels = y.axis$labels,
        expand = ggplot2::expansion(mult = padding),
        oob = scales::oob_keep
      ) +
      ggplot2::labs(x = "Unit weight", y = "Estimate with this control alone") +
      facet +
      look +
      ggplot2::coord_cartesian(clip = "off")
    return(add_panel_frames(p))
  }

  # The bars, for controls with weight, ordered within each panel from the
  # largest weight down. Their axis is padded by 8% rather than 2%, to make room
  # for the weights at the bars' ends.
  machinery <- tick_machinery(.08)
  bars <- data[data$weight > 0 & data$named, ]
  bars$id <- paste(bars$panel, bars$unit, sep = "\r")
  bars$id <- factor(bars$id, levels = bars$id[order(bars$panel, bars$weight)])
  bars$direction <- factor(
    ifelse(bars$own > bars$effect, "above the ATT", "below the ATT"),
    levels = c("above the ATT", "below the ATT")
  )
  bars$text <- formatC(bars$weight, format = "f", digits = 3)
  bars$right <- bars$own > bars$effect
  # Each panel's ticks are chosen from its controls' estimates and the effect,
  # and the effect gets a tick of its own, unlabeled, as its label sits above
  # its line.
  limits <- lapply(split(bars, bars$panel, drop = TRUE), function(b) {
    core <- range(c(b$own, b$effect[1]))
    base::range(c(machinery$choose.ticks(core[1], core[2])$at, core))
  })
  x.axis <- machinery$axis.functions("x", FALSE, NULL)
  effect.of <- function(range) {
    width <- diff(range) / (1 + 2 * .08)
    unpadded <- range + c(1, -1) * .08 * width
    for (f in names(limits)) {
      if (isTRUE(all.equal(limits[[f]], unpadded, tolerance = 1e-8))) {
        return(effects$effect[effects$panel == f])
      }
    }
    NULL
  }
  p <- ggplot2::ggplot(bars) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = effect,
        xend = own,
        y = id,
        yend = id,
        colour = direction
      ),
      linewidth = 4.2
    ) +
    ggplot2::geom_vline(
      ggplot2::aes(xintercept = effect),
      data = effects,
      colour = "red",
      linetype = "dashed"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = own,
        y = id,
        label = text,
        hjust = ifelse(right, -.15, 1.15)
      ),
      size = 3,
      colour = "grey30"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(x = effect, y = Inf, label = label),
      data = effects,
      colour = "red",
      vjust = -.6,
      size = 3.4
    ) +
    ggplot2::scale_colour_manual(
      values = c("above the ATT" = "grey35", "below the ATT" = "grey70"),
      drop = FALSE
    ) +
    ggplot2::scale_x_continuous(
      limits = function(range) {
        x.axis$limits(range)
      },
      breaks = function(range) {
        c(x.axis$breaks(range), effect.of(range))
      },
      labels = function(breaks) {
        effect <- vapply(
          breaks,
          function(b) {
            isTRUE(any(abs(b - effects$effect) < 1e-9))
          },
          TRUE
        )
        out <- rep("", length(breaks))
        out[!effect] <- format(breaks[!effect], trim = TRUE)
        out
      },
      expand = ggplot2::expansion(mult = .08),
      oob = scales::oob_keep
    ) +
    ggplot2::scale_y_discrete(labels = function(id) {
      sub("^.*\r", "", id)
    }) +
    ggplot2::labs(x = NULL, y = NULL) +
    facet +
    look +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(16, 5.5, 5.5, 5.5)
    ) +
    ggplot2::coord_cartesian(clip = "off")
  add_panel_frames(p)
}

# Open frames for the panels of a plot, one or several, each axis line running
# from the lowest to the highest of the panel's data and ticks, or over the
# positions of a discrete axis.
add_panel_frames <- function(p) {
  built <- ggplot2::ggplot_build(p)
  layout <- built$layout$layout
  finite <- function(v) {
    suppressWarnings(as.numeric(v))[is.finite(suppressWarnings(as.numeric(v)))]
  }
  frames <- list()
  for (i in seq_len(nrow(layout))) {
    in.panel <- lapply(built$data, function(d) {
      d[d$PANEL == layout$PANEL[i], , drop = FALSE]
    })
    params <- built$layout$panel_params[[i]]
    spans <- lapply(c("x", "y"), function(axis) {
      values <- finite(unlist(lapply(in.panel, function(d) {
        if (axis == "x") {
          c(d$x, d$xend, d$xintercept)
        } else {
          c(d$y, d$ymin, d$ymax, d$yend, d$yintercept)
        }
      })))
      discrete <- isTRUE(tryCatch(
        params[[axis]]$is_discrete(),
        error = function(e) FALSE
      ))
      if (discrete) {
        range(values)
      } else {
        range(c(values, finite(params[[axis]]$breaks)))
      }
    })
    panel <- if ("panel" %in% names(layout)) {
      layout$panel[i]
    } else {
      NULL
    }
    add <- function(x, y, group) {
      d <- data.frame(x = x, y = y, group = paste(group, i))
      if (!is.null(panel)) {
        d$panel <- panel
      }
      d
    }
    frames[[length(frames) + 1]] <- add(spans[[1]], -Inf, "x")
    frames[[length(frames) + 1]] <- add(-Inf, spans[[2]], "y")
  }
  p +
    ggplot2::theme(
      axis.line.x.bottom = ggplot2::element_blank(),
      axis.line.y.left = ggplot2::element_blank()
    ) +
    ggplot2::geom_path(
      ggplot2::aes(x = x, y = y, group = group),
      data = do.call(rbind, frames),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = .4
    )
}


#' Plot a synthdid object
#' @param x The object to plot
#' @param ... Additional arguments (currently ignored).
#' @method plot sdid_estimate
#' @export
#' @keywords internal
plot.sdid_estimate <- function(x, ...) {
  require_ggplot2()
  sdid_plot(x, ...)
}
