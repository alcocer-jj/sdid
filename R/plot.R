# ggplot2::aes() column names; declared to quiet R CMD check.
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

# Stop if ggplot2 is missing. Called via namespace, never attached.
require_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Plotting requires the package `ggplot2`. Install it to use this function."
    )
  }
}

# Draws a separate time-weight axis to the right of each stacked panel. A single
# secondary axis won't work because each panel's weights have their own scale.
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

# Cohort adoption-period label under its panel, ticked on the time-axis panel.
# That label shifts just enough to clear the time-axis labels in `neighbors`,
# or drops a line if no room. Resolved at draw time: overlap depends on size.
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

# Adoption-period label at `at` (share of panel width), optional tick, and
# the time-axis labels it must avoid.
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
  # Already labeled by a time-axis tick; skip.
  if (any(abs(neighbors - x$at) < 1e-9)) {
    return(grid::setChildren(x, children))
  }
  # Positions in pt from panel's left edge; pair spacing = half widths + 5.5 pt.
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
      # Room between neighbors: move minimally.
      position <- min(max(position, lowest), highest)
    } else {
      # No room: drop to a line below, under the tick.
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

# Frame for column-stacked panels: outcome axis lines (as add_range_frame),
# time axis on bottom panel, per-panel weight axis, staggered cohort labels
# (bottom one shifted to clear time-axis labels).
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
  # As the single-panel axis, but labels >= 12% of panel height apart (not 5%):
  # stacked panels are shorter.
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
    # Time-axis labels the bottom cohort label must clear.
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

# Shared tick logic: loose ticks as Stata picks them. `padding` = proportional
# axis padding. Returns choose.ticks() and axis.functions() (limits, breaks,
# labels; memoized so they agree).
tick_machinery <- function(padding = .02) {
  # Steps: 1, 5, 2, 2.5, 4, 3 x 10^k. Whole-number (e.g. year) axes drop 2.5 to
  # avoid ticks like 2002.5, and fall back to integer steps if still fractional.
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
    # User breaks used as-is; axis spans them and the data.
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
  # Memoize tick choices by range: ggplot2 asks for limits, breaks, labels
  # separately, passes padded ranges, and later sees range-frame lines as data.
  # Match by range with rounding tolerance; keyed per axis (x, y never share).
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
  # Limits/breaks/labels functions for an axis, with any user breaks.
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

# Single-panel range frame: each axis line spans min..max of data and breaks,
# not the whole panel. Weight axis keeps its own line.
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


#' Plot treated and synthetic control trajectories
#'
#' Plots the treated and synthetic control trajectories of an `sdid()` result,
#' along with the time weights, which show which pre-treatment periods were
#' most informative for predicting post-treatment outcomes among the controls.
#'
#' Two plotting styles are available: "minimal" and "effects". The default,
#' `style = "minimal"`, plots both trajectories in black, with the treated line
#' solid and the synthetic control dashed. The time weights are plotted as a
#' gray band with their scale on a right-hand axis, and a red dashed line marks
#' treatment. `style = "effects"` adds a diff-in-diff diagram of the treatment
#' effect and fades the trajectory lines.
#'
#' A staggered `sdid()` result is drawn as one panel per cohort, stacked in a
#' column and labeled by adoption period. A list of results is drawn the same
#' way, one panel each. To change the axis ticks, use `x.breaks` and `y.breaks`
#' rather than adding your own scale, which would drop the time-weight axis.
#'
#' Requires ggplot2.
#'
#' @param estimates The result of an sdid object or a list of sdid estimates.
#'        For a staggered sdid result, each adoption period is plotted in
#'        a stacked panel.
#' @param treated.name A string conveying the name of the treated curve that
#'        appears in the legend. Defaults to "treated".
#' @param control.name A string conveying the name of the control curve that
#'        appears in the legend. Defaults to 'synthetic control'.
#' @param spaghetti.units A list of units to plot individually alongside the
#'        treatment and control curves. Defaults to the empty list.
#' @param spaghetti.matrices A list of matrices, one per estimate, of extra
#'        trajectories to plot. Each row is a trajectory, one column per
#'        period, and row names are used as labels.
#' @param lambda.comparable Logical argument where if TRUE, the lambda
#'        weight axis is scaled to range from 0 to 1. Useful for when comparing
#'        across stacked panels. Defaults to FALSE.
#' @param overlay Between 0 and 1. Shifts the synthetic control toward the
#'        treated line by this fraction of their pre-treatment gap, which makes
#'        it easier to see whether they move in parallel. Defaults to 0 (no
#'        shift). A vector of values can be animated with plotly.
#' @param lambda.plot.scale Height of the time-weight band: the range of the
#'        trajectories divided by this number. Defaults to 3.
#' @param effect.curvature Curvature of the effect arrow in the effects style.
#' @param line.width Width of the trajectory lines.
#' @param point.size Size of the points in the effects style.
#' @param onset.alpha Opacity of the red dashed line that marks treatment.
#' @param spaghetti.line.width,spaghetti.label.size,spaghetti.line.alpha,spaghetti.label.alpha
#'        Width, label size and opacity of the spaghetti trajectories and their
#'        labels.
#' @param alpha.multiplier One number per estimate that scales its opacity,
#'        to emphasize some estimates over others.
#' @param style The plot's style, "minimal" (the default) or "effects". See
#'        Details.
#' @param x.breaks,y.breaks Optional tick positions for the time and outcome
#'        axes (numbers, or dates for a date axis). The axes are extended to
#'        the outermost ticks if needed. With stacked panels, `y.breaks`
#'        applies to every panel.
#' @param lambda.axis Logical argument where if TRUE, each panel with a time
#'        weights' band shows the band's scale on a right-hand axis. Defaults
#'        to TRUE.
#' @return A ggplot object, which can be extended with `+`.
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

    # use estimate-specific overlay attribute if present; SC estimates/overlays
    # are plotted differently
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

    # Treatment line midway between last pre-period and adoption period
    # (recorded with it).
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

  # add frame=frame to aes only with multiple frames: ggplotly animates on it,
  # ggplot2 warns 'ignoring unknown aesthetic'. returns
  # geom(aes, data=data, ...) where aes=base.aes if there's a single frame
  #					    aes=(base.aes, frame=frame) if there's more than one
  with.frame <- function(geom, base.aes, data, ...) {
    new.aes <- if (multiple.frames) {
      modifyList(base.aes, ggplot2::aes(frame = frame))
    } else {
      base.aes
    }
    do.call(geom, c(list(new.aes, data = data), list(...)))
  }

  # Alpha values are literal opacities (identity scale). Effects style:
  # half-opacity trajectories plus DiD diagram.
  ribbon.alpha <- .34
  trajectory.opacity <- if (style == "effects") {
    .5
  } else {
    1
  }
  ribbons <- no.sc(conc$ribbons)
  # Numeric axes as Stata labels them: ~5 loose ticks (Talbot, Lin & Hanrahan
  # 2010, ggplot2's default) covering the data, range widened to outer ticks.
  # 2% padding keeps content close to the axis lines.
  padding <- .02
  machinery <- tick_machinery(padding)
  choose.ticks <- machinery$choose.ticks
  axis.functions <- machinery$axis.functions
  y.axis <- axis.functions("y", FALSE, y.breaks)
  # Whole-number time axis iff all periods are integers.
  x.axis <- axis.functions(
    "x",
    all(abs(conc$lines$x - round(conc$lines$x)) < 1e-9),
    x.breaks
  )
  # Snap each band's baseline down to the nearest tick so zero weight lines up
  # with a labeled value. Moving changes the range, hence ticks, so iterate.
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
  # Keep positive-weight stretches plus a flanking zero so each rises from the
  # baseline; each stretch is its own group.
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
    # No band, e.g. SC (all time weights zero).
    ribbons <- no.sc(conc$ribbons)[0, ]
    ribbons$stretch <- character(0)
  }
  # Outline each positive stretch fully, including beneath.
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
    # DiD diagram in black: solid treated/SC trends, dotted counterfactual
    # parallel to SC, red effect arrow.
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
    # Stack panels in a column, in input order.
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
            # Dates are days since 1970-01-01; labels convert back.
            # R < 4.3.0 needs origin.
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
        # Keep ggplot2's date ticks: nice day counts aren't nice dates.
        ggplot2::scale_x_continuous(
          limits = x.axis$limits,
          breaks = x.axis$breaks,
          labels = x.axis$labels,
          expand = ggplot2::expansion(mult = padding)
        )
      }
    )

  # Right axis for time weights: ribbon height h above band bottom = weight
  # h / band.height * band.top. Single band/panel only.
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
      # Fewer ticks when weights fill little of the panel (labels ~5% apart).
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
  # Same vertical padding without a weight axis.
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
  # Single panel: open frame; outcome/time axes drawn as a range frame below,
  # weight axis stops at outer ticks. Stacked: same frame, time axis on bottom
  # only so gridlines run through; cohorts labeled by adoption period, others by
  # title. Side-by-side panels keep boxes.
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
    # Right margin for weight axes; bottom for cohort labels pushed below the
    # time-axis labels.
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
    # Only with a weight axis: otherwise ggplot2 >= 3.5 copies the left axis
    # to the right.
    if (
      utils::packageVersion("ggplot2") >= "3.5.0" && !is.null(weight.title.y)
    ) {
      p <- p + ggplot2::guides(y.sec = ggplot2::guide_axis(cap = "both"))
    }
  }
  # Weight axis title: vertical, centered on the axis span (0 to top tick), not
  # the panel. Placed in data coords past the ticks, clipping off. Widths left
  # to grid at draw time; measuring here would open a graphics device.
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

#' Plot an estimate with its in-time placebo
#'
#' Plots the treated and synthetic control trajectories of an `sdid()` result,
#' with an in-time placebo stacked below. Requires ggplot2.
#' @param estimate The result of `sdid()` with a single adoption period.
#' @param treated.fraction The share of the periods before treatment that
#'        the placebo treats, as in `sdid_placebo()`. Defaults to NULL.
#' @param style The plot's style, "minimal" (the default) or "effects", as
#'        in `sdid_plot()`.
#' @param placebo.time The period at which the placebo treatment begins, as
#'        in `sdid_placebo()`. Defaults to "auto".
#' @param x.breaks,y.breaks Optional ticks for the time and outcome axes, as
#'        in `sdid_plot()`. The outcome ticks apply to both panels.
#' @return A ggplot object, which can be extended with `+`.
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
#' estimate, with the overall effect as a red line. Controls with zero weight
#' are drawn as open circles. `type = "bars"` draws a bar from the effect to
#' each weighted control's estimate, largest weight first; bars to the right
#' pull the estimate up, bars to the left pull it down.
#'
#' Requires ggplot2. If ggrepel is installed, it is used to keep labels from
#' overlapping.
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
#' @param label.units Which controls to label: "all", "weighted" (positive
#'        weight only), or a number n for the n largest weights. Defaults to
#'        "all" for one panel and 10 for several. The bars only show labeled
#'        controls.
#' @param weighted.only TRUE to leave controls with no weight out of the
#'        scatter. Defaults to FALSE, drawing them as open circles.
#' @return A ggplot object, which can be extended with `+`.
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

  # Per-control DiD (treated vs that control alone); weighted mean = estimate.
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

  # Controls picked by label.units per panel (named in scatter, drawn as bars).
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
    # Effect label midway along its line; points spread under it keep control
    # labels clear.
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
    # No label border. ggplot2 4.0.0 renamed label.size -> linewidth (via ...),
    # so branch on version.
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

  # Bars for weighted controls, largest first per panel. 8% padding (not 2%)
  # leaves room for weight labels.
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
  # Ticks from controls' estimates and the effect; effect gets its own unlabeled
  # tick (its label sits above its line).
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

# Open frames for one or more panels; axis lines span data and ticks, or a
# discrete axis's positions.
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
